import sys
from os.path import basename, splitext, isfile, exists
from os import makedirs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
pd.options.display.mpl_style = 'default'
font = {'family': 'sans-serif', 'size': 7}
matplotlib.rc('font', **font)
from concurrent import futures
from peakachulib.library import Library
from peakachulib.wiggle import WiggleWriter
from peakachulib.deseq2 import RunDESeq2
from time import time
from collections import OrderedDict
from subprocess import Popen, PIPE


class PredefinedPeakApproach(object):
    '''
    This class is used for peak detection via predefining peaks based on shape
    and subsequent comparison to a control
    '''
    def __init__(self, replicon_dict, max_proc, padj_threshold, fc_cutoff,
                 output_folder):
        self._lib_dict = OrderedDict()
        self._replicon_dict = replicon_dict  # own copy of replicon_dict
        self._max_proc = max_proc
        self._padj_threshold = padj_threshold
        self._fc_cutoff = fc_cutoff
        self._output_folder = output_folder + "/predefined_peak_approach"
        if not exists(self._output_folder):
            makedirs(self._output_folder)
        
    def init_libraries(self, paired_end, max_insert_size, ctr_libs,
                       exp_libs):
        self._ctr_lib_list = [splitext(basename(lib_file))[0]
                              for lib_file in ctr_libs]
        self._exp_lib_list = [splitext(basename(lib_file))[0]
                              for lib_file in exp_libs]
        # add libs to lib_dict
        for lib_file in exp_libs + ctr_libs:
            if not isfile(lib_file):
                sys.stderr.write("ERROR: The library file %s does not exist.\n"
                                 % lib_file)
                sys.exit(1)
            self._lib_dict[splitext(basename(lib_file))[0]] = Library(
                paired_end, max_insert_size, lib_file, self._replicon_dict)
        self._lib_names_list = list(self._lib_dict.keys())
        print("The following libraries were initialized:\n"
              "# Experiment libraries\n{0}\n"
              "# Control libraries\n{1}".format(
                  '\n'.join(self._exp_lib_list),
                  '\n'.join(self._ctr_lib_list)))
              
    def generate_combined_bed_file(self):
        # execute read conversion in parallel
        print("** Converting reads to bed format for %s libraries..." % len(
            self._exp_lib_list), flush=True)
        exp_lib_dict = {lib_name: self._lib_dict[lib_name] for lib_name in
                        self._exp_lib_list}
        t_start = time()
        with futures.ProcessPoolExecutor(
                max_workers=self._max_proc) as executor:
            future_to_lib_name = {
                executor.submit(lib.merge_reads):
                lib.lib_name for lib in exp_lib_dict.values()}
        for future in futures.as_completed(future_to_lib_name):
            lib_name = future_to_lib_name[future]
            try:
                self._lib_dict[lib_name].replicon_dict = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (lib_name, exc),
                      flush=True)
        output_df = pd.DataFrame()
        for replicon in sorted(self._replicon_dict):
            self._replicon_dict[replicon]["reads"] = pd.Series()
            for lib_name, lib in exp_lib_dict.items():
                self._replicon_dict[replicon]["reads"] = self._replicon_dict[
                    replicon]["reads"].add(lib.replicon_dict[replicon][
                        "reads"], fill_value=0)
            self._replicon_dict[replicon]["reads"] = self._replicon_dict[
                replicon]["reads"].reset_index(name="count")
            split_index = pd.DataFrame(list(self._replicon_dict[replicon][
                "reads"]["index"].str.split(',')), columns=[
                    "start", "end", "strand"])
            split_index = split_index.convert_objects(convert_numeric=True)
            del self._replicon_dict[replicon]["reads"]["index"]
            self._replicon_dict[replicon]["reads"] = split_index.join(
                self._replicon_dict[replicon]["reads"]).sort(
                    ["strand", "start", "end"], ascending=[False, True, True])
            self._replicon_dict[replicon]["reads"]["replicon"] = replicon
            output_df = output_df.append(
                self._replicon_dict[replicon]["reads"], ignore_index=True)
        output_df["tag_id"] = (output_df.index + 1).map('tag_{:.0f}'.format)
        output_df = output_df.loc[:,
                                  ["replicon",
                                   "start",
                                   "end",
                                   "tag_id",
                                   "count",
                                   "strand"]]
        output_df.to_csv(
            "%s/sorted_reads_for_blockbuster.bed" % (self._output_folder),
            sep='\t', header=False, index=False, encoding='utf-8')
        t_end = time()
        print("Reads converted to bed format in %s seconds.\n" % (
            t_end-t_start), flush=True)
        
    def run_blockbuster(self):
        p = Popen(["blockbuster.x", "-minBlockHeight", "10", "-print",
                   "1", "-distance", "1",
                   "{}/sorted_reads_for_blockbuster.bed".format(
                       self._output_folder)], stdout=PIPE, stderr=PIPE,
                  universal_newlines=True)
        self._blockbuster_output, err = p.communicate()
        print("blockbuster exited with status {}".format(p.returncode))
        if not p.returncode == 0:
            print(err)
            sys.exit(1)
        with open("%s/blockbuster.txt" % (self._output_folder),
                  'w') as blockbuster_fh:
            blockbuster_fh.write(self._blockbuster_output)
        
    def generate_peaks_from_blockbuster(self, min_cluster_expr_frac,
                                        min_block_overlap,
                                        min_max_block_expr_frac):
        for replicon in self._replicon_dict:
            self._replicon_dict[replicon]["peak_df"] = pd.DataFrame()
        cluster = {}
        for line in self._blockbuster_output.rstrip().split('\n'):
            if line.startswith('>'):
                if cluster:
                    self._call_cluster_peaks(cluster, min_cluster_expr_frac,
                                             min_block_overlap,
                                             min_max_block_expr_frac)
                    cluster = {}
                cluster["header"] = line
                cluster["blocks"] = []
            else:
                cluster["blocks"].append(line)
        if cluster:
            self._call_cluster_peaks(cluster, min_cluster_expr_frac,
                                     min_block_overlap,
                                     min_max_block_expr_frac)
            
    def _call_cluster_peaks(self, cluster, min_cluster_expr_frac,
                            min_block_overlap, min_max_block_expr_frac):
        cluster_entries = cluster["header"].strip().split('\t')
        cluster_expr = float(cluster_entries[5])
        cluster_strand = cluster_entries[4]
        cluster_replicon = cluster_entries[1]
        peak_df = pd.DataFrame()
        
        if len(cluster["blocks"]) == 1:
            block_entries = cluster["blocks"][0].strip().split('\t')
            peak_start = int(block_entries[2]) + 1
            peak_end = int(block_entries[3])
            peak_df = peak_df.append(pd.Series([peak_start, peak_end], index=[
                "peak_start", "peak_end"]), ignore_index=True)
        else:
            blocks = [block.strip().split('\t') for block in cluster["blocks"]]
            block_df = pd.DataFrame(blocks, columns=[
                "blockNb", "blockChrom",
                "blockStart", "blockEnd", "blockStrand", "blockExpression",
                "readCount"])
            block_df = block_df.convert_objects(convert_numeric=True)
            peak_df = self._split_cluster_peaks(block_df, cluster_expr,
                                                peak_df, min_cluster_expr_frac,
                                                min_block_overlap,
                                                min_max_block_expr_frac)
        peak_df = peak_df.astype(np.int64)
        peak_df["peak_strand"] = cluster_strand
        self._replicon_dict[cluster_replicon]["peak_df"] = self._replicon_dict[
            cluster_replicon]["peak_df"].append(peak_df, ignore_index=True)
            
    def _split_cluster_peaks(self, block_df, cluster_expr, peak_df,
                             min_cluster_expr_frac, min_block_overlap,
                             min_max_block_expr_frac):
        if block_df.empty:
            return peak_df
        max_block_ix = block_df["blockExpression"].idxmax()
        max_block_expr = block_df.loc[max_block_ix, "blockExpression"]
        if max_block_expr/cluster_expr < min_cluster_expr_frac:
            return peak_df
        min_overlap = round(
            (block_df.loc[max_block_ix, "blockEnd"] -
                block_df.loc[max_block_ix, "blockStart"]) * min_block_overlap)
        overlaps_with_max_block = (block_df.loc[:, "blockEnd"].apply(
            min, args=(block_df.loc[
                max_block_ix, "blockEnd"],)) - block_df.loc[
                    :, "blockStart"].apply(
                        max, args=(block_df.loc[
                            max_block_ix, "blockStart"],))).apply(
                                max, args=(0,))
        peak_blocks = block_df.loc[overlaps_with_max_block >= min_overlap, :]
        peak_blocks = peak_blocks.loc[
            (peak_blocks["blockExpression"] /
                max_block_expr) >= min_max_block_expr_frac, :]
        peak_start = peak_blocks["blockStart"].min()
        peak_end = peak_blocks["blockEnd"].max()
        overlaps_with_peak = (block_df.loc[:, "blockEnd"].apply(min, args=(
            peak_end,)) - block_df.loc[:, "blockStart"].apply(max, args=(
                peak_start,))).apply(max, args=(0,))
        next_block_df = block_df.loc[overlaps_with_peak == 0, :].reset_index(
            drop=True)
        peak_df = peak_df.append(pd.Series([peak_start + 1, peak_end], index=[
            "peak_start", "peak_end"]), ignore_index=True)
        return self._split_cluster_peaks(next_block_df, cluster_expr, peak_df,
                                         min_cluster_expr_frac,
                                         min_block_overlap,
                                         min_max_block_expr_frac)
        
    def calculate_peak_expression(self):
        for lib in self._lib_dict.values():
            for replicon in self._replicon_dict:
                lib.replicon_dict[replicon]["peak_df"] = self._replicon_dict[
                    replicon]["peak_df"]
        self._generate_peak_counts()
        self._peak_df = pd.DataFrame()
        for replicon in sorted(self._replicon_dict):
            self._replicon_dict[replicon]["peak_df"]["replicon"] = replicon
            for lib_name, lib in self._lib_dict.items():
                self._replicon_dict[replicon][
                    "peak_df"][lib_name] = lib.replicon_dict[
                        replicon]["peak_counts"]
            # add pseudocounts
            # self._replicon_dict[
            #    replicon]["peak_df"].loc[:, self._lib_names_list] += 1.0
            self._peak_df = self._peak_df.append(self._replicon_dict[replicon][
                "peak_df"])
        
    def _generate_peak_counts(self):
        # execute read counting in parallel
        print("** Peak read counting started for %s libraries..." % len(
            self._lib_dict), flush=True)
        t_start = time()
        with futures.ProcessPoolExecutor(
                max_workers=self._max_proc) as executor:
            future_to_lib_name = {
                executor.submit(lib.count_reads_for_peaks):
                lib.lib_name for lib in self._lib_dict.values()}
        for future in futures.as_completed(future_to_lib_name):
            lib_name = future_to_lib_name[future]
            try:
                self._lib_dict[lib_name].replicon_dict = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (lib_name, exc),
                      flush=True)
        t_end = time()
        print("Peak read counting finished in %s seconds." % (t_end-t_start),
              flush=True)
    
    def run_deseq2_analysis(self, size_factors):
        count_df = self._peak_df.loc[:, self._exp_lib_list +
                                     self._ctr_lib_list]
        run_deseq2 = RunDESeq2(
            count_df, self._exp_lib_list, self._ctr_lib_list, size_factors)
        result_df, self._size_factors = run_deseq2.run_deseq2()
        # normalize counts
        self._peak_df[self._lib_names_list] = self._peak_df[
            self._lib_names_list].div(self._size_factors, axis='columns')
        # append DESeq2 output
        self._peak_df = pd.concat([self._peak_df, result_df], axis=1)
        # write initial peaks
        peak_columns = (["replicon",
                         "peak_start",
                         "peak_end",
                         "peak_strand"] +
                        [lib_name for lib_name in self._lib_dict] +
                        ["baseMean",
                         "log2FoldChange",
                         "lfcSE",
                         "stat",
                         "pvalue",
                         "padj"])
        self._peak_df.loc[:, peak_columns].to_csv(
            "%s/initial_peaks.csv" % (self._output_folder),
            sep='\t', index=False, encoding='utf-8')
        # filter peaks
        print("* Filtering peaks...", flush=True)
        sig_peak_df = self._filter_peaks(self._peak_df)
        unsig_peak_df = self._peak_df[~self._peak_df.index.isin(
            sig_peak_df.index)]
        self._plot_initial_peaks(unsig_peak_df.baseMean,
                                 np.power(2.0, unsig_peak_df.log2FoldChange),
                                 sig_peak_df.baseMean,
                                 np.power(2.0, sig_peak_df.log2FoldChange))
        self._peak_df = sig_peak_df
    
    def _filter_peaks(self, df):
        print("Removing peaks based on minimum fold change "
              "from DataFrame with %s rows..." % len(df), flush=True)
        t_start = time()
        log2_fc_cutoff = np.log2(self._fc_cutoff)
        df = df.query('log2FoldChange >= @log2_fc_cutoff')
        t_end = time()
        print("Removal took %s seconds. DataFrame contains now %s rows." % (
            (t_end-t_start), len(df)), flush=True)
        print("Removing peaks based on padj from DataFrame with %s rows..."
              % len(df), flush=True)
        t_start = time()
        df = df.query('padj < @self._padj_threshold')
        t_end = time()
        print("Removal took %s seconds. DataFrame contains now %s rows." % (
            (t_end-t_start), len(df)), flush=True)
        return df
        
    def write_output(self):
        peak_columns = (["replicon",
                         "peak_id",
                         "peak_start",
                         "peak_end",
                         "peak_strand"] +
                        [lib_name for lib_name in self._lib_dict] +
                        ["baseMean",
                         "log2FoldChange",
                         "lfcSE",
                         "stat",
                         "pvalue",
                         "padj"])
        feature_columns = ["feature_type",
                           "feature_start",
                           "feature_end",
                           "feature_strand",
                           "feature_locus_tag",
                           "feature_name",
                           "subfeature_type",
                           "feature_product",
                           "overlap_length"]
        for replicon in self._replicon_dict:
            self._replicon_dict[replicon]["peak_df"] = self._peak_df.loc[
                self._peak_df["replicon"] == replicon, :].copy()
            if self._replicon_dict[replicon]["peak_df"].empty:
                continue
            output_df = pd.DataFrame()
            self._replicon_dict[replicon]["peak_df"].sort(
                ["replicon", "peak_start"], inplace=True)
            self._replicon_dict[replicon]["peak_df"].reset_index(
                drop=True, inplace=True)
            self._replicon_dict[replicon]["peak_df"].loc[:, "peak_id"] = (
                self._replicon_dict[replicon]["peak_df"].index + 1)
            for peak in self._replicon_dict[replicon]["peak_df"].to_dict(
                    'records'):
                overlapping_features = self._find_overlapping_features(peak)
                for match in overlapping_features:
                    output_df = output_df.append(pd.Series(peak).append(
                        pd.Series(match)), ignore_index=True)
            output_df = output_df.loc[:, peak_columns + feature_columns]
            output_df.to_csv(
                "%s/peaks_%s.csv" % (self._output_folder, replicon),
                sep='\t', index=False, encoding='utf-8')

            self._write_gff_file(replicon, self._replicon_dict[replicon]
                                 ["peak_df"])
        
    def _find_overlapping_features(self, peak):
        overlapping_features = []
        for feature in self._replicon_dict[peak["replicon"]]["features"]:
            if not peak["peak_strand"] == feature["strand"]:
                continue
            overlap = self._get_overlap(peak["peak_start"], peak["peak_end"],
                                        feature["start"], feature["end"])
            if not overlap > 0:
                continue
            overlapping_features.append({
                "feature_type": feature["type"],
                "feature_start": feature["start"],
                "feature_end": feature["end"],
                "feature_strand": feature["strand"],
                "feature_locus_tag": feature["locus_tag"],
                "feature_name": feature["Name"],
                "subfeature_type": feature["subfeature_type"] if (
                    "subfeature_type" in feature) else None,
                "feature_product": feature["product"],
                "overlap_length": overlap})
        if not overlapping_features:
            overlapping_features.append({
                "feature_type": "intergenic",
                "feature_start": None,
                "feature_end": None,
                "feature_strand": None,
                "feature_locus_tag": None,
                "feature_name": None,
                "subfeature_type": None,
                "feature_product": None,
                "overlap_length": None})
        return overlapping_features

    def generate_normalized_wiggle_files(self):
        wiggle_folder = "%s/normalized_coverage" % (self._output_folder)
        if not exists(wiggle_folder):
            makedirs(wiggle_folder)
        
        # Generate coverage files in parallel
        print("** Generating normalized coverage files for %s libraries..." %
              (len(self._lib_dict)), flush=True)
        t_start = time()
        with futures.ProcessPoolExecutor(
                max_workers=self._max_proc) as executor:
            future_to_lib_name = {
                executor.submit(
                    self._generate_normalized_wiggle_file_for_lib,
                    lib, size_factor, wiggle_folder):
                lib.lib_name for lib, size_factor in zip(
                    self._lib_dict.values(), self._size_factors)}
        for future in futures.as_completed(future_to_lib_name):
            lib_name = future_to_lib_name[future]
            print("* Coverage files for library %s generated." % lib_name)
        t_end = time()
        print("Coverage file generation finished in %s seconds." % (
            t_end-t_start), flush=True)
    
    def _generate_normalized_wiggle_file_for_lib(self, lib, size_factor,
                                                 wiggle_folder):
        """Perform the coverage calculation for a given library."""
        strand_dict = {"+": "forward", "-": "reverse"}
        wiggle_writers = dict([(strand, WiggleWriter(
            "%s_%s" % (lib.lib_name, strand),
            open("%s/%s_div_by_%s_%s.wig" % (
                wiggle_folder, lib.lib_name, size_factor, strand),
                "w"))) for strand in strand_dict.values()])
        for replicon in self._replicon_dict:
            for strand in strand_dict:
                if strand == "-":
                    factor = -1.0/size_factor
                else:
                    factor = 1.0/size_factor
                try:
                    wiggle_writers[strand_dict[
                        strand]].write_replicons_coverages(
                            replicon, lib.replicon_dict[replicon]["coverages"][
                                strand], factor=factor)
                except Exception as exc:
                    print("Library %s, replicon %s, %s strand generated an "
                          "exception during coverage file generation: %s" %
                          (lib.lib_name, replicon, strand, exc), flush=True)
        for strand in strand_dict:
            wiggle_writers[strand].close_file()
    
    def _write_gff_file(self, replicon, df):
        with open("%s/peaks_%s.gff" % (self._output_folder, replicon),
                  'w') as out_gff_fh:
            out_gff_fh.write("##gff-version 3\n"
                             "#!gff-spec-version 1.20\n"
                             "##sequence-region %s %s %s\n"
                             "%s%s"
                             "###\n" % (
                                 replicon,
                                 self._replicon_dict[replicon]
                                 ['seq_start_pos'] + 1,
                                 self._replicon_dict[replicon]['seq_end_pos'],
                                 '\n'.join(df.apply(
                                     self._write_gff_entry, axis=1)),
                                 '\n' if not df.empty else ""))
    
    def _write_gff_entry(self, peak):
        return "%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=%s:peak_%s" % (
            peak["replicon"],
            "PEAKachu",
            "peak_region",
            peak["peak_start"],
            peak["peak_end"],
            peak["peak_strand"],
            peak["replicon"],
            peak.name + 1)
        
    def _plot_initial_peaks(self, unsig_base_means, unsig_fcs,
                            sig_base_means, sig_fcs):
        # MA plot
        plt.plot(np.log10(unsig_base_means),
                 np.log2(unsig_fcs), ".",
                 markersize=2.0, alpha=0.3)
        plt.plot(np.log10(sig_base_means),
                 np.log2(sig_fcs), ".",
                 markersize=2.0, color="red", alpha=0.3)
        plt.axhline(y=np.median(np.log2(unsig_fcs.append(sig_fcs))))
        plt.axvline(x=np.median(np.log10(unsig_base_means.append(
                                         sig_base_means))))
        plt.title("MA_plot")
        plt.xlabel("log10 base mean")
        plt.ylabel("log2 fold-change")
        plt.savefig("%s/MA_plot.png" % (self._output_folder), dpi=600)
        plt.close()
        # HexBin plot
        df = pd.DataFrame({'log10 base mean': np.log10(unsig_base_means.append(
            sig_base_means)), 'log2 fold-change': np.log2(unsig_fcs.append(
                sig_fcs))})
        df.plot(kind='hexbin', x='log10 base mean',
                y='log2 fold-change', gridsize=50, bins='log')
        plt.axhline(y=np.median(np.log2(unsig_fcs.append(sig_fcs))))
        plt.axvline(x=np.median(np.log10(unsig_base_means.append(
                                         sig_base_means))))
        plt.title("HexBin_plot")
        plt.savefig("%s/HexBin_plot.pdf" % (self._output_folder))
        plt.close()
        
    def _get_overlap(self, peak_start, peak_end, feature_start, feature_end):
        return max(
            0, min(peak_end, feature_end) - max(peak_start, feature_start) + 1)
