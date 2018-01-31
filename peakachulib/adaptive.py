import sys
from os.path import basename, splitext, isfile, exists
from os import makedirs
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.robust.scale import mad
import pandas as pd
import json
from concurrent import futures
from peakachulib.library import Library
from peakachulib.deseq2 import DESeq2Runner
from peakachulib.intersection import Intersecter, Interval
from time import time
from collections import OrderedDict
from subprocess import Popen, PIPE
from copy import deepcopy


class AdaptiveApproach(object):
    '''
    This class is used for peak detection via predefining peaks based on shape
    and subsequent comparison to a control
    '''
    def __init__(self, replicon_dict, max_proc, padj_threshold, mad_multiplier,
                 fc_cutoff, output_folder):
        self._lib_dict = OrderedDict()
        self._replicon_dict = replicon_dict  # own copy of replicon_dict
        self._max_proc = max_proc
        self._padj_threshold = padj_threshold
        self._mad_multiplier = mad_multiplier
        self._fc_cutoff = fc_cutoff
        self._output_folder = output_folder
        if not exists(self._output_folder):
            makedirs(self._output_folder)

    def init_libraries(self, paired_end, max_insert_size, ctr_libs,
                       exp_libs):
        self._paired_end = paired_end
        self._max_insert_size = max_insert_size
        self._ctr_lib_list = [splitext(basename(lib_file))[0]
                              for lib_file in ctr_libs]
        self._exp_lib_list = [splitext(basename(lib_file))[0]
                              for lib_file in exp_libs]
        # add libs to lib_dict
        for lib_file in exp_libs + ctr_libs:
            if not isfile(lib_file):
                sys.stderr.write(
                    "ERROR: The library file {} does not exist.\n".format(
                        lib_file))
                sys.exit(1)
            self._lib_dict[splitext(basename(lib_file))[0]] = Library(
                paired_end, max_insert_size, lib_file,
                deepcopy(self._replicon_dict))
        self._lib_names_list = list(self._lib_dict.keys())
        print("The following libraries were initialized:\n"
              "# Experiment libraries\n{0}\n"
              "# Control libraries\n{1}".format(
                  '\n'.join(self._exp_lib_list),
                  '\n'.join(self._ctr_lib_list)))

    def generate_combined_bed_file(self):
        # execute read conversion in parallel
        print("** Converting reads to bed format for {} libraries...".format(
            len(self._exp_lib_list)), flush=True)
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
                print("{} generated an exception: {}".format(lib_name, exc),
                      flush=True)
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
            split_index.loc[:, ["start", "end"]] = split_index.loc[
                :, ["start", "end"]].apply(pd.to_numeric)
            del self._replicon_dict[replicon]["reads"]["index"]
            self._replicon_dict[replicon]["reads"] = split_index.join(
                self._replicon_dict[replicon]["reads"]).sort_values(
                    ["strand", "start", "end"], ascending=[False, True, True])
            self._replicon_dict[replicon]["reads"]["replicon"] = replicon
            self._replicon_dict[replicon]["reads"]["tag_id"] = (
                self._replicon_dict[replicon]["reads"].index + 1).map(
                'tag_{:.0f}'.format)
            self._replicon_dict[replicon]["reads"] = self._replicon_dict[
                replicon]["reads"].loc[:,
                                       ["replicon",
                                        "start",
                                        "end",
                                        "tag_id",
                                        "count",
                                        "strand"]]
            # create blockbuster input folder if it does not exist
            self._blockbuster_input_folder = "{}/blockbuster_input".format(
                    self._output_folder)
            if not exists(self._blockbuster_input_folder):
                makedirs(self._blockbuster_input_folder)
            self._replicon_dict[replicon]["reads"].to_csv(
                "{}/{}_sorted_reads_for_blockbuster.bed".format(
                    self._blockbuster_input_folder, replicon),
                sep='\t', header=False, index=False, encoding='utf-8')
        t_end = time()
        print("Reads converted to bed format in {} seconds.\n".format(
            t_end-t_start), flush=True)

    def run_blockbuster(self, blockbuster_path):
        self._blockbuster_output = ""
        for replicon in sorted(self._replicon_dict):
            self._replicon_dict[replicon][
                "blockbuster"] = self._blockbuster_worker(blockbuster_path,
                                                          replicon)
            print("blockbuster for replicon {} exited with status {}".format(
                replicon, self._replicon_dict[replicon]["blockbuster"][
                    "returncode"]))
            if not self._replicon_dict[replicon]["blockbuster"][
                    "returncode"] == 0:
                print(self._replicon_dict[replicon]["blockbuster"]["error"])
                sys.exit(1)
            self._blockbuster_output += self._replicon_dict[replicon][
                "blockbuster"]["output"]
            del self._replicon_dict[replicon]["blockbuster"]
        with open("{}/blockbuster.txt".format(self._output_folder),
                  'w') as blockbuster_fh:
            blockbuster_fh.write(self._blockbuster_output)

    def _blockbuster_worker(self, blockbuster_path, replicon):
        p = Popen(
            ["blockbuster.x", "-minBlockHeight", "10", "-print", "1",
             "-distance", "1",
             "{}/{}_sorted_reads_for_blockbuster.bed".format(
                 self._blockbuster_input_folder, replicon)], stdout=PIPE,
            stderr=PIPE, universal_newlines=True)
        output, error = p.communicate()
        returncode = p.returncode
        return {"output": output, "error": error, "returncode": returncode}

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
            block_df = pd.DataFrame(
                blocks, columns=["blockNb", "blockChrom", "blockStart",
                                 "blockEnd", "blockStrand", "blockExpression",
                                 "readCount"])
            block_df[["blockNb", "blockStart", "blockEnd", "blockExpression",
                      "readCount"]] = block_df[
                    ["blockNb", "blockStart", "blockEnd", "blockExpression",
                     "readCount"]].apply(pd.to_numeric)
            peak_df = self._split_cluster_peaks(block_df, cluster_expr,
                                                peak_df, min_cluster_expr_frac,
                                                min_block_overlap,
                                                min_max_block_expr_frac)
        if peak_df.empty:
            return
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
        self._generate_peak_counts()
        self._peak_df = pd.DataFrame()
        for replicon in sorted(self._replicon_dict):
            if self._replicon_dict[replicon]["peak_df"].empty:
                continue
            self._replicon_dict[replicon]["peak_df"]["replicon"] = replicon
            for lib_name, lib in self._lib_dict.items():
                self._replicon_dict[replicon][
                    "peak_df"][lib_name] = lib.replicon_dict[
                        replicon]["peak_counts"]
                del lib.replicon_dict[replicon]["peak_counts"]
            # add pseudocounts
            # self._replicon_dict[
            #    replicon]["peak_df"].loc[:, self._lib_names_list] += 1.0
            self._peak_df = self._peak_df.append(self._replicon_dict[replicon][
                "peak_df"], ignore_index=True)

    def _generate_peak_counts(self):
        print("** Peak read counting started for {} libraries...".format(
            len(self._lib_dict)), flush=True)
        t_start = time()
        for lib_name, lib in self._lib_dict.items():
            print(lib_name, flush=True)
            for replicon in self._replicon_dict:
                lib.replicon_dict[replicon]["peak_df"] = self._replicon_dict[
                    replicon]["peak_df"]
            lib.count_reads_for_peaks()
        t_end = time()
        print("Peak read counting finished in {} seconds.".format(
            t_end-t_start), flush=True)

    def run_deseq2_analysis(self, size_factors, pairwise_replicates):
        count_df = self._peak_df.loc[:, self._exp_lib_list +
                                     self._ctr_lib_list]
        deseq2_runner = DESeq2Runner(count_df)
        result_df, self._size_factors = deseq2_runner.run_deseq2(
            self._exp_lib_list, self._ctr_lib_list, size_factors,
            pairwise_replicates)
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
            "{}/initial_peaks.csv".format(self._output_folder),
            sep='\t', na_rep='NA', index=False, encoding='utf-8')
        # filter peaks
        print("* Filtering peaks...", flush=True)
        sig_peak_df = self._filter_peaks(self._peak_df)
        unsig_peak_df = self._peak_df[~self._peak_df.index.isin(
            sig_peak_df.index)]
        # plot peaks
        print("* Plotting initial peaks...", flush=True)
        t_start = time()
        self._plot_initial_peaks(unsig_peak_df.baseMean,
                                 np.power(2.0, unsig_peak_df.log2FoldChange),
                                 sig_peak_df.baseMean,
                                 np.power(2.0, sig_peak_df.log2FoldChange))
        t_end = time()
        print("Plotting took {} seconds.".format(t_end-t_start), flush=True)
        self._peak_df = sig_peak_df

    def run_analysis_without_replicates(self, size_factors):
        # check if size factors were defined and otherwise calculate them based
        # on DESeq normalization
        if size_factors is None:
            deseq2_runner = DESeq2Runner(self._peak_df[self._lib_names_list])
            self._size_factors = deseq2_runner.calc_size_factors()
        else:
            self._size_factors = size_factors
        # normalize counts
        self._peak_df[self._lib_names_list] = self._peak_df[
            self._lib_names_list].div(self._size_factors, axis='columns')
        # calculate base means for all peaks
        self._peak_df["base_means"] = self._peak_df.loc[
            :, self._lib_names_list].mean(axis=1)
        # calculate fcs for all peaks
        self._peak_df["fold_change"] = (
            self._peak_df.loc[:, self._exp_lib_list].sum(axis=1) /
            self._peak_df.loc[:, self._ctr_lib_list].sum(axis=1))
        # write initial peaks
        peak_columns = (["replicon",
                         "peak_start",
                         "peak_end",
                         "peak_strand"] +
                        [lib_name for lib_name in self._lib_dict] +
                        ["base_means",
                         "fold_change"])
        self._peak_df.loc[:, peak_columns].to_csv(
            "{}/initial_peaks.csv".format(self._output_folder),
            sep='\t', na_rep='NA', index=False, encoding='utf-8')
        # filter peaks
        print("* Filtering peaks...", flush=True)
        sig_peak_df = self._filter_peaks_without_replicates(self._peak_df)
        unsig_peak_df = self._peak_df[~self._peak_df.index.isin(
            sig_peak_df.index)]
        # plot peaks
        print("* Plotting initial peaks...", flush=True)
        t_start = time()
        self._plot_initial_peaks(unsig_peak_df.base_means,
                                 unsig_peak_df.fold_change,
                                 sig_peak_df.base_means,
                                 sig_peak_df.fold_change)
        t_end = time()
        print("Plotting took {} seconds.".format(t_end-t_start), flush=True)
        self._peak_df = sig_peak_df

    def run_analysis_without_control(self, size_factors):
        # check if size factors were defined and otherwise calculate them based
        # on DESeq normalization
        if size_factors is None:
            deseq2_runner = DESeq2Runner(self._peak_df[self._lib_names_list])
            self._size_factors = deseq2_runner.calc_size_factors()
        else:
            self._size_factors = size_factors
        # normalize counts
        self._peak_df[self._lib_names_list] = self._peak_df[
            self._lib_names_list].div(self._size_factors, axis='columns')
        # calculate base means for all peaks
        self._peak_df["base_means"] = self._peak_df.loc[
            :, self._lib_names_list].mean(axis=1)
        # write initial peaks
        peak_columns = (["replicon",
                         "peak_start",
                         "peak_end",
                         "peak_strand"] +
                        [lib_name for lib_name in self._lib_dict] +
                        ["base_means"])
        self._peak_df.loc[:, peak_columns].to_csv(
            "{}/initial_peaks.csv".format(self._output_folder),
            sep='\t', na_rep='NA', index=False, encoding='utf-8')
        # filter peaks
        print("* Filtering peaks...", flush=True)
        self._peak_df = self._filter_peaks_without_control(self._peak_df)

    def _filter_peaks(self, df):
        # calculate mad for original data frame
        median_abs_dev_from_zero = mad(df.loc[:, self._exp_lib_list].mean(
            axis=1), center=0.0)
        # padj filter
        print("Removing peaks based on padj from DataFrame with {} rows...".
              format(len(df)), flush=True)
        t_start = time()
        df = df.query('padj < @self._padj_threshold')
        t_end = time()
        print("Removal took {} seconds. DataFrame contains now {} rows.".
              format((t_end-t_start), len(df)), flush=True)
        if df.empty:
            return df
        # minimum expression cutoff based on mean over experiment libraries
        print("Removing peaks based on mad cutoff from DataFrame "
              "with {} rows...".format(len(df)), flush=True)
        t_start = time()
        min_expr = (self._mad_multiplier * median_abs_dev_from_zero)
        print("Minimal peak expression based on mean over RIP/CLIP "
              "libraries:" "{} (MAD from zero: {})".format(
                  min_expr, median_abs_dev_from_zero), flush=True)
        df = df.loc[df.loc[:, self._exp_lib_list].mean(axis=1) >= min_expr, :]
        t_end = time()
        print("Removal took {} seconds. DataFrame contains now {} rows.".
              format((t_end-t_start), len(df)), flush=True)
        if df.empty:
            return df
        # minimum fold change
        print("Removing peaks based on minimum fold change "
              "from DataFrame with {} rows...".format(len(df)), flush=True)
        t_start = time()
        log2_fc_cutoff = np.log2(self._fc_cutoff)
        df = df.query('log2FoldChange >= @log2_fc_cutoff')
        t_end = time()
        print("Removal took {} seconds. DataFrame contains now {} rows.".
              format((t_end-t_start), len(df)), flush=True)
        return df

    def _filter_peaks_without_replicates(self, df):
        # calculate mad for original data frame
        median_abs_dev_from_zero = mad(df.loc[:, self._exp_lib_list].mean(
            axis=1), center=0.0)
        # minimum expression cutoff based on mean over experiment libraries
        print("Removing peaks based on mad cutoff from DataFrame "
              "with {} rows...".format(len(df)), flush=True)
        t_start = time()
        min_expr = (self._mad_multiplier * median_abs_dev_from_zero)
        print("Minimal peak expression based on mean over RIP/CLIP "
              "libraries:" "{} (MAD from zero: {})".format(
                  min_expr, median_abs_dev_from_zero), flush=True)
        df = df.loc[df.loc[:, self._exp_lib_list].mean(axis=1) >= min_expr, :]
        t_end = time()
        print("Removal took {} seconds. DataFrame contains now {} rows.".
              format((t_end-t_start), len(df)), flush=True)
        if df.empty:
            return df
        # minimum fold change
        print("Removing windows based on minimum fold change from DataFrame "
              "with {} rows...".format(len(df)), flush=True)
        t_start = time()
        df = df.query('fold_change >= @self._fc_cutoff')
        t_end = time()
        print("Removal took {} seconds. DataFrame contains now {} rows.".
              format((t_end-t_start), len(df)), flush=True)
        return df

    def _filter_peaks_without_control(self, df):
        # calculate mad for original data frame
        median_abs_dev_from_zero = mad(df.loc[:, self._exp_lib_list].mean(
            axis=1), center=0.0)
        # minimum expression cutoff based on mean over experiment libraries
        print("Removing peaks based on mad cutoff from DataFrame "
              "with {} rows...".format(len(df)), flush=True)
        t_start = time()
        min_expr = (self._mad_multiplier * median_abs_dev_from_zero)
        print("Minimal peak expression based on mean over RIP/CLIP "
              "libraries:" "{} (MAD from zero: {})".format(
                  min_expr, median_abs_dev_from_zero), flush=True)
        df = df.loc[df.loc[:, self._exp_lib_list].mean(axis=1) >= min_expr, :]
        t_end = time()
        print("Removal took {} seconds. DataFrame contains now {} rows.".
              format((t_end-t_start), len(df)), flush=True)
        return df

    def write_output(self):
        # write parameters to json file
        self._write_parameters()
        # create peak table folder if it does not exist
        peak_table_folder = "{}/peak_tables".format(self._output_folder)
        if not exists(peak_table_folder):
            makedirs(peak_table_folder)
        peak_columns = (["replicon",
                         "peak_id",
                         "peak_start",
                         "peak_end",
                         "peak_strand"] +
                        [lib_name for lib_name in self._lib_dict])
        if len(self._exp_lib_list) > 1 and len(self._ctr_lib_list) > 1:
            peak_columns += ["baseMean",
                             "log2FoldChange",
                             "lfcSE",
                             "stat",
                             "pvalue",
                             "padj"]
        elif self._ctr_lib_list:
            peak_columns += ["base_means",
                             "fold_change"]
        else:
            peak_columns += ["base_means"]
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
            self._replicon_dict[replicon]["peak_df"].sort_values(
                ["replicon", "peak_start"], inplace=True)
            self._replicon_dict[replicon]["peak_df"].reset_index(
                drop=True, inplace=True)
            self._replicon_dict[replicon]["peak_df"].loc[:, "peak_id"] = (
                self._replicon_dict[replicon]["peak_df"].index + 1)
            feature_tree = self._build_feature_tree(replicon)
            for peak in self._replicon_dict[replicon]["peak_df"].to_dict(
                    'records'):
                overlapping_features = self._find_overlapping_features(
                        peak,
                        feature_tree)
                for match in overlapping_features:
                    entry_dict = peak.copy()
                    entry_dict.update(match)
                    output_df = output_df.append(
                        pd.DataFrame(entry_dict, index=[0],
                                     columns=peak_columns+feature_columns),
                        ignore_index=True)
            # write peak table for replicon
            output_df.to_csv(
                "{}/peaks_{}.csv".format(peak_table_folder, replicon),
                sep='\t', na_rep='NA', index=False, encoding='utf-8')
            # write peak gff file for replicon
            self._write_gff_file(replicon, self._replicon_dict[replicon]
                                 ["peak_df"])
        del self._peak_df

    def _build_feature_tree(self, replicon):
        interval_tree = Intersecter()
        for feature in self._replicon_dict[replicon]["features"]:
            interval_tree.add_interval(
                Interval(feature["start"],
                         feature["end"],
                         value=feature,
                         strand=feature["strand"]))
        return interval_tree

    def _find_overlapping_features(self, peak, feature_tree):
        overlapping_features = []
        all_overlapping = feature_tree.find(peak["peak_start"],
                                            peak["peak_end"])
        for feature in all_overlapping:
            if not peak["peak_strand"] == feature.strand:
                continue
            overlap = self._get_overlap(peak["peak_start"], peak["peak_end"],
                                        feature.start, feature.end)
            overlapping_features.append({
                "feature_type": feature.value["type"],
                "feature_start": feature.start,
                "feature_end": feature.end,
                "feature_strand": feature.strand,
                "feature_locus_tag": feature.value["locus_tag"],
                "feature_name": feature.value["Name"],
                "subfeature_type": feature.value["subfeature_type"] if (
                    "subfeature_type" in feature.value) else None,
                "feature_product": feature.value["product"],
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

    def _write_gff_file(self, replicon, df):
        # create peak annotation folder if it does not exist
        peak_anno_folder = "{}/peak_annotations".format(self._output_folder)
        if not exists(peak_anno_folder):
            makedirs(peak_anno_folder)
        with open("{}/peaks_{}.gff".format(
                peak_anno_folder, replicon), 'w') as out_gff_fh:
            out_gff_fh.write("##gff-version 3\n"
                             "#!gff-spec-version 1.20\n"
                             "##sequence-region {} {} {}\n"
                             "{}{}"
                             "###\n".format(
                                 replicon,
                                 self._replicon_dict[replicon]
                                 ['seq_start_pos'] + 1,
                                 self._replicon_dict[replicon]['seq_end_pos'],
                                 '\n'.join(
                                     df.apply(self._write_gff_entry, axis=1)),
                                 '\n' if not df.empty else ""))

    def _write_gff_entry(self, peak):
        return "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\tID={}:peak_{}".format(
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
        # create plot folder if it does not exist
        plot_folder = "{}/plots".format(self._output_folder)
        if not exists(plot_folder):
            makedirs(plot_folder)
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
        plt.title("Initial_peaks_MA_plot")
        plt.xlabel("log10 base mean")
        plt.ylabel("log2 fold-change")
        plt.savefig("{}/Initial_peaks_MA_plot.png".format(plot_folder),
                    dpi=600)
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
        plt.title("Initial_peaks_HexBin_plot")
        plt.savefig("{}/Initial_peaks_HexBin_plot.pdf".format(plot_folder))
        plt.close()

    def _get_overlap(self, peak_start, peak_end, feature_start, feature_end):
        return max(
            0, min(peak_end, feature_end) - max(peak_start, feature_start) + 1)

    def _write_parameters(self):
        parameter_dict = {"max_insert_size": self._max_insert_size,
                          "paired_end": self._paired_end,
                          "libraries": {}}
        for lib, size_factor in zip(
                self._lib_dict.values(), self._size_factors):
            parameter_dict["libraries"][lib.lib_name] = {
                "bam_file": lib.bam_file,
                "size_factor": size_factor}
        with open("{}/parameters.json".format(self._output_folder),
                  'w') as parameter_fh:
            json.dump(parameter_dict, parameter_fh, indent=4, sort_keys=True)
