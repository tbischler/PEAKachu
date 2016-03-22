import sys
from os.path import basename, splitext, isfile, exists
from os import makedirs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.robust.scale import mad
from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd
matplotlib.style.use('seaborn-colorblind')
font = {'family': 'sans-serif', 'size': 7}
matplotlib.rc('font', **font)
from concurrent import futures
from peakachulib.library import Library
from peakachulib.tmm import TMM
from peakachulib.gtest import GTest
from peakachulib.wiggle import WiggleWriter
from time import time
from collections import OrderedDict
from copy import deepcopy


class WindowApproach(object):
    '''
    This class is used for peak detection via a sliding window approach
    '''
    def __init__(self, w_size, step_size, replicon_dict, max_proc, norm_method,
                 size_factors, het_p_val_threshold, rep_pair_p_val_threshold,
                 padj_threshold, mad_multiplier, fc_cutoff,
                 pairwise_replicates, output_folder):
        self._lib_dict = OrderedDict()
        self._replicon_dict = replicon_dict  # own copy of replicon_dict
        self._max_proc = max_proc
        self._w_size = w_size
        self._step_size = step_size
        self._norm_method = norm_method
        self._size_factors = size_factors
        self._het_p_val_threshold = het_p_val_threshold
        self._rep_pair_p_val_threshold = rep_pair_p_val_threshold
        self._padj_threshold = padj_threshold
        self._mad_multiplier = mad_multiplier
        self._fc_cutoff = fc_cutoff
        self._pairwise_replicates = pairwise_replicates
        self._output_folder = output_folder + "/window_approach"
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

    def generate_window_counts(self):
        self._generate_windows()
        # execute read counting in parallel
        print("** Window read counting started for %s libraries..." % len(
            self._lib_dict), flush=True)
        t_start = time()
        with futures.ProcessPoolExecutor(
                max_workers=self._max_proc) as executor:
            future_to_lib_name = {
                executor.submit(lib.count_reads_for_windows):
                lib.lib_name for lib in self._lib_dict.values()}
        for future in futures.as_completed(future_to_lib_name):
            lib_name = future_to_lib_name[future]
            try:
                self._lib_dict[lib_name].replicon_dict = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (lib_name, exc),
                      flush=True)
        t_end = time()
        print("Window read counting finished in %s seconds.\n" % (
            t_end-t_start), flush=True)
        print("** Generating data frames and filtering windows...", flush=True)
        t_start = time()
        self._convert_to_data_frame()
        t_end = time()
        print("Data frame generation and filtering finished in %s seconds.\n"
              % (t_end-t_start), flush=True)

    def _generate_windows(self):
        for replicon in self._replicon_dict:
            self._replicon_dict[replicon]["window_list"] = []
            for w_start in range(
                    self._replicon_dict[replicon]['seq_start_pos'],
                    self._replicon_dict[replicon]['seq_end_pos'],
                    self._step_size):
                w_end = w_start + self._w_size
                if w_end > self._replicon_dict[replicon]['seq_end_pos']:
                    w_end = self._replicon_dict[replicon]['seq_end_pos']
                    self._replicon_dict[replicon]["window_list"].append(
                        (w_start, w_end))
                    break
                self._replicon_dict[replicon]["window_list"].append(
                    (w_start, w_end))

    def perform_g_test_with_repl_for_windows(self):
        if self._window_df.empty:
            return
        p_values = self._window_df.apply(self._single_g_test, axis=1)
        padj_values = p_values.loc[
            :, ["pooled_G_p_value", "total_G_p_value"]].apply(
                self._correct_p_values, axis=0)
        padj_values.columns = [col_name.replace("p_value", "padj_value")
                               for col_name in padj_values.columns]
        padj_values.index = p_values.index
        significance = pd.concat([p_values, padj_values], axis=1).apply(
            self._check_significance_with_repl, axis=1)
        significance.name = "significant"
        self._window_df = pd.concat([self._window_df, p_values, padj_values,
                                     significance], axis=1)

    def perform_g_test_without_repl_for_windows(self):
        if self._window_df.empty:
            return
        p_values = self._window_df.apply(self._single_g_test, axis=1)
        padj_values = p_values.apply(self._correct_p_values, axis=0)
        padj_values.columns = [
            col_name.replace("p_value", "padj_value")
            for col_name in padj_values.columns]
        padj_values.index = p_values.index
        significance = pd.concat(
            [p_values, padj_values],
            axis=1).apply(self._check_significance_without_repl, axis=1)
        significance.name = "significant"
        self._window_df = pd.concat([self._window_df, p_values, padj_values,
                                     significance], axis=1)

    def _perform_g_test_with_repl_for_peaks(self, df):
        if df.empty:
            return df
        p_values = df.apply(self._single_g_test, axis=1)
        padj_values = p_values.loc[
            :, ["pooled_G_p_value", "total_G_p_value"]].apply(
                self._correct_p_values, axis=0)
        padj_values.columns = [
            col_name.replace("p_value", "padj_value")
            for col_name in padj_values.columns]
        padj_values.index = p_values.index
        df = pd.concat([df, p_values, padj_values], axis=1)
        return df

    def _perform_g_test_without_repl_for_peaks(self, df):
        if df.empty:
            return df
        p_values = df.apply(self._single_g_test, axis=1)
        padj_values = p_values.apply(self._correct_p_values, axis=0)
        padj_values.columns = [
            col_name.replace("p_value", "padj_value")
            for col_name in padj_values.columns]
        padj_values.index = p_values.index
        df = pd.concat([df, p_values, padj_values], axis=1)
        return df

    def combine_peaks_and_recalculate_values(self):
        for replicon in self._replicon_dict:
            # forward strand
            peak_df_forward = self._combine_windows(
                self._window_df[
                    (self._window_df["replicon"] == replicon) & (
                        self._window_df["strand"] == "+")])
            if not peak_df_forward.empty:
                peak_df_forward["peak_strand"] = '+'
            # reverse strand
            peak_df_reverse = self._combine_windows(self._window_df[
                (self._window_df["replicon"] == replicon) & (
                    self._window_df["strand"] == "-")])
            if not peak_df_reverse.empty:
                peak_df_reverse["peak_strand"] = '-'
            # combined
            self._replicon_dict[replicon]["peak_df"] = pd.concat([
                peak_df_forward, peak_df_reverse], axis=0, ignore_index=True)
        for lib in self._lib_dict.values():
            for replicon in self._replicon_dict:
                lib.replicon_dict[replicon]["peak_df"] = self._replicon_dict[
                    replicon]["peak_df"]
        self._generate_peak_counts()

        for replicon in sorted(self._replicon_dict):
            for lib_name, lib in self._lib_dict.items():
                self._replicon_dict[
                    replicon]["peak_df"][lib_name] = lib.replicon_dict[
                        replicon]["peak_counts"]
            # add pseudocounts
            self._replicon_dict[replicon]["peak_df"][
                self._lib_names_list] += 1.0
            # normalize counts
            self._replicon_dict[
                replicon]["peak_df"][
                    self._lib_names_list] = self._replicon_dict[
                    replicon]["peak_df"][
                    self._lib_names_list].div(self._size_factors,
                                              axis='columns')
            if len(self._exp_lib_list) > 1:
                self._replicon_dict[
                    replicon][
                        "peak_df"] = self._perform_g_test_with_repl_for_peaks(
                        self._replicon_dict[replicon]["peak_df"])
            else:
                self._replicon_dict[replicon][
                    "peak_df"] = self._perform_g_test_without_repl_for_peaks(
                        self._replicon_dict[replicon]["peak_df"])
            # calculate base means for all peaks
            self._replicon_dict[replicon][
                "peak_df"]["base_means"] = self._replicon_dict[
                    replicon]["peak_df"].loc[
                        :, self._lib_names_list].mean(axis=1)
            # calculate fcs for all peaks
            self._replicon_dict[replicon]["peak_df"]["fold_change"] = (
                self._replicon_dict[
                    replicon]["peak_df"].loc[:, self._exp_lib_list].sum(
                        axis=1) / self._replicon_dict[replicon][
                            "peak_df"].loc[:, self._ctr_lib_list].sum(axis=1))

    def _generate_peak_counts(self):
        print("** Peak read counting started for %s libraries..." % len(
            self._lib_dict), flush=True)
        t_start = time()
        for lib_name, lib in self._lib_dict.items():
            print(lib_name)
            lib.count_reads_for_peaks()
        t_end = time()
        print("Peak read counting finished in %s seconds." % (t_end-t_start),
              flush=True)

    def write_output(self):
        peak_columns = (["replicon",
                         "peak_id",
                         "peak_start",
                         "peak_end",
                         "peak_strand"] +
                        [lib_name for lib_name in self._lib_dict] +
                        ["base_means",
                         "fold_change"])
        if len(self._exp_lib_list) > 1:
            peak_columns += ["heterogenous_G_p_value",
                             "replicate_G_p_values",
                             "pooled_G_padj_value",
                             "total_G_padj_value"]
        else:
            peak_columns.append("single_G_padj_value")
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
            if self._replicon_dict[replicon]["peak_df"].empty:
                continue
            output_df = pd.DataFrame()
            self._replicon_dict[replicon]["peak_df"]["replicon"] = replicon
            self._replicon_dict[replicon]["peak_df"].sort_values(
                ["replicon", "peak_start"], inplace=True)
            self._replicon_dict[replicon]["peak_df"].reset_index(
                drop=True, inplace=True)
            self._replicon_dict[replicon]["peak_df"]["peak_id"] = (
                self._replicon_dict[replicon]["peak_df"].index + 1)
            for peak in self._replicon_dict[replicon]["peak_df"].to_dict(
                    'records'):
                overlapping_features = self._find_overlapping_features(peak)
                for match in overlapping_features:
                    entry_dict = peak.copy()
                    entry_dict.update(match)
                    output_df = output_df.append(pd.DataFrame(entry_dict,
                        index=[0], columns=peak_columns+feature_columns),
                        ignore_index=True)
            output_df.to_csv(
                "%s/peaks_%s.csv" % (self._output_folder, replicon),
                sep='\t', na_rep='NA', index=False, encoding='utf-8')

            self._write_gff_file(replicon, self._replicon_dict[replicon]
                                 ["peak_df"])
        # plot windows
        print("Plotting normalized windows...", flush=True)
        t_start = time()
        unsig_window_df = self._initial_window_df[
            ~self._initial_window_df.index.isin(self._window_df.index)]
        self._plot_initial_windows(unsig_window_df.base_means,
                                   unsig_window_df.fold_change,
                                   self._window_df.base_means,
                                   self._window_df.fold_change)
        t_end = time()
        print("Plotting took %s seconds." % (t_end-t_start), flush=True)

        if self._window_df.empty:
            return
        self._window_df.to_csv(
            "%s/windows_after_prefiltering.csv" % (self._output_folder),
            sep='\t', index=False, encoding='utf-8')

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
            future_to_lib_name = {executor.submit(
                self._generate_normalized_wiggle_file_for_lib, lib,
                size_factor, wiggle_folder
            ): lib.lib_name for lib, size_factor in zip(
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
            "%s_%s" % (lib.lib_name,
                       strand),
            open("%s/%s_div_by_%s_%s.wig" % (wiggle_folder, lib.lib_name,
                                             size_factor, strand), "w")))
            for strand in strand_dict.values()])
        for replicon in self._replicon_dict:
            for strand in strand_dict.keys():
                if strand == "-":
                    factor = -1.0/size_factor
                else:
                    factor = 1.0/size_factor
                try:
                    wiggle_writers[
                        strand_dict[strand]].write_replicons_coverages(
                            replicon, lib.replicon_dict[
                                replicon]["coverages"][strand], factor=factor)
                except Exception as exc:
                    print("Library %s, replicon %s, %s strand generated an"
                          "exception during coverage file generation: %s" %
                          (lib.lib_name, replicon, strand, exc), flush=True)
        for strand in strand_dict.values():
            wiggle_writers[strand].close_file()

    def _write_gff_file(self, replicon, df):
        with open("%s/peaks_%s.gff" % (
                self._output_folder, replicon), 'w') as out_gff_fh:
            out_gff_fh.write("##gff-version 3\n"
                             "#!gff-spec-version 1.20\n"
                             "##sequence-region %s %s %s\n"
                             "%s%s"
                             "###\n" % (
                                 replicon,
                                 self._replicon_dict[replicon]
                                 ['seq_start_pos'] + 1,
                                 self._replicon_dict[replicon]['seq_end_pos'],
                                 '\n'.join(
                                     df.apply(self._write_gff_entry, axis=1)),
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

    def _convert_to_data_frame(self):
        self._window_df = pd.DataFrame()
        for replicon in sorted(self._replicon_dict):
            for strand in ["+", "-"]:
                # add window positions to data frame
                row_number = len(self._replicon_dict[replicon]["window_list"])
                df = pd.concat([
                    pd.Series([replicon] * row_number),
                    pd.Series([strand] * row_number),
                    pd.Series([window[0]+1 for window in
                               self._replicon_dict[
                                   replicon]["window_list"]]),
                    pd.Series([window[1] for window in
                               self._replicon_dict[
                        replicon]["window_list"]])], axis=1)
                df.columns = ["replicon", "strand", "w_start", "w_end"]
                # add library counts to data frame
                for lib_name, lib in self._lib_dict.items():
                    df[lib_name] = (pd.Series(lib.replicon_dict[
                        replicon]["window_counts"].loc[:, strand]))
                self._window_df = self._window_df.append(df,
                                                         ignore_index=True)
        # remove windows without expression in any library
        print("Removing empty windows from DataFrame with %s rows..." % len(
            self._window_df.index), flush=True)
        t_start = time()
        self._window_df = self._window_df.loc[
            (self._window_df.loc[:, self._lib_names_list].sum(axis=1) > 0), :]
        t_end = time()
        print("Removal took %s seconds. DataFrame contains now %s rows." % (
            (t_end-t_start), len(self._window_df.index)), flush=True)
        if self._window_df.empty:
            print("**Dataframe empty**", flush=True)
            return
        # define size factors
        print("Calculating size factors and normalizing libraries...",
              flush=True)
        t_start = time()
        if self._norm_method == "tmm":
            # calc size factors based on tmm using windows with expression
            # in control
            tmm_df = self._window_df.loc[
                self._window_df.loc[:, self._ctr_lib_list].max(axis=1) > 0,
                self._lib_names_list]
            # if data frame with reads in the control is empty skip
            # normalization
            if tmm_df.empty:
                self._size_factors = pd.Series([1.0] * len(
                    self._lib_names_list),
                    index=self._lib_names_list)
            else:
                norm = TMM(tmm_df)
                self._size_factors = norm.calc_size_factors()
        elif self._norm_method == "count":
            # calc size factors based on library counts using windows with
            # expression in control
            count_df = self._window_df.loc[
                self._window_df.loc[:, self._ctr_lib_list].max(axis=1) > 0,
                self._lib_names_list]
            # if data frame with reads in the control is empty skip
            # normalization
            if count_df.empty:
                self._size_factors = pd.Series([1.0] * len(
                    self._lib_names_list),
                    index=self._lib_names_list)
            else:
                lib_sums = count_df.sum(axis=0)
                self._size_factors = lib_sums/lib_sums.max()
        else:
            self._size_factors = pd.Series(self._size_factors,
                                           index=self._lib_names_list)
        print("Size factors used for normalization\n%s" % (pd.Series.to_string(
            self._size_factors)), flush=True)
        # add pseudocounts
        self._window_df[self._lib_names_list] += 1.0
        # normalize counts
        self._window_df[self._lib_names_list] = self._window_df[
            self._lib_names_list].div(
                self._size_factors, axis='columns')
        t_end = time()
        print("Normalization took %s seconds." % (t_end-t_start), flush=True)
        # calculate base means for all windows
        print("Calculating base means and fold changes...", flush=True)
        t_start = time()
        self._window_df["base_means"] = self._window_df.loc[
            :, self._lib_names_list].mean(axis=1)
        # calculate fcs for all windows
        self._window_df["fold_change"] = (
            self._window_df.loc[:, self._exp_lib_list].sum(axis=1) /
            self._window_df.loc[:, self._ctr_lib_list].sum(axis=1))
        t_end = time()
        print("Calculation took %s seconds." % (t_end-t_start), flush=True)
        # write raw windows to file
        print("Writing normalized windows to file...", flush=True)
        t_start = time()
        self._window_df.to_csv("%s/raw_windows.csv" % (self._output_folder),
                               sep='\t', index=False, encoding='utf-8')
        t_end = time()
        print("Writing took %s seconds." % (t_end-t_start), flush=True)
        # filter windows
        print("* Filtering windows...", flush=True)
        self._initial_window_df = self._window_df.copy()
        self._window_df = self._filter_windows_df(self._window_df)

    def _plot_initial_windows(self, unsig_base_means, unsig_fcs,
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

    def _filter_windows_df(self, df):
        ''' This function filters the windows in a data frame by minimum
            expression based on a MAD cutoff and requires higher expression
            in the experiment libs than in the controls
        '''
        # remove windows where not all experiment libs show expression:
        #   expression = 1/size_factor ( = pseudocount)
        print("Removing windows where not all experiment libs show "
              "expression from DataFrame with %s rows..." % len(df),
              flush=True)
        t_start = time()
        for exp_lib in self._exp_lib_list:
            exp_lib_zero_count = 1/self._size_factors[exp_lib]
            df = df.loc[(df.loc[:, exp_lib] > exp_lib_zero_count), :]
        t_end = time()
        print("Removal took %s seconds. DataFrame contains now %s rows." % (
            (t_end-t_start), len(df)), flush=True)
        if df.empty:
            return df
        # minimum expression cutoff based on mean over experiment libraries
        print("Removing windows based on mad cutoff from DataFrame "
              "with %s rows..." % len(df), flush=True)
        t_start = time()
        median_abs_dev_from_zero = mad(df.loc[:, self._exp_lib_list].mean(
            axis=1), center=0.0)
        min_expr = (self._mad_multiplier * median_abs_dev_from_zero)
        print("Minimal window expression based on mean over RIP/CLIP "
              "libraries:" "%s (MAD from zero: %s)" % (
                  min_expr, median_abs_dev_from_zero), flush=True)
        df = df.loc[df.loc[:, self._exp_lib_list].mean(axis=1) >= min_expr, :]
        t_end = time()
        print("Removal took %s seconds. DataFrame contains now %s rows." % (
            (t_end-t_start), len(df)), flush=True)
        if df.empty:
            return df
        print("Removing windows where experiment expression is lower than "
              "control expression from DataFrame with %s rows..." % len(df),
              flush=True)
        t_start = time()
        if self._pairwise_replicates:
            # experiment expression must be larger than respective control
            # for each library pair
            for exp_lib, ctr_lib in zip(
                    self._exp_lib_list, self._ctr_lib_list):
                df = df.loc[(df.loc[:, exp_lib] > df.loc[:, ctr_lib]), :]
        else:
            # minimum experiment expression larger than maximum
            # control expression
            df = df.loc[df.loc[:, self._exp_lib_list].min(
                axis=1) > df.loc[:, self._ctr_lib_list].max(axis=1), :]
        t_end = time()
        print("Removal took %s seconds. DataFrame contains now %s rows." % (
            (t_end-t_start), len(df)), flush=True)
        if df.empty:
            return df
        # minimum fold change
        print("Removing windows based on minimum fold change from DataFrame "
              "with %s rows..." % len(df), flush=True)
        t_start = time()
        df = df.query('fold_change >= @self._fc_cutoff')
        t_end = time()
        print("Removal took %s seconds. DataFrame contains now %s rows." % (
            (t_end-t_start), len(df)), flush=True)
        return df

    def _single_g_test(self, counts):
        ctr_counts = counts[self._ctr_lib_list]
        ctr_counts = ctr_counts.reset_index(drop=True)
        exp_counts = counts[self._exp_lib_list]
        exp_counts = exp_counts.reset_index(drop=True)
        g_test = GTest(ctr_counts, exp_counts, self._pairwise_replicates)
        if len(exp_counts) > 1:
            return pd.Series(g_test.run_with_repl())
        else:
            return pd.Series(g_test.run_without_repl())

    def _correct_p_values(self, p_values):
        return pd.Series(multipletests(p_values, method="fdr_bh")[1])

    def _check_significance_with_repl(self, p_and_padj_values):
        replicate_G_p_values = pd.Series(p_and_padj_values[
            "replicate_G_p_values"].split('/')).astype('float')
        if (p_and_padj_values.loc["heterogenous_G_p_value"]
                >= self._het_p_val_threshold and
                p_and_padj_values.loc["pooled_G_padj_value"]
                < self._padj_threshold):
            return True
        if ((p_and_padj_values.loc[
                "total_G_padj_value"] < self._padj_threshold)
                and ((replicate_G_p_values <
                      self._rep_pair_p_val_threshold).all())):
            return True
        return False

    def _check_significance_without_repl(self, p_and_padj_values):
        if (p_and_padj_values.loc["single_G_padj_value"]
                < self._padj_threshold):
            return True
        return False

    def _combine_windows(self, df):
        peak_list = []
        peak = {"peak_start": None, "peak_end": None}
        for index, window in df.iterrows():
            # significant window
            if window.loc["significant"]:
                # start new peak region if no peak region was started before
                if peak["peak_start"] is None:
                    peak["peak_start"] = window.loc["w_start"]
                    peak["peak_end"] = window.loc["w_end"]
                # start new peak region while still inside previous peak region
                #   *this is due to gaps in the windows caused by pre-filtering
                #   *add previous peak region to output
                # +1 -> combine adjacent peaks
                elif window.loc["w_start"] > peak["peak_end"] + 1:
                    peak_list.append(deepcopy(peak))
                    peak["peak_start"] = window.loc["w_start"]
                    peak["peak_end"] = window.loc["w_end"]
                # elongate peak if window overlaps
                else:
                    peak["peak_end"] = window.loc["w_end"]
            # non-significant window
            else:
                # jump to next window if outside of peak region
                # or current position upstream of peak end
                # +1 -> combine adjacent peaks
                if (peak["peak_start"] is None or window.loc[
                        "w_start"] <= peak["peak_end"] + 1):
                    continue
                # otherwise end peak region and append to output list
                peak_list.append(deepcopy(peak))
                peak["peak_start"] = None
                peak["peak_end"] = None
        # append peak if last window in data frame was significant
        if peak["peak_start"] is not None:
            peak_list.append(deepcopy(peak))
        peak_df = pd.DataFrame(peak_list)
        if peak_df.shape[0] > 1:
            peak_df = peak_df.loc[:, ["peak_start", "peak_end"]]
        return peak_df

    def _get_overlap(self, peak_start, peak_end, feature_start, feature_end):
        return max(0,
                   min(peak_end, feature_end)
                   - max(peak_start, feature_start) + 1)
