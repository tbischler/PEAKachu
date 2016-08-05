import sys
import matplotlib
import matplotlib.pyplot as plt
from peakachulib.replicons import Replicons
from peakachulib.window import WindowApproach
from peakachulib.adaptive import AdaptiveApproach
from peakachulib.coverage import generate_normalized_wiggle_files
from peakachulib.consensus_peak import ConsensusPeakGenerator
from time import time

matplotlib.use('Agg')
plt.style.use('seaborn-colorblind')
font = {'family': 'sans-serif', 'size': 7}
matplotlib.rc('font', **font)


class Controller(object):

    """Manage the actions of the subcommands.

    The Controller take care of providing the argumentes like path
    names and the parallel processing of tasks.

    """

    def __init__(self, args):
        """Create an instance."""
        self._args = args

    def window_approach(self):
        '''
        This function initiates the peak detection via the sliding window
        approach
        '''
        # Count libraries
        lib_count = len(self._args.exp_libs) + len(self._args.ctr_libs)
        assert(len(self._args.exp_libs) == len(self._args.ctr_libs))
        # Step size must be lower than or equal to window size
        assert(self._args.step_size <= self._args.window_size)
        # Select normalization method and initialize size factors
        size_factors = self._select_normalization_window(lib_count)
        # Initialize replicons
        replicons = self._init_replicons()
        # Initialize class object for the window approach
        window = WindowApproach(
                self._args.window_size,
                self._args.step_size,
                replicons.replicon_dict,
                self._args.max_proc,
                self._args.stat_test,
                self._args.norm_method,
                size_factors,
                self._args.het_p_val_threshold,
                self._args.rep_pair_p_val_threshold,
                self._args.padj_threshold,
                self._args.mad_multiplier,
                self._args.fc_cutoff,
                self._args.pairwise_replicates,
                self._args.output_folder)
        # Initialize libraries
        print("** Initializing libraries...", flush=True)
        t_start = time()
        window.init_libraries(
                self._args.paired_end,
                self._args.max_insert_size,
                self._args.ctr_libs,
                self._args.exp_libs)
        t_end = time()
        print("Finished library initialization in {} seconds.\n".format(
            t_end-t_start), flush=True)
        window.generate_window_counts()
        # Perform statistical test
        if self._args.stat_test == "gtest":
            print("** Performing G-test for all windows...", flush=True)
            t_start = time()
            self._perform_g_test_window(window, lib_count)
            t_end = time()
            print("G-test finished in {} seconds.\n".format(t_end-t_start),
                  flush=True)
        elif self._args.stat_test == "deseq":
            print("** Running DESeq2 for all windows...", flush=True)
            t_start = time()
            self._run_deseq2_window(window)
            t_end = time()
            print("DESeq2 finished in {} seconds.\n".format(t_end-t_start),
                  flush=True)
        # Merge windows into peaks
        print("** Merging windows to peaks and recalculating values...",
              flush=True)
        t_start = time()
        window.combine_peaks_and_recalculate_values()
        t_end = time()
        print("Peak generation finished in {} seconds.\n".format(
            t_end-t_start), flush=True)
        # Generate output files
        print("** Writing peak output files...", flush=True)
        t_start = time()
        window.write_output()
        t_end = time()
        print("Writing output files took {} seconds.\n".format(t_end-t_start),
              flush=True)

    def adaptive_approach(self):
        '''
        This function initiates the peak detection via the adaptive approach
        '''
        # Count libraries
        lib_count = len(self._args.exp_libs) + len(self._args.ctr_libs)
        # Initialize replicons
        replicons = self._init_replicons()
        # Initialize class object for the adaptive approach
        adaptive = AdaptiveApproach(
                replicons.replicon_dict,
                self._args.max_proc,
                self._args.padj_threshold,
                self._args.mad_multiplier,
                self._args.fc_cutoff,
                self._args.output_folder)
        # Initialize libraries
        print("** Initializing libraries...", flush=True)
        t_start = time()
        adaptive.init_libraries(
                self._args.paired_end,
                self._args.max_insert_size,
                self._args.ctr_libs,
                self._args.exp_libs)
        t_end = time()
        print("Finished library initialization in {} seconds.\n".format(
            t_end-t_start), flush=True)
        # Combine reads from experiment libraries in one BED file
        adaptive.generate_combined_bed_file()
        # Run blockbuster
        print("** Running blockbuster...", flush=True)
        t_start = time()
        adaptive.run_blockbuster()
        t_end = time()
        print("blockbuster finished in {} seconds.\n".format(t_end-t_start),
              flush=True)
        # Generate peaks from blockbuster output
        print("** Generating peaks from blockbuster output...",
              flush=True)
        t_start = time()
        adaptive.generate_peaks_from_blockbuster(
            self._args.min_cluster_expr_frac, self._args.min_block_overlap,
            self._args.min_max_block_expr_frac)
        t_end = time()
        print("Peak generation finished in {} seconds.\n".format(
            t_end-t_start), flush=True)
        # Calculate peak expression
        print("** Calculating peak expression...",
              flush=True)
        t_start = time()
        adaptive.calculate_peak_expression()
        t_end = time()
        print("Peak expression calculated in {} seconds.\n".format(
            t_end-t_start), flush=True)
        # Select normalization method and initialize size factors
        size_factors = self._select_normalization_adaptive(lib_count)
        # Calculate significant peaks
        self._calc_sig_peaks_adaptive(adaptive, size_factors)
        # Generate output files
        print("** Writing peak output files...", flush=True)
        t_start = time()
        adaptive.write_output()
        t_end = time()
        print("Writing output files took {} seconds.\n".format(t_end-t_start),
              flush=True)

    def coverage(self):
        '''
        This function generates normalized coverage files in wiggle format for
        each library and strand
        '''
        generate_normalized_wiggle_files(
            self._args.project_folder, self._args.max_proc)

    def consensus_peak(self):
        '''
        This function generates a consensus peak per library based on
        previously called peaks and normalized coverage files
        '''
        consensus_peak_generator = ConsensusPeakGenerator(
            self._args.project_folder, self._args.consensus_length)
        consensus_peak_generator.plot_consensus_peak()

    def _init_replicons(self):
        print("** Initializing replicons and reading annotations from .gff "
              "files if present...", flush=True)
        t_start = time()
        replicons = Replicons(self._args.ctr_libs, self._args.exp_libs,
                              self._args.gff_folder, self._args.features,
                              self._args.sub_features)
        replicons.init_replicons()
        t_end = time()
        print("Finished replicon initialization in {} seconds.\n".format(
            t_end-t_start), flush=True)
        return replicons

    def _select_normalization_window(self, lib_count):
        print("Selected normalization method is: {}".format(
            self._args.norm_method))
        if self._args.norm_method == "manual":
            if not (len(self._args.size_factors) == (lib_count)):
                sys.stderr.write(
                        "Normalization factors do not match library number!\n")
                sys.exit(1)
            size_factors = [value/max(self._args.size_factors)
                            for value in self._args.size_factors]
        elif self._args.norm_method == "none":
            if self._args.size_factors:
                sys.stderr.write("Specified size factors were ignored!\n")
            size_factors = [1.0] * (lib_count)
        else:
            if self._args.size_factors:
                sys.stderr.write("Specified size factors were ignored!\n")
            size_factors = None
        return size_factors

    def _select_normalization_adaptive(self, lib_count):
        print("Selected normalization method is: {}".format(
            self._args.norm_method))
        if self._args.norm_method == "manual":
            if not (len(self._args.size_factors) == (lib_count)):
                sys.stderr.write(
                        "Normalization factors do not match library number!\n")
                sys.exit(1)
            size_factors = [value/max(self._args.size_factors)
                            for value in self._args.size_factors]
        elif self._args.norm_method == "none" or not self._args.ctr_libs:
            if self._args.size_factors:
                sys.stderr.write("Specified size factors were ignored!\n")
            size_factors = [1.0] * (lib_count)
        else:
            if self._args.size_factors:
                sys.stderr.write("Specified size factors were ignored!\n")
            size_factors = None
        return size_factors

    def _calc_sig_peaks_adaptive(self, adaptive, size_factors):
        # Run DESeq2 only if >1 libraries are available for both, experiment
        # and control libraries
        if len(self._args.exp_libs) > 1 and len(self._args.ctr_libs) > 1:
            print("** Calculating peak significance with DESeq2...",
                  flush=True)
            t_start = time()
            adaptive.run_deseq2_analysis(size_factors,
                                         self._args.pairwise_replicates)
            t_end = time()
            print("DESeq2 finished in {} seconds.\n".format(t_end-t_start),
                  flush=True)
        # If at least one control is available use fold change and MAD to
        # define significant peaks
        elif self._args.ctr_libs:
            print("** Calculating peak significance for insufficient "
                  "replicates based on fold change and MAD cutoff...",
                  flush=True)
            t_start = time()
            adaptive.run_analysis_without_replicates(size_factors)
            t_end = time()
            print("Peak calculation finished in {} seconds.\n".format(
                t_end-t_start), flush=True)
        # If no controls are available return initial peaks based on MAD cutoff
        else:
            print("** Calculating peaks without control based on MAD "
                  "cutoff...", flush=True)
            t_start = time()
            adaptive.run_analysis_without_control(size_factors)
            t_end = time()
            print("Peak calculation finished in {} seconds.\n".format(
                t_end-t_start), flush=True)

    def _perform_g_test_window(self, window, lib_count):
        print("Testing differential expression via G-test")
        if lib_count > 2:
            print("* Running in replicate-mode...", flush=True)
            if self._args.pairwise_replicates:
                print("  Experiment and control libraries are treated as pairs"
                      " according to input order.")
            else:
                print("  All combinations of experiment and control library "
                      "pairs that include all libraries are tested and the "
                      "result with lowest significance is selected.")
            window.perform_g_test_with_repl_for_windows()
        else:
            print("* Running in without-replicate-mode...", flush=True)
            window.perform_g_test_without_repl_for_windows()

    def _run_deseq2_window(self, window):
        print("Testing differential expression via DESeq2")
        if self._args.pairwise_replicates:
            print("  Experiment and control libraries are treated as pairs"
                  " according to input order.")
        else:
            print("  All combinations of experiment and control library "
                  "pairs that include all libraries are tested and the "
                  "result with lowest significance is selected.")
        window.run_deseq2_for_windows()
