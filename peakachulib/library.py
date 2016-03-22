from os.path import basename, splitext
import pandas as pd
import numpy as np
from copy import deepcopy
from peakachulib.coverage import Coverage
from peakachulib.bam_to_bed import BamToBed


class Library(object):
    '''
    This class reads the alignment file for a library and counts and stores
    the reads mapping to different annotations
    '''
    def __init__(self, paired_end, max_insert_size, bam_file, replicon_dict):
        self._paired_end = paired_end
        self.bam_file = bam_file
        self._max_insert_size = max_insert_size
        self.lib_name = splitext(basename(bam_file))[0]
        self.replicon_dict = replicon_dict
        self._coverage_calculated = False

    def _calc_coverage(self):
        coverage = Coverage(self._paired_end, self._max_insert_size)
        for replicon, coverages in coverage.calc_coverages(self.bam_file):
            self.replicon_dict[replicon]["coverages"] = deepcopy(coverages)

    def _calc_window_expr(self):
        for replicon in self.replicon_dict:
            window_count_list = []
            for window in self.replicon_dict[replicon]["window_list"]:
                window_count_list.append(self._count_reads_per_window(
                    replicon, window[0], window[1]))
            self.replicon_dict[replicon]['window_counts'] = pd.DataFrame(
                [[window_expr[strand] for strand in ['+', '-']]
                 for window_expr in window_count_list], columns=['+', '-'])

    def _calc_peak_expr(self):
        for replicon in self.replicon_dict:
            peak_count_list = []
            for peak in self.replicon_dict[replicon]["peak_df"].to_dict(
                    'records'):
                peak_count_list.append(self._count_reads_per_peak(
                    replicon, peak["peak_start"]-1, peak["peak_end"],
                    peak["peak_strand"]))
            self.replicon_dict[replicon]["peak_counts"] = peak_count_list

    def _count_reads_per_window(self, replicon, start, end):
        window_expr = {}
        for strand in ["+", "-"]:
            window_expr[strand] = np.max(self.replicon_dict[replicon][
                "coverages"][strand][start:end])
        return window_expr

    def _count_reads_per_peak(self, replicon, start, end, strand):
        return np.max(self.replicon_dict[replicon]["coverages"][strand]
                       [start:end])

    def count_reads_for_windows(self):
        if not self._coverage_calculated:
            self._calc_coverage()
            self._coverage_calculated = True
        self._calc_window_expr()
        return self.replicon_dict  # it seems that a copy is returned!

    def merge_reads(self):
        bam_to_bed = BamToBed(self._paired_end, self._max_insert_size)
        for replicon, reads in bam_to_bed.generate_bed_format(self.bam_file):
            self.replicon_dict[replicon]["reads"] = pd.Series(reads)
        return self.replicon_dict  # it seems that a copy is returned!

    def count_reads_for_peaks(self):
        if not self._coverage_calculated:
            self._calc_coverage()
            self._coverage_calculated = True
        self._calc_peak_expr()
