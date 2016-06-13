from os.path import basename, splitext
import pandas as pd
from copy import deepcopy
from peakachulib.bam_to_bed import BamToBed
from peakachulib.count import ReadCounter


class Library(object):
    '''
    This class reads the alignment file for a library and counts and stores
    the reads mapping to different annotations
    '''
    def __init__(self, paired_end, max_insert_size, bam_file, replicon_dict):
        self.paired_end = paired_end
        self.bam_file = bam_file
        self.max_insert_size = max_insert_size
        self.lib_name = splitext(basename(bam_file))[0]
        self.replicon_dict = replicon_dict

    def _calc_window_expr(self):
        read_counter = ReadCounter(self.paired_end, self.max_insert_size,
                                   self.bam_file)
        for replicon in self.replicon_dict:
            self.replicon_dict[replicon]['window_counts'] = pd.DataFrame()
            for strand in ['+', '-']:
                window_counts = read_counter.count_reads_for_windows(
                            replicon,
                            strand,
                            self.replicon_dict[replicon]["window_list"])
                self.replicon_dict[replicon]['window_counts'][
                    strand] = window_counts
        read_counter.close_bam()

    def _calc_peak_expr(self):
        read_counter = ReadCounter(self.paired_end, self.max_insert_size,
                                   self.bam_file)
        for replicon in self.replicon_dict:
            peak_counts = read_counter.count_reads_for_peaks(
                replicon,
                self.replicon_dict[replicon]["peak_df"].to_dict('records'))
            del self.replicon_dict[replicon]["peak_df"]
            self.replicon_dict[replicon]["peak_counts"] = peak_counts
        read_counter.close_bam()

    def count_reads_for_windows(self):
        self._calc_window_expr()

    def merge_reads(self):
        bam_to_bed = BamToBed(self.paired_end, self.max_insert_size)
        for replicon, reads in bam_to_bed.generate_bed_format(self.bam_file):
            self.replicon_dict[replicon]["reads"] = pd.Series(reads)
        return self.replicon_dict  # it seems that a copy is returned!

    def count_reads_for_peaks(self):
        self._calc_peak_expr()
