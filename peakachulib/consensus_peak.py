from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from peakachulib.wiggle import WiggleParser


class ConsensusPeakGenerator(object):

    def __init__(self, project_folder, consensus_length):
        self._project_folder = project_folder
        self._consensus_length = consensus_length
        self._replicon_peak_dict = defaultdict(lambda: defaultdict(set))

    def _store_peaks(self):
        peak_table_folder = "{}/peak_tables".format(self._project_folder)
        peak_files = [join(peak_table_folder, f) for f in listdir(
            peak_table_folder) if isfile(join(peak_table_folder, f))]
        for peak_file in peak_files:
            peak_df = pd.read_table(peak_file, sep='\t')
            for peak in peak_df.to_dict("records"):
                self._replicon_peak_dict[peak["replicon"]][
                    peak["peak_strand"]].add(
                        (peak["peak_start"], peak["peak_end"]))

    def _get_peak_coverage(self):
        norm_coverage_folder = "{}/normalized_coverage".format(
            self._project_folder)
        coverage_files = [join(norm_coverage_folder, f) for f in listdir(
            norm_coverage_folder) if isfile(join(norm_coverage_folder, f))]
        wiggle_parser = WiggleParser()
        cons_value_dict = defaultdict(dict)
        for coverage_file in coverage_files:
            cons_values = np.zeros(self._consensus_length)
            with open(coverage_file, 'r') as cov_fh:
                for wiggle_entry in wiggle_parser.entries(cov_fh):
                    lib_name_and_strand = wiggle_entry.track_name
                    lib_name = '_'.join(lib_name_and_strand.split('_')[:-1])
                    lib_strand = '+' if lib_name_and_strand.split(
                        '_')[-1] == "forward" else '-'
                    replicon = wiggle_entry.replicon
                    pos_value_pairs = dict(wiggle_entry.pos_value_pairs)
                    self._get_coverage_for_replicon_peaks(
                        replicon, lib_strand, pos_value_pairs, cons_values)
            cons_value_dict[lib_name][lib_strand] = cons_values
        # combine strands
        comb_cons_value_dict = {}
        for lib in cons_value_dict:
            comb_cons_value_dict[lib] = np.zeros(self._consensus_length)
            for strand in cons_value_dict[lib]:
                comb_cons_value_dict[lib] += cons_value_dict[lib][strand]
        return comb_cons_value_dict

    def _get_coverage_for_replicon_peaks(self, replicon, lib_strand,
                                         pos_value_pairs, cons_values):
        for peak in self._replicon_peak_dict[replicon][lib_strand]:
            value_list = []
            peak_center = int((peak[0] + peak[1]) / 2)
            if self._consensus_length % 2 == 0:
                cons_start = peak_center - (
                    int(self._consensus_length / 2) - 1)
                cons_end = peak_center + int(self._consensus_length / 2)
            else:
                cons_start = peak_center - (
                    int(self._consensus_length / 2) - 1)
                cons_end = peak_center + (int(self._consensus_length / 2) + 1)
            for pos in range(cons_start, cons_end + 1):
                value_list.append(abs(pos_value_pairs.get(pos, 0.0)))
            cons_values += np.array(value_list)

    def plot_consensus_peak(self):
        self._store_peaks()
        comb_cons_value_dict = self._get_peak_coverage()
        df = pd.DataFrame(comb_cons_value_dict, columns=sorted(
            comb_cons_value_dict))
        ax = df.plot(title="Consensus peak per library")
        ax.set_xlabel("Nucleotide position")
        ax.set_ylabel("Relative expression")
        plt.savefig("{}/plots/consensus_peaks.pdf".format(
            self._project_folder))
