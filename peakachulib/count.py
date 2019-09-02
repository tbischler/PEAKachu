import pysam
import numpy as np
from collections import defaultdict
from peakachulib.intersection import Intersecter, Interval


class ReadCounter(object):

    def __init__(self, paired_end, max_insert_size, bam_file):
        self._paired_end = paired_end
        self._max_insert_size = max_insert_size
        self._bam_fh = pysam.AlignmentFile(bam_file, "rb")

    def count_reads_for_windows(self, replicon, strand, window_list):
        self._interval_tree = Intersecter()
        self._counts = np.array([0] * len(window_list))
        for ind, window in enumerate(window_list):
            self._interval_tree.add_interval(
                Interval(window[0],
                         window[1],
                         value=ind,
                         strand=strand))
        if self._paired_end:
            self._cache_read2(replicon)
        self._count_reads(replicon)
        return self._counts

    def count_reads_for_peaks(self, replicon, peak_list):
        self._interval_tree = Intersecter()
        self._counts = np.array([0] * len(peak_list))
        for ind, peak in enumerate(peak_list):
            self._interval_tree.add_interval(
                Interval(peak["peak_start"]-1,
                         peak["peak_end"],
                         value=ind,
                         strand=peak["peak_strand"]))
        if self._paired_end:
            self._cache_read2(replicon)
        self._count_reads(replicon)
        return self._counts

    def close_bam(self):
        self._bam_fh.close()

    def _count_reads(self, replicon):
        for aligned_read in self._bam_fh.fetch(replicon):
            if self._paired_end:
                self._add_paired_end_read(aligned_read)
            else:
                self._add_single_end_read(aligned_read)

    def _add_single_end_read(self, aligned_read):
        # Note: No translation from SAMParser coordinates to python
        # list coorindates is needed.
        start = aligned_read.reference_start
        end = aligned_read.reference_end
        overlap = self._interval_tree.find(start, end)
        for interval in overlap:
            if aligned_read.is_reverse and interval.strand == "-":
                self._counts[interval.value] += 1
            elif not aligned_read.is_reverse and interval.strand == "+":
                self._counts[interval.value] += 1

    def _add_paired_end_read(self, aligned_read):
        if not aligned_read.is_read1:
            return
        if not (aligned_read.is_paired and aligned_read.is_proper_pair):
            return
        if aligned_read.mate_is_unmapped:
            return
        if aligned_read.is_reverse == aligned_read.mate_is_reverse:
            return
        if not (aligned_read.reference_id == aligned_read.next_reference_id):
            return
        try:
            mate = self._read2_dict[aligned_read.query_name.split()[0]][
                aligned_read.next_reference_start]
        except KeyError:
            # Invalid mate (usually post-filtered)
            print("Mate not found for read {}".format(aligned_read))
            return
        if not mate.is_read2:
            return
        # Note: No translation from SAMParser coordinates to python
        # list coorindates is needed.
        if aligned_read.is_reverse:
            start = mate.reference_start
            end = aligned_read.reference_end
        else:
            start = aligned_read.reference_start
            end = mate.reference_end
        if end <= start:
            return
        if end - start > self._max_insert_size:
            return
        overlap = self._interval_tree.find(start, end)
        for interval in overlap:
            if aligned_read.is_reverse and interval.strand == "-":
                self._counts[interval.value] += 1
            elif not aligned_read.is_reverse and interval.strand == "+":
                self._counts[interval.value] += 1

    def _cache_read2(self, replicon):
        self._read2_dict = defaultdict(dict)
        for aligned_read in self._bam_fh.fetch(replicon):
            if not aligned_read.is_read2:
                continue
            if not (aligned_read.is_paired and aligned_read.is_proper_pair):
                continue
            if aligned_read.mate_is_unmapped:
                continue
            self._read2_dict[aligned_read.query_name.split()[0]][
                aligned_read.reference_start] = aligned_read
        self._bam_fh.seek(0)
