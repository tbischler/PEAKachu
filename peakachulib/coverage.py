import pysam
import numpy as np


class Coverage(object):

    def __init__(self, paired_end, max_insert_size):
        self._paired_end = paired_end
        self._max_insert_size = max_insert_size
        self._coverages = {}

    def calc_coverages(self, bam_file):
        self._bam_file = bam_file
        bam_fh = pysam.Samfile(self._bam_file, "rb")
        for replicon, length in sorted(zip(bam_fh.references, bam_fh.lengths)):
            self._init_coverage_list(length)
            if self._paired_end:
                self._cache_read2(replicon, bam_fh)
            self._calc_coverage(replicon, bam_fh)
            yield(replicon, self._coverages)

    def _init_coverage_list(self, length):
        for strand in ["+", "-"]:
            self._coverages[strand] = np.array([0.0] * length)

    def _calc_coverage(self, replicon, bam_fh):
        for aligned_read in bam_fh.fetch(replicon):
            if self._paired_end:
                self._add_paired_end_coverage(aligned_read)
            else:
                self._add_single_end_coverage(aligned_read)

    def _add_single_end_coverage(self, aligned_read):
        # Note: No translation from SAMParser coordinates to python
        # list coorindates is needed.
        start = aligned_read.pos
        end = aligned_read.aend
        if aligned_read.is_reverse:
            self._coverages["-"][start:end] += 1.0
        else:
            self._coverages["+"][start:end] += 1.0

    def _add_paired_end_coverage(self, aligned_read):
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
            mate = self._read2_dict[aligned_read.query_name.split()[0]]
            # print("Mate found: %s" % mate)
        except KeyError:
            # Invalid mate (usually post-filtered)
            print("Mate not found for read %s" % aligned_read)
            return
        if not mate.is_read2:
            return
        # Note: No translation from SAMParser coordinates to python
        # list coorindates is needed.
        start = aligned_read.pos
        end = mate.aend
        if end <= start:
            return
        if end - start > self._max_insert_size:
            return
        if aligned_read.is_reverse:
            self._coverages["-"][start:end] += 1.0
        else:
            self._coverages["+"][start:end] += 1.0

    def _cache_read2(self, replicon, bam_fh):
        self._read2_dict = {}
        for aligned_read in bam_fh.fetch(replicon):
            if not aligned_read.is_read2:
                continue
            if not (aligned_read.is_paired and aligned_read.is_proper_pair):
                continue
            if aligned_read.mate_is_unmapped:
                continue
            self._read2_dict[aligned_read.query_name.split()[0]] = aligned_read
        bam_fh.seek(0)
