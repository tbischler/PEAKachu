import pysam
from collections import defaultdict


class BamToBed(object):

    def __init__(self, paired_end, max_insert_size):
        self._paired_end = paired_end
        self._max_insert_size = max_insert_size

    def generate_bed_format(self, bam_file):
        self._bam_fh = pysam.AlignmentFile(bam_file, "rb")
        for replicon in self._bam_fh.references:
            self._bed_read_dict = defaultdict(int)
            if self._paired_end:
                self._cache_read2(replicon)
            self._calc_bed_reads(replicon)
            yield(replicon, self._bed_read_dict)
        self._bam_fh.close()

    def _calc_bed_reads(self, replicon):
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
        if aligned_read.is_reverse:
            self._bed_read_dict["{0},{1},{2}".format(start, end, '-')] += 1
        else:
            self._bed_read_dict["{0},{1},{2}".format(start, end, '+')] += 1

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
            mate = self._read2_dict[aligned_read.query_name.split()[0]]
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
        if aligned_read.is_reverse:
            self._bed_read_dict["{0},{1},{2}".format(start, end, '-')] += 1
        else:
            self._bed_read_dict["{0},{1},{2}".format(start, end, '+')] += 1

    def _cache_read2(self, replicon):
        self._read2_dict = {}
        for aligned_read in self._bam_fh.fetch(replicon):
            if not aligned_read.is_read2:
                continue
            if not (aligned_read.is_paired and aligned_read.is_proper_pair):
                continue
            if aligned_read.mate_is_unmapped:
                continue
            self._read2_dict[aligned_read.query_name.split()[0]] = aligned_read
        self._bam_fh.seek(0)
