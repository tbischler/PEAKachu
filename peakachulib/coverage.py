import pysam
import numpy as np
from os import makedirs
from os.path import exists
from concurrent import futures
from time import time
import json
from peakachulib.wiggle import WiggleWriter


class Coverage(object):

    def __init__(self, paired_end, max_insert_size):
        self._paired_end = paired_end
        self._max_insert_size = max_insert_size
        self._coverages = {}

    def calc_coverages(self, bam_file):
        bam_fh = pysam.Samfile(bam_file, "rb")
        for replicon, length in sorted(zip(bam_fh.references, bam_fh.lengths)):
            self._init_coverage_list(length)
            if self._paired_end:
                self._cache_read2(replicon, bam_fh)
            self._calc_coverage(replicon, bam_fh)
            yield(replicon, self._coverages)
        bam_fh.close()

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


def generate_normalized_wiggle_files(project_folder, max_proc):
    parameter_dict = _read_parameters(project_folder)
    # create normalized coverage folder if it does not exist
    wiggle_folder = "{}/normalized_coverage".format(project_folder)
    if not exists(wiggle_folder):
        makedirs(wiggle_folder)
    # Generate coverage files in parallel
    print("** Generating normalized coverage files for {} libraries...".format(
          len(parameter_dict["libraries"])), flush=True)
    t_start = time()
    with futures.ProcessPoolExecutor(
            max_workers=max_proc) as executor:
        future_to_lib_name = {
            executor.submit(
                _generate_normalized_wiggle_file_for_lib, lib_name,
                lib["bam_file"], parameter_dict["paired_end"],
                parameter_dict["max_insert_size"], lib["size_factor"],
                wiggle_folder): lib_name for lib_name, lib
            in parameter_dict["libraries"].items()}
    for future in futures.as_completed(future_to_lib_name):
        lib_name = future_to_lib_name[future]
        print("* Coverage files for library {} generated.".format(lib_name),
              flush=True)
    t_end = time()
    print("Coverage file generation finished in {} seconds.".format(
        t_end-t_start), flush=True)


def _read_parameters(project_folder):
    with open("{}/parameters.json".format(project_folder),
              'r') as parameter_fh:
        parameter_dict = json.load(parameter_fh)
    return parameter_dict


def _generate_normalized_wiggle_file_for_lib(lib_name, bam_file, paired_end,
                                             max_insert_size, size_factor,
                                             wiggle_folder):
    """Perform the coverage calculation for a given library."""
    strand_dict = {"+": "forward", "-": "reverse"}
    wiggle_writers = dict([(strand, WiggleWriter(
        "{}_{}".format(lib_name, strand),
        open("{}/{}_div_by_{}_{}.wig".format(
            wiggle_folder, lib_name, size_factor, strand),
            "w"))) for strand in strand_dict.values()])
    coverage = Coverage(paired_end, max_insert_size)
    for replicon, coverages in coverage.calc_coverages(bam_file):
        for strand in strand_dict:
            if strand == "-":
                factor = -1.0/size_factor
            else:
                factor = 1.0/size_factor
            try:
                wiggle_writers[strand_dict[
                    strand]].write_replicons_coverages(
                        replicon, coverages[strand], factor=factor)
            except Exception as exc:
                print("Library {}, replicon {}, {} strand generated an "
                      "exception during coverage file generation: {}".format(
                          lib_name, replicon, strand, exc), flush=True)
    for strand in strand_dict.values():
        wiggle_writers[strand].close_file()
