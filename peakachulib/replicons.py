import sys
import pysam
from os import listdir
from os.path import isfile, join
from BCBio import GFF

class Replicons(object):
    def __init__(self, control_libs, tagged_libs, gff_folder, features,
                 sub_features):
        self._libs = control_libs + tagged_libs
        self._gff_files = [ join(gff_folder, f) for f in listdir(
                gff_folder) if isfile(join(gff_folder, f)) ]
        self._limit_info = dict(
                gff_type = features)
        self._sub_features = sub_features
        self.replicon_dict = {}
    
    def init_replicons(self):
        for bam_file in self._libs:
            bam_fh = pysam.Samfile(bam_file, "rb")
            replicon_dict = dict([[replicon, {'seq_start_pos': 0,
                    'seq_end_pos': length}] for replicon, length in
                    zip(bam_fh.references, bam_fh.lengths)])
            if self.replicon_dict and not (self.replicon_dict == replicon_dict):
                sys.stderr.write("Mismatch in replicon content of bam files!\n")
                sys.exit(1)
            self.replicon_dict = replicon_dict
        if not self._gff_files:
            sys.stderr.write("No .gff file found in specified folder")
            #sys.exit(1)
        for gff_file in self._gff_files:
            self._store_annotations(gff_file)
        
        print("Peak detection will be conducted for the following sequence "
              "regions:", flush=True)
        for seq_id, replicon in sorted(self.replicon_dict.items()):
            print("%s: %s %s" % (
                    seq_id,
                    replicon['seq_start_pos'] + 1,
                    replicon['seq_end_pos']), flush=True)
            if not 'features' in replicon:
                replicon['features'] = []
    
    def _store_annotations(self, gff_file):
        with open(gff_file) as gff_fh:
            for rec in GFF.parse(gff_fh, limit_info=self._limit_info):
                if not rec.id in self.replicon_dict:
                    sys.stderr.write(
                            "Annotations for sequence ID %s skipped as sequence"
                            "is not present in alignment files!\n" % (rec.id))
                    continue
                if not 'features' in self.replicon_dict[rec.id]:
                    self.replicon_dict[rec.id]['features'] = []
                for feature in rec.features:
                    feature_entry = {
                            'type': feature.type,
                            'start': feature.location.start.position+1,
                            'end': feature.location.end.position,
                            'strand': '+' if feature.strand == 1 else '-',
                            'locus_tag': '/'.join(feature.qualifiers[
                                'locus_tag']) if ('locus_tag' in 
                                feature.qualifiers) else None,
                            'Name': '/'.join(feature.qualifiers[
                                'Name']) if ('Name' in 
                                feature.qualifiers) else None,
                            'product': '/'.join(feature.qualifiers[
                                'product']) if ('product' in 
                                feature.qualifiers) else None
                            }
                    for subfeature in feature.sub_features:
                        if not subfeature.type in self._sub_features:
                            continue
                        feature_entry['subfeature_type'] = subfeature.type
                        if feature_entry['product'] is not None:
                            continue
                        feature_entry['product'] = '/'.join(
                                subfeature.qualifiers['product']) if (
                                'product' in subfeature.qualifiers
                                ) else None
                    self.replicon_dict[rec.id]['features'].append(
                            feature_entry)
        
