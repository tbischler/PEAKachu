import csv

class WiggleParser(object):
    """

    Warning - this does not implement the full specification!

    """

    def entries(self, input_fh):
        track_name = None
        replicon = None
        span = None
        pos_value_pairs = []
        for line in input_fh:
            row = line[:-1].split()
            if len(row) == 0:
                continue
            if row[0].startswith("track"):
                track_name = self._track_name(row)
            elif row[0].startswith("variableStep"):
                if replicon:
                    prev_replicon = replicon
                    prev_span = span
                    prev_pos_value_pairs = pos_value_pairs
                    replicon = self._replicon(row)
                    span = None
                    pos_value_pairs = []
                    yield WiggleEntry(
                            track_name, prev_replicon, prev_span,
                            prev_pos_value_pairs)
                else:
                    replicon = self._replicon(row)
            else:
                pos_value_pairs.append([int(row[0]), float(row[1])])
        yield WiggleEntry(track_name, replicon, span, pos_value_pairs)

    def _replicon(self, row):
        return self._attrs_and_values(row)["chrom"]

    def _track_name(self, row):
        return self._attrs_and_values(row)["name"]

    def _attrs_and_values(self, row):
        attrs_and_values = {}
        for attr_and_value in row:
            if not "=" in attr_and_value:
                continue
            attr, value = attr_and_value.split("=")
            value = value.replace("\"", "")
            attrs_and_values[attr] = value
        return attrs_and_values


class WiggleEntry(object):

    def __init__(self, track_name, replicon, span, pos_value_pairs):
        self.track_name = track_name
        self.replicon = replicon
        self.span = span
        self.pos_value_pairs = pos_value_pairs


class WiggleWriter(object):

    def __init__(self, track_str, fh):
        self._fh = fh
        self._fh.write(("track type=wiggle_0 name=\"%s\"\n" % (track_str)))

    def write_replicons_coverages(
            self, replicon_str, coverages, discard_zeros=True, factor=1.0):
            self._fh.write("variableStep chrom=%s span=1\n" % (replicon_str))
            # Filter values of 0 and multiply the remaining ones by
            # the given factor. pos is increased by 1 as a translation
            # from a 0-based sysem (Python list) to a 1 based system
            # (wiggle) takes place.
            self._fh.write(
                "\n".join(["%s %s" % (pos + 1, coverage * factor)
                           for pos, coverage in
                           filter(lambda pos_and_cov: pos_and_cov[1] != 0.0,
                                  enumerate(coverages))]) + "\n")

    def close_file(self):
        self._fh.close()
