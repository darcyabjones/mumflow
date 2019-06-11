#!/usr/bin/env python3

import argparse
import sys


class BedGraph(object):

    def __init__(
        self,
        seqid,
        start,
        end,
        values,
        repeat=None,
        repeat_prop=0.3
    ):
        # assert isinstance(seqid, str)
        # assert isinstance(start, int)
        # assert isinstance(end, int)
        # assert all_element_isinstance(values, int)

        self.seqid = seqid
        self.start = start
        self.end = end
        self.values = values

        self.length = abs(end - start)

        self.repeat_prop = repeat_prop
        if repeat is None:
            self.repeat = self.is_repeat(repeat_prop)
        else:
            self.repeat = repeat
        return

    def is_repeat(self, prop):
        """ Returns true if many coverage values are greater that 1. """

        from statistics import mean
        return mean(map(lambda x: x > 1, self.values)) > prop

    def flatten_repeat(self):
        """ Converts repeat coverage to simple presence absence """

        return self.__class__(
            self.seqid,
            self.start,
            self.end,
            [int(v >= 1) for v in self.values],
            repeat=self.repeat,
            repeat_prop=self.repeat_prop,
        )

    def all_present(self):
        """ Checks if all coverages >= 1."""
        return all(map(lambda v: v >= 1, self.values))

    def merge_with(self, other):
        """ Combine two BedGraph objects.

        Assumes that both have the same seqid, repeat, and values.
        """

        return self.__class__(
            self.seqid,
            self.start,
            other.end,
            self.values,
            repeat=self.repeat,
            repeat_prop=self.repeat_prop
        )

    @classmethod
    def from_line(cls, line, **kwargs):
        """ Parse a single bedgraph file line into an object. """

        sline = line.strip().split('\t')

        seqid = sline[0]
        start = int(sline[1])
        end = int(sline[2])

        values = list(map(int, sline[3:]))
        return cls(seqid, start, end, values, **kwargs)

    @classmethod
    def from_handle(cls, handle, comment="#", **kwargs):
        """ Parse bedgraph lines from a file handle. """

        for line in handle:
            if line.startswith(comment):
                continue
            yield cls.from_line(line, **kwargs)
        return

    def __str__(self):
        template = "{}\t{}\t{}\t{}\n"
        values = "\t".join(map(str, self.values))
        return template.format(self.seqid, self.start, self.end, values)


def all_element_isinstance(obj, typ):
    """ Reports if all elements of a collection is a type. """

    return all(isinstance(o, typ) for o in obj)


def are_adjacent(bg1, bg2):
    """ Tells us if we can merge two blocks. """

    same_chr = bg1.seqid == bg2.seqid
    same_profile = bg1.values == bg2.values
    either_is_repeat = bg1.repeat or bg2.repeat
    return same_chr and same_profile and not either_is_repeat


def flanks_are_same(bg_left, bg_mid, bg_right, tolerance):
    """ Tells us if we can merge three blocks based on flanks. """

    same_chr = bg_left.seqid == bg_mid.seqid == bg_right.seqid
    is_short = bg_mid.length <= tolerance
    same_boundaries = bg_left.values == bg_right.values
    boundaries_are_repeats = bg_left.repeat or bg_right.repeat

    return (same_chr and is_short and same_boundaries
            and not boundaries_are_repeats)


def should_report(bg, length):
    """ Decide if we should report this block to the results. """

    long_enough = bg.length > length
    return long_enough and not bg.repeat and not bg.all_present()


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Finds blocks in a reference genome that have presence-absence
        variations.

        Takes as input a bedgraph file containing alignment coverage, for
        arbitrarily many samples. It simplifies coverage (optionally filtering
        out repeat regions), merges adjacent identical blocks and
        short blocks flanked by identical blocks, and reports all with PAVs
        over a minimum length.

        For PAVs below about 50 bp I'd suggest a traditional read-alignment
        based pipeline.

        Assumes that the input coverage bedfile covers the whole genome with no
        gaps between adjacent lines or overlap between blocks.

        Also assumes that the input is sorted.
        If it isn't, the unix command `sort -k1,1 -k2,2n my.bedgraph` will do.
        """
    )

    parser.add_argument(
        "-i", "--infile",
        required=True,
        type=argparse.FileType('r'),
        help="Input bedgraph file. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output fasta file path. Default stdout.",
    )

    parser.add_argument(
        "-n", "--no-header",
        dest="header",
        action="store_false",
        default=True,
        help="Use if the input file doesn't have a header line.",
    )

    parser.add_argument(
        "-t", "--tol",
        default=20,
        type=int,
        help="Tolerance in bp to merge aligned blocks.",
    )

    parser.add_argument(
        "-m", "--min-length",
        default=50,
        type=int,
        help="Minimum length of blocks to report.",
    )

    parser.add_argument(
        "-r", "--proportion-repeats",
        dest="proportion_repeats",
        default=1.0,
        type=float,
        help=(
            "Filter out blocks where coverage is > 1 in greater than this "
            "proportion of samples. Float in range [0.0, 1.0]. "
            "Default 1.0 (i.e. no filtering)"
        ),
    )

    return parser.parse_args(args)


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    # Skip over the header line and print it.
    if args.header:
        header = next(args.infile)
        args.outfile.write(header)

    prev = None
    current = None

    bg_iter = BedGraph.from_handle(
        args.infile,
        repeat_prop=args.proportion_repeats
    )

    for next_ in bg_iter:
        next_ = next_.flatten_repeat()

        if prev is None:
            prev = next_
            continue
        if current is None:
            current = next_
            continue

        # Merge adjacent blocks if necessary
        if are_adjacent(prev, current):
            prev = prev.merge_with(current)
            current = None
            continue

        # Merge short deviations flanked by same.
        if flanks_are_same(prev, current, next_, args.tol):
            prev = prev.merge_with(next_)
            current = None
        else:
            if should_report(prev, args.min_length):
                args.outfile.write(str(prev))
            prev = current
            current = next_

    # After the loop we need to handle the trailing blocks.

    # Merge adjacent blocks if necessary
    if current is not None and prev is not None:
        if are_adjacent(prev, current):
            prev = prev.merge_with(current)
            current = None

    if (prev is not None) and should_report(prev, args.min_length):
        args.outfile.write(str(prev))
    if (current is not None) and should_report(current, args.min_length):
        args.outfile.write(str(current))

    return


if __name__ == "__main__":
    main()
