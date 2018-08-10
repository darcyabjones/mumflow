#!/usr/bin/env python3

import csv
import operator
import argparse
import os
import sys
import copy


def parse_fasta(handle):
    """ Simple fasta parser, returns dict of strings """

    output = dict()

    current_id = None
    current_seq = []
    for line in handle:
        if line.startswith(">"):
            if current_id is not None:
                output[current_id] = ''.join(current_seq)

            current_id = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())

    if current_id is not None:
        output[current_id] = ''.join(current_seq)

    return output




#with open("${snp}") as handle:
def read_snp(handle):
    """ Reads a mummer snp file """

    dialect = csv.Sniffer().sniff(handle.read(1024))
    handle.seek(0)

    reader = csv.DictReader(
        handle,
        fieldnames=[
            "pos_ref",
            "allele_ref",
            "allele_query",
            "pos_query",
            "buff", # Distance from this SNP to the nearest mismatch
            "dist", # Distance from this SNP to the nearest sequence end.
            "len_ref",
            "len_query",
            "strand_ref",
            "strand_query",
            "id_ref",
            "id_query",
        ],
        dialect=dialect,
    )

    return reader


def convert_snp(line, sample_name):
    output = {}
    output["#CHROM"] = line["id_ref"]
    output["POS"] = int(line["pos_ref"])
    output["ID"] = "."
    output["REF"] = line["allele_ref"] # We handle INDELS later
    output["ALT"] = line["allele_query"]
    output["QUAL"] = "."
    output["FILTER"] = "PASS"
    output["INFO"] = "."
    output["FORMAT"] = "GT"
    output[str(sample_name)] = "1" # This may need to change for diploid orgs.
    return output


def snp_to_vcf(snps, sample_name):
    """ Simple wrapper around a generator to convert snps """
    return (convert_snp(l, sample_name) for l in snps)


def fold_indels(lines, ref):
    """ """

    # These should already be sorted but just to be sure.
    lines = sorted(lines, key=operator.itemgetter("#CHROM", "POS"))

    output = []

    last_line = None
    block = []

    for line in lines:
        if should_extend(block, line):
            block.append(line)
        else:
            if len(block) > 0:
                output.append(join_block(block, ref))
                block = []

            if snp_type(line) == "snp":
                output.append(line)
            else:
                block.append(line)

    if len(block) > 0:
        output.append(join_block(block, ref))

    return output


def join_block(block, ref):
    chrom = block[0]["#CHROM"]
    first_pos = block[0]["POS"]
    ref_base = ref[chrom][first_pos - 2]

    ref_allele = [ref_base] + [b["REF"] for b in block if b["REF"] != '.']
    alt_allele = [ref_base] + [b["ALT"] for b in block if b["ALT"] != '.']

    new_line = copy.copy(block[0])
    new_line["POS"] -= 1
    new_line["REF"] = ''.join(ref_allele)
    new_line["ALT"] = ''.join(alt_allele)
    return new_line


def snp_type(line):
    if line["REF"] == ".":
        type_ = "insertion"
    elif line["ALT"] == ".":
        type_ = "deletion"
    else:
        type_ = "snp"
    return type_


def should_extend(block, line):
    if len(block) == 0:
        return False

    last_chrom = block[-1]["#CHROM"]
    last_pos = block[-1]["POS"]
    if ((snp_type(line) == snp_type(block[-1]))
            and (last_pos == line["POS"] - 1)
            and (last_chrom == line["#CHROM"])):
        return True
    else:
        return False


#with open("${snp.baseName}.vcf") as handle:
def write_vcf(lines, handle, sample_name):
    """ Writes a vcf file from a list of dictionaries. """

    handle.write("##fileformat=VCFv4.1\n")
    handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fieldnames = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                  'INFO', 'FORMAT', str(sample_name)]

    writer = csv.DictWriter(handle, fieldnames=fieldnames, dialect="excel-tab")
    writer.writeheader()

    for line in lines:
        writer.writerow(line)

    return


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Project description."
    )

    parser.add_argument(
        "-s", "--snp",
        required=True,
        type=argparse.FileType('r'),
        help="SNP file. Use '-' for stdin.",
    )

    parser.add_argument(
        "-r", "--reference",
        required=True,
        type=argparse.FileType('r'),
        help="FASTA file. Use '-' for stdin.",
    )

    parser.add_argument(
        "-n", "--name",
        required=True,
        type=str,
        help="name for the snp in vcf",
    )

    parser.add_argument(
        "-o", "--output",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output fasta file path. Default stdout.",
    )

    return parser.parse_args(args)


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    ref = parse_fasta(args.reference)
    snps = snp_to_vcf(read_snp(args.snp), args.name)

    snps = fold_indels(snps, ref)

    write_vcf(snps, args.output, args.name)
    return


if __name__ == "__main__":
    main()
