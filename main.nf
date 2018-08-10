#!/usr/bin/env nextflow

params.fastas = "$baseDir/data/genomes/*.fasta"

reference_file = file(params.reference)
fasta_files = Channel.fromPath( params.fastas )

process align {
    container "quay.io/biocontainers/mummer4:4.0.0beta2--pl526hfc679d8_2"
    publishDir "alignments"

    input:
    file reference from reference_file
    file query from fasta_files

    output:
    file "${query.baseName}.delta" into deltas

    """
    nucmer --prefix=${query.baseName} ${reference} ${query}
    """
}

process dnadiff {
    container "quay.io/biocontainers/mummer4:4.0.0beta2--pl526hfc679d8_2"

    publishDir "diffs"

    input:
    file delta from deltas

    output:
    file "${delta.baseName}.snps" into snps
    file "${delta.baseName}.1coords" into onecoords
    file "${delta.baseName}.1delta" into onedeltas
    file "${delta.baseName}.mcoords" into mcoords
    file "${delta.baseName}.mdelta" into mdeltas
    file "${delta.baseName}.qdiff" into qdiff
    file "${delta.baseName}.rdiff" into rdiff
    file "${delta.baseName}.report" into dnadiffReports
    file "${delta.baseName}.unqry" into unqrys


    """
    dnadiff --prefix=${delta.baseName} --delta ${delta}
    """
}


process snp2vcf {

    container "python:3.6-alpine"

    publishDir "vcfs"

    input:
    file snp from snps
    file reference from reference_file

    output:
    file "${snp.baseName}.vcf" into vcfs

    """
    snps2vcf.py --snp ${snp} --reference ${reference} --name ${snp.baseName} --output ${snp.baseName}.vcf
    """
}


process coordsToBED {
    publishDir "bed_alignments"

    input:
    file coord from mcoords

    output:
    file "${coord.baseName}.bed" into coordsBEDs

    """
    tail -n +5 ${coord} \
    | awk '{ printf "%s\\t%s\\t%s\\t\\.\\t+\\t%s\\tfiltered\\n", \$12, \$1-1, \$2, \$7 }' \
    | sort -k1,1 -k2,2n \
    > ${coord.baseName}.bed
    """
}


process genomeIndex {
    container "quay.io/biocontainers/samtools"
    input:
    file reference from reference_file

    output:
    file "${reference}.fai" into referenceIndex

    """
    samtools faidx ${reference} > ${reference}.fai
    """
}


process bedCoverage {
    container "quay.io/biocontainers/bedtools:2.27.1--1"
    publishDir "bedgraph_coverages"

    input:
    file bed from coordsBEDs
    file index from referenceIndex

    output:
    file "${bed.baseName}.bedgraph" into coverageBEDs

    """
    bedtools genomecov -bga -g ${index} -i ${bed} > ${bed.baseName}.bedgraph
    """
}


process combinedBEDCoverage {
    container "quay.io/biocontainers/bedtools:2.27.1--1"
    publishDir "bedgraph_coverages"

    input:
    file "*.bedgraph" from coverageBEDs.collect()

    output:
    file "combined.bedgraph" into combinedCoverage

    """
    bedtools unionbedg -i *.bedgraph > combined.bedgraph
    """
}


