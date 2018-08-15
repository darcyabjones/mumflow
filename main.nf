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
    file "${delta.baseName}.unqry" optional true into unqrys


    """
    dnadiff --prefix=${delta.baseName} --delta ${delta}
    """
}


process snp2vcf {

    //container "python:3.6-alpine"

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
    container "quay.io/biocontainers/samtools:1.9--h46bd0b3_0"
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
    bedtools unionbedg -header -names *.bedgraph -i *.bedgraph > combined.bedgraph
    """
}

percentages = [1, 2, 3, 4, 5, 7, 10, 20, 30, 40, 50, 60, 70, 80, 90, 93, 95, 96, 97, 98, 99]

process compareCoverage {
    publishDir "compare_coverages"

    input:
    file bg from combinedCoverage
    val perc from percentages

    output:
    file "core_${perc}pc_regions.bedgraph" into coreRegions
    file "duplicate_${perc}pc_regions.bedgraph" into duplicateRegions
    file "single_${perc}pc_regions.bedgraph" into singleRegions

    """
    /usr/bin/env python3

    import os
    import pandas as pd

    table = pd.read_table("${bg}")
    table.columns = [os.path.split(os.path.splitext(c)[0])[1] for c in table.columns]

    # Cheating a little bit here
    table.drop(columns=["15FG031_contigs", "15FG039_contigs", "15FG046_contigs", "203FG217_contigs"], inplace=True) 
    pc = ${perc} / 100

    table2 = table[(table.iloc[:, 3:] > 0).mean(axis=1) > pc]
    table2.to_csv("core_${perc}pc_regions.bedgraph", index=False, sep="\t")

    table3 = table2[(table2.iloc[:, 3:] > 1).mean(axis=1) > pc]
    table3.to_csv("duplicate_${perc}pc_regions.bedgraph", index=False, sep="\t")

    table4 = table2[(table2.iloc[:, 3:] == 1).mean(axis=1) > pc]
    table4.to_csv("single_${perc}pc_regions.bedgraph", index=False, sep="\t")
    """
}


covFiles = coreRegions.spread("core").concat( duplicateRegions.spread("duplicate"), singleRegions.spread("single") )

covFiles.subscribe { println it }
