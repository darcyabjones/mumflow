#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.genomes = false
params.references = false


if ( params.references ) {
    references = Channel.fromPath(
        params.references,
        checkIfExists: true,
        type: "file"
    )
} else {
    log.info "Hey I need some reference genomes please!"
    exit 1
}


if ( params.genomes ) {
    genomes = Channel.fromPath(
        params.genomes,
        checkIfExists: true,
        type: "file"
    )
} else {
    log.info "Hey I need some genomes to align please!"
    exit 1
}


references.into {
    references4Cross;
    references4Snp2Vcf;
    references4GenomeIndex;
}


pairs = references4Cross
    .combine(genomes)
    .filter { r, g -> r.name != g.name }


process align {
    label "mummer"
    label "medium_task"
    publishDir "${params.outdir}/alignments/${ref.baseName}"
    tag { "${ref.baseName} - ${query.baseName}" }

    input:
    set file(ref), file(query) from pairs

    output:
    set val(ref.baseName), file("${query.baseName}.delta") into deltas

    """
    nucmer --threads=${task.cpus} --prefix=${query.baseName} ${ref} ${query}
    """
}


process dnadiff {
    label "mummer"
    label "small_task"
    publishDir "${params.outdir}/diffs/${ref}"

    tag { "${ref} - ${delta.baseName}" }

    input:
    set val(ref), file(delta) from deltas

    output:
    set val(ref), file("${delta.baseName}.snps") into snps
    set val(ref), file("${delta.baseName}.1coords") into onecoords
    set val(ref), file("${delta.baseName}.1delta") into onedeltas
    set val(ref), file("${delta.baseName}.mcoords") into mcoords
    set val(ref), file("${delta.baseName}.mdelta") into mdeltas
    set val(ref), file("${delta.baseName}.qdiff") into qdiff
    set val(ref), file("${delta.baseName}.rdiff") into rdiff
    set val(ref), file("${delta.baseName}.report") into dnadiffReports
    set val(ref), file("${delta.baseName}.unqry") optional true into unqrys

    """
    dnadiff --prefix=${delta.baseName} --delta ${delta}
    """
}


process snp2vcf {
    label "python3"
    label "small_task"
    publishDir { "${params.outdir}/vcfs/${ref.baseName}" }

    tag { "${ref.baseName} - ${snp.baseName}" }

    input:
    set file(ref), file(snp) from references4Snp2Vcf
        .map { [it.baseName, it] }
        .combine( snps, by: 0 )
        .map { rn, r, s -> [r, s] }

    output:
    set val(ref.baseName), file("${snp.baseName}.vcf") into vcfs

    """
    snps2vcf.py \
      --snp ${snp} \
      --reference ${ref} \
      --name ${snp.baseName} \
      --output ${snp.baseName}.vcf
    """
}


process coordsToBED {
    label "posix"
    label "small_task"
    publishDir "${params.outdir}/diffs/${ref}"

    tag { "${ref} - ${coord.baseName}" }

    input:
    set val(ref), file(coord) from mcoords

    output:
    set val(ref), file("${coord.baseName}.bed") into coordsBEDs

    """
    tail -n +5 ${coord} \
    | awk '{ printf "%s\\t%s\\t%s\\t\\.\\t+\\t%s\\tfiltered\\n", \$12, \$1-1, \$2, \$7 }' \
    | sort -k1,1 -k2,2n \
    > ${coord.baseName}.bed
    """
}


process genomeIndex {
    label "samtools"
    label "small_task"

    tag "${ref}"

    input:
    file ref from references4GenomeIndex

    output:
    set file(ref), file("${ref}.fai") into referenceIndex

    """
    samtools faidx ${ref} > ${ref}.fai
    """
}


process bedCoverage {
    label "bedtools"
    label "small_task"
    publishDir "${params.outdir}/bedgraphs/${ref.baseName}"

    tag "${ref.baseName} - ${bed.baseName}"

    input:
    set file(ref), file(index), file(bed) from referenceIndex
        .map { r, i -> [r.baseName, r, i] }
        .combine( coordsBEDs, by: 0 )
        .map { rn, r, i, b -> [r, i, b] }

    output:
    set val(ref.baseName), file("${bed.baseName}.bedgraph") into coverageBEDs

    """
    bedtools genomecov -bga -g ${index} -i ${bed} > ${bed.baseName}.bedgraph
    """
}


process combinedBEDCoverage {
    label "bedtools"
    label "small_task"
    publishDir "${params.outdir}/bedgraphs"

    tag "${ref}"

    input:
    set val(ref), file("*") from coverageBEDs.groupTuple(by: 0)

    output:
    file "${ref.baseName}.bedgraph" into combinedCoverage

    """
    NAMES=( *.bedgraph )

    bedtools unionbedg \
      -header \
      -names \${NAMES[@]%.bedgraph} \
      -i *.bedgraph \
    > ${ref.baseName}.bedgraph
    """
}


/*
percentages = [1, 2, 3, 4, 5, 7, 10, 20, 30, 40, 50,
               60, 70, 80, 90, 93, 95, 96, 97, 98, 99]

process compareCoverage {
    label "python3"
    label "small_task"

    publishDir "${params.outdir}/core_accessory/${bg.baseName}/${perc}"

    tag "${bg.baseName} - ${perc}"

    input:
    file bg from combinedCoverage
    each val perc from percentages

    output:
    file "core.bedgraph" into coreRegions
    file "duplicate.bedgraph" into duplicateRegions
    file "single.bedgraph" into singleRegions

    """
    /usr/bin/env python3

    import os
    import pandas as pd

    table = pd.read_csv("${bg}", sep="\t")
    table.columns = [os.path.split(os.path.splitext(c)[0])[1] for c in table.columns]

    pc = ${perc} / 100

    table2 = table[(table.iloc[:, 3:] > 0).mean(axis=1) > pc]
    table2.to_csv("core.bedgraph", index=False, sep="\t")

    table3 = table2[(table2.iloc[:, 3:] > 1).mean(axis=1) > pc]
    table3.to_csv("duplicate.bedgraph", index=False, sep="\t")

    table4 = table2[(table2.iloc[:, 3:] == 1).mean(axis=1) > pc]
    table4.to_csv("single.bedgraph", index=False, sep="\t")
    """
}
*/
