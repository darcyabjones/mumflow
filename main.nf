#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.genomes = false
params.references = false
params.genes = false
params.promer = false


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

if ( params.genes ) {
    genes = Channel.fromPath(
        params.genes,
        checkIfExists: true,
        type: "file"
    ).map { f -> [f.baseName, f ] }
} else {
    genes = Channel.empty()
}

references.into {
    references4Cross;
    references4Snp2Vcf;
    references4GenomeIndex;
    references4MergeGenes;
}


// Exit if gff provided that can't be matched.
// Filter out refs without gff provided.
// Note we could use join instead of combine, which would allow us
// to get rid of the filter, but we wouldn't be able fail on unmatched gff.
refWithGenes = references4MergeGenes
    .map { f -> [ f.baseName, f ] }
    .combine( genes, by: 0 )
    .map { n, f, g ->
        if ( f == null || f == '' ) {
            log.error "The annotation file ${g.name} specified in --genes could " +
                "not be matched to any reference genome names."

            log.error "Please make sure the annotation file and reference genomes " +
                "have the same basename (up to the last extension)."

            exit 1
        };
        [n, f, g]
    }
    .filter { n, f, g -> (g == null || g == '') }


pairs = references4Cross
    .combine(genomes)
    .filter { r, g -> r.name != g.name }


process mumalign {
    label "mummer"
    label "medium_task"
    publishDir "${params.outdir}/alignments/${ref.baseName}"
    tag { "${ref.baseName} - ${query.baseName}" }

    input:
    set file(ref), file(query) from pairs

    output:
    set val(ref.baseName), file("${query.baseName}.delta") into deltas

    script:
    if (params.promer) {
        exe = "promer"
    } else {
        exe = "nucmer"
    }

    """
    ${exe} \
      --threads=${task.cpus} \
      --maxmatch \
      --prefix=${query.baseName} \
      ${ref} \
      ${query}
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


/*
 * Convert coords to awk.
 * Note that the output of show-coords and dnadiff are slightly different.
 * This is for mcoords from dnadiff.
 */
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
    awk '{ printf "%s\\t%s\\t%s\\t\\.\\t+\\t%s\\tfiltered\\n", \$12, \$1-1, \$2, \$7 }' ${coord} \
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

referenceIndex.into {
    referenceIndex4BedCoverage;
    referenceIndex4MakeWindows;
    referenceIndex4MakePBWindows;
    referenceIndex4PlotCoverages;
}

process bedCoverage {
    label "bedtools"
    label "small_task"

    tag "${ref.baseName} - ${bed.baseName}"

    input:
    set file(ref), file(index), file(bed) from referenceIndex4BedCoverage
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
    publishDir "${params.outdir}/bedgraphs/${ref}"

    tag "${ref}"

    input:
    set val(ref), file("*") from coverageBEDs.groupTuple(by: 0)

    output:
    set val(ref), file("${ref}.bedgraph") into combinedCoverage

    """
    NAMES=( *.bedgraph )

    bedtools unionbedg \
      -header \
      -names \${NAMES[@]%.bedgraph} \
      -i *.bedgraph \
    > ${ref}.bedgraph
    """
}

combinedCoverage.into {
    combinedCoverage4GetPBCoverage;
    combinedCoverage4FindPAV;
}


window_sizes = [10000, 50000, 100000]

process makeWindows {

    label "bedtools"
    label "small_task"

    tag "${window_size}"

    input:
    set file(ref), file(index) from referenceIndex4MakeWindows
    each window_size from window_sizes

    output:
    set val(ref.baseName), val(window_size),
        file("windows_${window_size}.bed") into windows

    """
    bedtools makewindows -g "${index}" -w "${window_size}" > windows.tmp.bed

    # Doing this in 2 steps (rather than piping) is necessary to avoid
    # overlayfs requirements in singularity.
    sort -k1,1 -k2,2n windows.tmp.bed > "windows_${window_size}.bed"

    rm -f windows.tmp.bed
    """
}

process makePBWindows {
    label "bedtools"
    label "small_task"

    tag "${ref.baseName}"

    input:
    set file(ref), file(index) from referenceIndex4MakePBWindows

    output:
    set val(ref.baseName), file("windows.bed") into pbWindows

    """
    bedtools makewindows -g "${index}" -w 1 > windows.tmp.bed

    # Doing this in 2 steps is necessary to avoid overlayfs
    # requirements in singularity. (because sort uses tmp files).
    sort -k1,1 -k2,2n windows.tmp.bed > windows.bed

    rm -f windows.tmp.bed
    """
}

process getPBCoverage {

    label "bedtools"
    label "small_task"

    tag "${ref}"

    input:
    set val(ref), file(pb_window), file(coverage) from pbWindows
        .combine( combinedCoverage4GetPBCoverage, by: 0 )

    output:
    set val(ref), file("pbcov.bedgraph") into pbCoverages

    """
    bedtools intersect \
      -a "${coverage}" \
      -b "${pb_window}" \
      -sorted \
      -header \
    > "pbcov.bedgraph"
    """
}

process getMeanWindowedCoverage {

    label "bedtools"
    label "small_task"

    tag "${ref} - ${window_size}"

    publishDir "${params.outdir}/bedgraphs/${ref}"

    input:
    set val(ref), val(window_size), file("windows.bed"),
        file("pbcov.bedgraph") from windows
            .combine( pbCoverages, by: 0 )

    output:
    set val(ref), val(window_size),
        file("cov_mean_${window_size}.bedgraph") into meanWindowedCoverage

    """
    NCOLS=\$(head -n 1 pbcov.bedgraph | wc -w)
    COLS=\$(seq 4 \${NCOLS} | tr '\n' ',' | sed 's/,\$//')

    head -n 1 "pbcov.bedgraph" > "cov_mean_${window_size}.bedgraph"

    bedtools map \
      -a "windows.bed" \
      -b "pbcov.bedgraph" \
      -c "\${COLS}" \
      -o mean \
    >> "cov_mean_${window_size}.bedgraph"
    """
}


process plotCoverages {

    label "r"
    label "small_task"

    tag "${ref} - ${window_size}"

    publishDir "${params.outdir}/coverage_plots/${ref}"

    input:
    set val(ref), file(index), val(window_size), file(bg) from referenceIndex4PlotCoverages
        .map { f, i -> [f.baseName, i] }
        .combine( meanWindowedCoverage, by: 0 )

    output:
    set val(ref), val(window_size), file("${window_size}") into coveragePlots

    """
    plot_circos.R --bedgraph "${bg}" --faidx "${index}" --outdir "${window_size}"
    """
}


process findPAV {

    label "python3"
    label "small_task"
    tag "${ref}"

    publishDir "${params.outdir}/pav/${ref}"

    input:
    set val(ref), file(bg) from combinedCoverage4FindPAV

    output:
    set val(ref), file("pavs.bedgraph") into foundPAVs

    """
    find_pavs.py \
      --infile "${bg}" \
      --outfile pavs.bedgraph \
      --tol 20 \
      --min-length 50 \
      --proportion-repeats 1.0
    """
}


process pavGenes {

    label "bedtools"
    label "small_task"
    tag "${ref}"

    publishDir "${params.outdir}/pav/${ref}"

    input:
    set val(ref), file("pavs.bedgraph"), file("genome.fasta"),
        file(genes) from foundPAVs
            .join( refWithGenes, remainder: false, by: 0 )

    output:
    set val(ref), file("gene_pavs.bedgraph") into genePAVs

    """
    bedtools intersect \
      -a pavs.bedgraph \
      -b "${genes}" \
      -header \
      -u \
    > gene_pavs.bedgraph
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
