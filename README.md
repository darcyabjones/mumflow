# mumflow
A nextflow pipeline to run the mummer pipelines (align genomes, call snps) and process outputs into sane formats.

Mumflow pairwise aligns a number of assembled genomes to a set of reference genomes, and outputs common
results for presence absence variant (PAV), snp, and repeat finding in standard formats.

## Requirements

The pipeline uses [Nextflow](https://www.nextflow.io/) to execute tasks.
Main steps require samtools, mummer (version >= 4), and bedtools.
Scripts require python (Version >=3.5) and R (version >=3.5).
Required R packages are R `ape`, `tidyr`, and `circlize`.

A [singularity](https://www.sylabs.io/singularity/) container providing all dependencies is available at <https://cloud.sylabs.io/library/darcyabjones/default/mumflow#>.
To install singularity see the [latest docs](https://www.sylabs.io/docs/) or see if your package manager has it available.
To build the container yourself you will need version >=3.2.
Running the build container should work with any version >= 3.0.

To pull the container run...

```bash
#NB its usually a good idea to take note of the actual latest version number.
singularity pull library://darcyabjones/default/mumflow:latest
```

## Usage

```bash
nextflow run mumflow -resume --references "refs/*.fasta" --genomes "genomes/*.fasta"
```

To use a pulled singularity container...

```bash
singularity pull library://darcyabjones/default/mumflow:v0.0.2
nextflow run mumflow -resume -with-singularity ./mumflow_v0.0.2.sif --references "refs/*.fasta" --genomes "genomes/*.fasta"
```
