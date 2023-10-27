# Minigraph-Cactus Pangenome Construction and Downstream Applications

## Table of Contents

* [Abstract](#abstract)
* [Key Reference Material](#key-reference-material)
* [Part 1: Pangenome Graph Construction](#part-1-pangenome-graph-construction)
     * [Cactus Setup](#cactus-setup)
     * [Input Data](#input-data)
     * [Build and Index the Pangenome Graph](#build-and-index-the-pangenome-graph)
* [Part 2: Pangenome Graph Properties](#part-2-pangenome-graph-properties)

## Abstract

This is a tutorial written to support the **Reference Graph Pangenome Data Analysis Hackathon 2023 Nov. 13-17 in Cape Town, South Africa**. The aim is to provide detailed instructions on how to create a pangenome reference graph with Minigraph-Cactus then use it for some downstream analysis like variant calling and genotyping.

Unlike some previous workshops, and most of the existing Cactus documentation, this tutorial will focus on whole-genome human data. As such, it will need to be run over a period of time longer than a typical workshop session. The running times and memory usage of each command will be given. 

Slack (`#refgraph_hackathon_2023`) will probably be the best place to reach out to me (Glenn Hickey) for support. 

## Key Reference Material

Please visit these links for related material and background information before proceeding further. **The first link is essential and should absolutely be consulted before continuing and the rest are highly recommended.** 

* [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md): This is essential to read, and includes several small examples (with data) that should be run before tackling whole-genomes.
* [Minigraph-Cactus Paper](https://doi.org/10.1038/s41587-023-01793-w): The methods are described in detail here.
* [HPRC v1.1 Minigraph-Cactus Instructions](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/mc-pangenomes/hprc-v1.1-mc.md): Commands and explanations in order to exactly reproduce the latest released HPRC graphs. The commands themselves assume a SLURM cluster but can be trivially modified to run on a single computer (remove `--batchSystem slurm`). 
* [HPRC Paper](https://doi.org/10.1038/s41586-023-05896-x): Detailed analysis of the HPRC graph, and examples of many downstream applications of Minigraph-Cactus pangenomes. 
* @jeizenga's [2023 Memphis Workshop](https://github.com/pangenome/MemPanG23/blob/main/lessons/Day_3a_vg_mapping_and_calling.md), which served as an inspiration for this tutorial.

## Part 1: Pangenome Graph Construction

### Cactus Setup

**Important:** We will be using [Cactus v2.6.9](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.6.9) for this tutorial. Be warned that it may not work for newer or older versions.

For simplicity, all cactus will be run in "single-machine" mode via its [docker](https://www.docker.com/) image.  Cactus also supports distributed computing environments via [slurm](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster) and [AWS/Mesos](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/running-in-aws.md).

In order to make sure docker is working, try running the following and verify that you do not get an error. If this step does not work, you will need to consult your local sysadmin. 
```
docker run hello-world
```

You can then pull the Cactus image onto your coputer
```
docker pull quay.io/comparative-genomics-toolkit/cactus:v2.6.9
```

### Input Data

As you've seen in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) (please go back and read it if you haven't already), the input is a list of sample name and genome assembly pairs.  For diploid assemblies, the convention of `SAMPLE.1 / SAMPLE.2` must be used (and dots avoided in sample names otherwise).

In addition to your samples of interest, you should include at least one reference genome. This will allow you to use reference coordinates to, for example, project variants on.  In this example, which is based on a small subset of 4 samples of the [HPRC]((https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/mc-pangenomes/hprc-v1.1-mc.md)) data, we will use GRCh38 and CHM13.

Please copy-paste the following data into `hprc10.seqfile`

```
GRCh38	https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
CHM13	https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
HG00438.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.paternal.f1_assembly_v2_genbank.fa.gz
HG00438.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.maternal.f1_assembly_v2_genbank.fa.gz
HG00621.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.paternal.f1_assembly_v2_genbank.fa.gz
HG00621.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.maternal.f1_assembly_v2_genbank.fa.gz
HG00673.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.paternal.f1_assembly_v2_genbank.fa.gz
HG00673.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.maternal.f1_assembly_v2_genbank.fa.gz
HG00733.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.paternal.f1_assembly_v2_genbank.fa.gz
HG00733.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.maternal.f1_assembly_v2_genbank.fa.gz

```

**If you are making a pangenome graph with your own data, this input listing should be the only part you need to change, but do see the explanation of the options below as some may require adjustments for different data sizes**. Also, nothing changes if you want to use haploid assemblies -- just do not use the `.1` and `.2` suffixes (see `CHM13` and `GRCh38` above).

### Build and Index the Pangenome Graph

I am going to run on 32-cores in order to simulate my understanding of an "average" node on your cluster. As you've seen in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) (please go back and read it if you haven't already), the simplest way to build the graph is with the `cactus-pangenome` command.

Here it is, with an explanation of each option following below.

```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
cactus-pangenome /data/js /data/hprc10.seqfile --outDir /data/hprc10 --outName hprc10 --reference GRCh38 CHM13 \
--filter 2 --haplo --giraffe clip filter --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full \
--gfa clip full --vcf --vcfReference GRCh38 CHM13 --logFile /data/hprc10.log \
--consCores 8
```

For `docker run`:
* `-it`: interactive / tty.  Boilerplate for most command line tools
* `--rm`: save space by removing the container after running
* `-v $(pwd):/data`: mount the current directory to `/data` in the container
* `--user $UID:$GID`: run as the current user in the container
* `quay.io/comparative-genomics-toolkit/cactus:v2.6.9`: the cactus docker image

For `cactus-pangenome`:
* `/data/js`: Scratch directory that will be created for Toil's jobstore
* `/data/hprc10.seqfile`: The input samples and assemblies. This file was created above.
* `--outDir /data/hprc10`: The output directory. All results will be here. Remember anything relative to `/data` in the docker container will end up in your current working directory that you're running `docker run` from.
* `--outName hprc10`: This will be the prefix of all output files.
* `--reference GRCh38 CHM13`: Specify these two samples as reference genomes. Reference samples are indexed a little different in `vg` to make their coordinates easier to use. Also, the first reference given (GRCh38 in this case), is used to anchor the entire graph and is treated differently than the other samples.  Please see [here for more details](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#reference-sample).
* `--filter 2`: Create an Allele-Frequency filtered (AF) graph that contains only nodes and edges supported by at least 2 haplotypes. This can lead to better mapping performance. `--filter 9` was used for the 90-assembly HPRC graph.  
* `--haplo`: We are actually phasing out the Allele-Frequency filtering as described above in favour or dynamic creation of personal pangenomes. Using this option will create the necessary indexes for this functionality.
* `--giraffe clip filter`: Make giraffe indexes for the Allele-Frequency filtered graph in addition to the (default) clipped graph.
* `--viz`: Make an ODGI 1D visualization image for each chromosome.
* `--odgi`: Make an ODGI formatted whole-genome graph
* `--chrom-vg clip filter`: Make VG formatted chromosome graphs for the both the AF filtered and (default) clipped pangenome.
* `--chrom-og`: Make ODGI formatted chromosome graphs for the full (unclipped) graph. Useful for visualization. 
* `--gbz clip filter full`: Make GBZ formatted whole-genome graphs for AF filtered, (default) clipped and full (containing unaligned centromeres) graphs.
* `--gfa clip full`: Make GFA formatted whole-genome graphs for (default) clipped and full graphs.
* `--vcf`: Make a VCF (based on the first reference) version of the graph
* `--vcfReference GRCh38 CHM13`: Specify that we want two VCFs, one for each reference
* `--logFile /data/hprc10.log`: All logging information will end up here in addition to `stderr`.  Important to save!
* `--consCores 8`: Specify 8 threads for each core cactus job (`cactus_consolidated`). By default it will use all cores available on your system.  By reducing to `8`, we attempt to run up to 4 chromosomes at once to save time (assuming 32 cores total). Note that this will increase peak memory usage. 

All of the above is explained in more detail in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md). We are erring on the side of producing lots of different indexes, but it's usually easier than going back and regenerating any forgotten ones.

Here are some details about the resources used. I'm running on 32 cores of a shared server with a slow network drive, so your times could be a bit faster. They are take from `hpr10.log` which lists the wall time and memory usage of each command in the pipeline. The log for the full 90-way HPRC graph can be found [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.log).

1) Minigraph Construction : 106Gi, ~3 hours
2) Minigraph Mapping : ~200Gi (max per-job 40Gi) ~2 hour
3) Cactus Alignment : ~65Gi (max per-job 16Gi) ~1 hour
4) Normalization and Indexing : ~64 Gi ~3 hours 
5) Overall : 11 hours

Here are the output files:
```
4.0K	chrom-alignments
4.0K	chrom-subproblems
278M	hprc10.CHM13.raw.vcf.gz
1.6M	hprc10.CHM13.raw.vcf.gz.tbi
224M	hprc10.CHM13.vcf.gz
1.6M	hprc10.CHM13.vcf.gz.tbi
4.0K	hprc10.chroms
773M	hprc10.d2.dist
1.8G	hprc10.d2.gbz
35G	hprc10.d2.min
65M	hprc10.d2.snarls
1.1G	hprc10.dist
2.2G	hprc10.full.gbz
1.5G	hprc10.full.gfa.gz
9.7G	hprc10.full.hal
9.1G	hprc10.full.og
57M	hprc10.full.snarls
104M	hprc10.gaf.gz
1.8G	hprc10.gbz
769M	hprc10.gfa.fa.gz
1.4G	hprc10.gfa.gz
1.1G	hprc10.hapl
35G	hprc10.min
450M	hprc10.paf
7.3M	hprc10.paf.filter.log
117M	hprc10.paf.unfiltered.gz
279M	hprc10.raw.vcf.gz
1.6M	hprc10.raw.vcf.gz.tbi
907M	hprc10.ri
1.7K	hprc10.seqfile
52M	hprc10.snarls
406K	hprc10.stats.tgz
729M	hprc10.sv.gfa.gz
226M	hprc10.vcf.gz
1.6M	hprc10.vcf.gz.tbi
4.0K	hprc10.viz
```

There are four versions of the graph produced (please see [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#clipping-filtering-and-indexing) for more details), denoted by these preffixes:

* `hprc.sv` : This is the output of `minigraph` and contains only structural variants. The input haplotypes are not embedded as paths
* `hprc.full` : This is a base-level graph containing all sequence that could be assigned to a reference chromosome. Centromeres are included but are unaligned.
* `hprc.` : This is a subgraph of `hprc.full` but with centromeres removed.  This is usually the most relevant graph for analysis.
* `hprc.d2` : This is a subgraph of `hprc.` but with nodes and edges supported by fewer than 2 haplotypes removed. This graph yields better results for read mapping with the original `giraffe` pipeline. We used allele frequency filtering in the HPRC paper but have recently changed `giraffe` so that it is no longer necessary (more details later in the mapping section).

The graphs themselves are present in `.gfa.gz` (standard, text-based), `.gbz` (highly compressed, `vg`) and `.og` (`odgi`) formats. `vg giraffe` mapping requires the `.gbz` and `.hapl` index (or `.gbz`, `.dist` and `.min` for the original pipeline).

There are four VCF files, two each for GRCh38 and CHM13:

* `hprc10.raw.vcf.gz` and `hprc10.CHM13.raw.vcf.gz` : Output of `vg deconstruct` for the GRCh38- and CHM13-based graphs, respectively. These VCFs contain nested variants, and need special attention when using.
* `hprc10.vcf.gz` and `hprc10.CHM13.vcf.gz` : "Flattened" versions of the above VCFs (using `vcfbub`) that do not contain nested variants and will be more useful for standard tools.

For the HPRC, we took some extra normalization steps using `vcfwave` to realign the variants.  This gives a slightly cleaner VCF in some regions.  See [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/hprc-v1.1/doc/mc-pangenomes/hprc-v1.1-mc.md#vcf-postprocessing) for details.

There are four directories:

* `chrom-subproblems` : The per-chromosome *inputs* to `cactus`.  This directory has some useful statistics about the chromosome decomposition such as `contig_sizes.tsv` which shows the amount of sequence from each sample for each reference chromosome as well as `minigraph.split.log` which lists which chromosome each contig gets assigned to and why, along with all contigs excluded from further analysis (counted `_AMBIGUOUS_`) because they didn't align anywhere.  
* `chrom-alignments` : The raw, per-chromosome output of `cactus` including `.vg` and `.hal` files.  Only useful for debugging and / or re-running the last part of the pipeline (indexing and normalization).
* `hprc10.chroms` : Chromosome graphs in `.vg` and `.og` format.  Useful for debugging and visualization. If you are using `GRCh38` as a reference, the unplaced contigs will all get lumped into the `chrOther` graph.  
* `hprc10.viz` : ODGI 1-D visualizations for each chromosome.

## Part 2: Pangenome Graph Properties

### Basic Statistics

The very first thing to check is the size of your graph.  You can do this with `vg stats -lz`:

```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
vg stats -lz /data/hprc10/hprc10.gbz

```

It will show the number of nodes and edges in the graph along with the total sequence length over all nodes:
```
nodes	28195250
edges	38322112
length	3145521882
```

You can compare that to the length of GRCh38 in the graph
```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
vg paths -x /data/hprc10/hprc10.gbz -S GRCh38 -E | awk '{sum += $2} END {print sum}'
```
which is `3099011216`. There is `46510666`bp of additional sequence added to the pangenome from CHM13 and the four samples. Something on the order of a few megabases per sample is reasonable.  If your results are much different, then that is a definite warning sign that something went very wrong.

Looking at the `.full` graph shows how much additional sequence is added by the centromeres.

```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
vg stats -lz /data/hprc10/hprc10.full.gbz
```
An extra gigabase in this case. We cannot effectively index or map to such graphs (centromere alignment is something we are actively working on, though!)
```
nodes	29372041
edges	39555335
length	4213877926
```

You can use `vg paths` to inspect the amount of sequence (total length of all embedded paths) of any given sample of haplotype.  For example

```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
vg paths -x /data/hprc10/hprc10.gbz -S HG00438 -E | awk '{sum += $2} END {print sum}'

```
```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
vg paths -x /data/hprc10/hprc10.gbz -Q HG00438#1 -E | awk '{sum += $2} END {print sum}'

```
```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
vg paths -x /data/hprc10/hprc10.gbz -Q HG00438#2 -E | awk '{sum += $2} END {print sum}'
```

Show that there is `5679580423`bp for `HG00438`, with `2841204110` and `` in its first (paternal) and second (maternal) haplotype, respectively. 

The aforementioned `hprc10/chrom-subproblems/contig_sizes.tsv` gives a breakdown of the length of each haplotype in each chromosome. Can be useful to load into a spreadsheet and/or graph in order to check that all input haplotypes are properly represented in the graph.

`minigraph-cactus` graphs are linearized along the reference genome (GRCh38 in this case).  There is exactly one graph component for each contig in GRCh38.  And each component has exactly two tips or stubs (nodes with zero edges at one of their ends). You can count the tips in the graph with

```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
vg stats -HT /data/hprc10/hprc10.gbz | sed -e 's/heads//g' -e 's/tails//g' | wc -w
```
Giving a result of `354`. This is two times the number of contigs in GRCh38, ` `, which can be inspected with

```
wget -q https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz | grep '>' | wc -l
```
