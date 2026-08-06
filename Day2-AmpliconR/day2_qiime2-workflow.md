
------------------------------------------------------------------------

# Day 2: Amplicon analysis using QIIME 2

**Author**: Veronica Roman-Reyna

**Topic:** 16S rRNA amplicon sequencing analysis using QIIME 2.

**Data type:** Demultiplexed paired-end Illumina reads (2×250, 16S rRNA V4 region).

**Source data:** Subset of the Mothur MiSeq SOP tutorial dataset.

**Qiime2 version:** qiime2-2026.7

**Reference workflow:** Adapted from the [QIIME 2 Moving Pictures tutorial](https://amplicon-docs.qiime2.org/en/stable/tutorials/moving-pictures/)

------------------------------------------------------------------------

## What You'll Need

- This document for instructions and commands. Preferably as a separate tab.

- The terminal.

- A second browser tab open with [QIIME 2 View](https://view.qiime2.org/) — this is where you'll upload `.qzv` visualization files to inspect your results.

## Assumptions

This workflow assumes that your sequencing data:

- Is for 16S rRNA, not for ITS (fungal).

- Samples have been demultiplexed (split into individual per-sample fastq files).

- Non-biological nucleotides have been removed (e.g., primers, adapters, linkers, etc.)

- Paired-end sequencing data: the forward and reverse fastq files contain reads in matched #order.

If these criteria are not true for your data, Qiime2 offers other commands and workflows to demultiplex, clean reads, and use other non-bacterial databases.

## Qiime2 file types

QIIME2 stores data and results in special zipped archive formats: `.qza` and `.qzv`.

`.qza` ("**Q**IIME **Z**ipped **A**rtifact") holds data/results for further processing.

`.qzv` ("**Q**IIME **Z**ipped **V**isualization") holds interactive visualizations you view in a web browser. You'll upload `.qzv` files to [view.qiime2.org](https://view.qiime2.org/) to see plots, tables, and interactive summaries.

------------------------------------------------------------------------

## Starting an Interactive Session

Request compute resources from the cluster's job scheduler (SLURM) so you have a dedicated space to run QIIME 2:

``` bash
salloc -N 1 -n 2 -t 4:00:00 --account=one_sc_default
```

**What this does:**

| Flag            | Meaning                            |
|-----------------|------------------------------------|
| `-N 1`          | Request 1 compute node             |
| `-n 2`         | Request 2 CPU tasks/cores         |
| `-t 4:00:00`    | Reserve the session for 4 hours    |
| `--account=one_sc_default` | Bill the job to the `one_sc_default` account |

> ⏲️ **About the `-t` (time) flag:** 3 hours is just a suggestion, since it comfortably covers the first part of today's session before lunch break. Feel free to request more or less time. The important part is what happens at the end of that window: **once the time limit is reached, your session ends automatically and any jobs still running at that moment are interrupted**.
>
> ⚠️ **Important:** Once you're in an interactive session, **do not close the terminal window**. Closing it ends the session and erases everything held in memory, including any variables and activated environments you've set up. If it closes, you'll need to start over from this step.

------------------------------------------------------------------------

## Background

**Raw data:** This is a subset of paired-end Illumina data (fastq files) from the Mothur MiSeq SOP tutorial, using 2×250 base pair sequencing of the 16S rRNA V4 region.

**Commands:** The steps below are a modified version of the official QIIME 2 tutorial commands.

> **Note:** This script processes **paired-end** reads (`sample_R1.fastq` and `sample_R2.fastq`), where each sample has two read files (forward and reverse). The same general workflow also works for **single-end** reads (`sample.fastq`), with minor changes to the import step.

------------------------------------------------------------------------

## Workflow Steps:

0.  Define Variables
1.  Load Qiime2
2.  Import data
3.  Summarize the Imported Data
4.  Denoise Reads with DADA2 plugin
5.  Visualize the denoised data
6.  Build a Phylogenetic Tree
7.  Alpha and Beta Diversity
8.  Taxonomic Classification
9.  ANCOM-BC2

## Step 0: Define Variables

Setting a variable makes long file paths easier to reuse throughout the script instead of retyping them every time.

``` bash
DAY2=/storage/group/one/default/workshop/2026/day2_16S/qiime2
```

- `$DAY2` points to the shared folder containing the workshop's input files (metadata, manifest, reference classifiers).
- `$WORK` is a variable that indicates the path to your **personal work directory** on the cluster. This a default variable defined in the supercomputer. This is where all the files you generate will be saved.

> **Remember:** User-defined variables only exist for the current session. Every time you start a new interactive session, you'll need to define `$DAY2` again.

------------------------------------------------------------------------

## Step 1: Load QIIME 2

Load the software environment that contains QIIME 2 and its dependencies:

``` bash
module load anaconda/2023.09
```
``` bash
conda activate /storage/group/one/default/workshop/2026/envs/rachis-qiime2-2026.7
```

- `module load anaconda` makes the Anaconda package manager available.
- `conda activate` switches into the specific pre-built QIIME 2 environment for this workshop.

> **Note:** This may take a minute to load. You'll know it worked when you see `(rachis-qiime2-2026.7)` appear at the beginning of your command line prompt. That confirms you're now working inside the QIIME 2 environment.

------------------------------------------------------------------------

## Step 2: Import Reads and Metadata

QIIME 2 needs your raw sequencing data and sample information converted into its own artifact format before it can do anything with them. The data is already in the supercomputer, and the manifest has the required information.

**Two supporting files are used here:**

- **`data-manifest.txt`** : tells QIIME 2 exactly where each sample's `fastq` files are located on the cluster.

- **`metadata.tsv`** : contains information about your experiment (e.g., sample names, treatment groups, collection site).

The first command *imports* the demultiplexed paired-end reads into a `.qza` artifact (`pe-demux.qza`).

The second command creates a visualization of your metadata table so you can double-check it looks correct.

``` bash
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $DAY2/data-manifest.txt \
    --output-path $WORK/pe-demux.qza \
    --input-format PairedEndFastqManifestPhred33V2
```
``` bash
qiime metadata tabulate \
  --m-input-file $DAY2/metadata.tsv \
  --o-visualization $WORK/metadata-viz.qzv
```

------------------------------------------------------------------------

## Step 3: Summarize the Imported Data

This generates `pe-demux.qzv`, an interactive summary showing read counts per sample and sequence quality plots.

``` bash
qiime demux summarize \
    --i-data $WORK/pe-demux.qza \
    --o-visualization $WORK/pe-demux.qzv
```

👀 **Download the file** `pe-demux.qzv` **and upload it to [view.qiime2.org](https://view.qiime2.org/)**.

🔍 What to look for in the forward and reverse reads:

- How is the read count distribution?

- What is the quality score across a sequence?

🔍 The quality plots help to:

- Get an overview of the read distribution and the quality of data.

- Decide if and what type of cleaning is required (informed step 4).

------------------------------------------------------------------------

## Step 4: Denoise Reads with DADA2

DADA2 filters out noisy/low-quality reads, corrects sequencing errors, removes chimeras, and merges paired-end reads into **ASVs (Amplicon Sequence Variants)**.

``` bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $WORK/pe-demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 240 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 160 \
  --o-representative-sequences $WORK/asv-seqs.qza \
  --o-table $WORK/asv-table.qza \
  --o-denoising-stats $WORK/denoising-stats.qza \
  --o-base-transition-stats $WORK/base-transition-stats.qza
```

**Key parameters:**

- `--p-trim-left-f` / `--p-trim-left-r`: bases to trim from the *start* of forward/reverse reads (e.g., primers or low-quality leading bases).

- `--p-trunc-len-f` / `--p-trunc-len-r`: length to *truncate* forward/reverse reads to, cutting off low-quality trailing bases (chosen based on the quality plots from Step 3).

**Outputs:**

- `asv-table.qza`: a feature table of how many reads of each ASV appear in each sample.

- `asv-seqs.qza`: the representative DNA sequence for each ASV.

- `denoising-stats.qza`: a log of how many reads survived filtering, merging, and chimera removal at each stage (this is the one we'll focus on today).

- `base-transition-stats.qza`: diagnostic output showing the error model DADA2 used to predict nucleotide substitutions (not covered in this workshop, but useful for troubleshooting).

> ⏱️ This step may take **\~6 minutes**. The terminal will appear to hang; that's normal, it's still working.

Tabulate the denoising stats so you can view them:

``` bash
qiime metadata tabulate \
    --m-input-file $WORK/denoising-stats.qza \
    --o-visualization $WORK/denoising-stats.qzv
```

👀 **Download the file** `denoising-stats.qzv` **and upload it to [view.qiime2.org](https://view.qiime2.org/)**.

Understanding `denoising-stats.qzv`:

- *Denoised* : reads (forward and reverse, separately) after DADA2 has corrected sequencing errors. DADA2 uses a statistical error model to figure out which small differences between reads are true biological variation versus sequencing mistakes, and "corrects" the likely mistakes.

- *Merged* : reads where the corrected forward and reverse read for a given fragment have been successfully overlapped and joined into one full-length sequence. Pairs that don't overlap well enough (e.g., too much low-quality trimming) fail to merge and are dropped here. This is often where you lose the most reads, so it's a good column to check if your final counts look low.

- *Non-chimeric* : merged reads remaining after chimera removal. A *chimera* is an artificial sequence accidentally created during PCR, when the amplification stops partway through one template and finishes on a different one, "stitching" together two unrelated sequences into one read. Because chimeras look like plausible-but-fake ASVs, DADA2 identifies and removes them so they don't get counted as real organisms.

🔍 What to look for:

- How many reads were retained at each processing step.

- Did you lose a lot of reads or are you retaining more than 80%?

------------------------------------------------------------------------

## Step 5: Visualize the denoised data

Qiime2 then creates a table to summarize the denoised and merged data. `summary.qzv`, indicates how many reads/features each sample has. This is useful for spotting samples that may need to be excluded due to low depth (low read number).

``` bash
qiime feature-table summarize \
    --i-table $WORK/asv-table.qza \
    --o-feature-frequencies $WORK/feature-frequencies.qza \
    --o-sample-frequencies $WORK/sample-frequencies.qza \
    --o-summary $WORK/summary.qzv \
    --m-metadata-file $DAY2/metadata.tsv
```

👀 **Download the file** `summary.qzv` **and upload it to [view.qiime2.org](https://view.qiime2.org/)**.

🔍 What to look for:

- How many unique features are in the sample set.

- How many samples have low read frequency

- What is the cut-off for sampling depth.

`asv-seqs.qzv`, is a browsable, BLAST-linkable list of the actual ASV sequences.

``` bash
qiime feature-table tabulate-seqs \
    --i-data $WORK/asv-seqs.qza \
    --o-visualization $WORK/asv-seqs.qzv
```

👀 **Download the file** `asv-seqs.qzv` **and upload it to [view.qiime2.org](https://view.qiime2.org/)**.

🔍 Why is this useful: having the sequences can help to identify potential hits and whether the taxonomy classifiers are proving false-positives or negatives.

------------------------------------------------------------------------

## Step 6: Build a Phylogenetic Tree

A phylogenetic tree gives structure to the data by showing the evolutionary relationships between the ASVs. This is required for phylogeny-based diversity metrics in the next steps.

``` bash
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences $WORK/asv-seqs.qza \
    --o-alignment $WORK/aligned-asv-seqs.qza \
    --o-masked-alignment $WORK/masked-aligned-asv-seqs.qza \
    --o-tree $WORK/unrooted-tree.qza \
    --o-rooted-tree $WORK/rooted-tree.qza
```

This single command performs four actions in sequence: multiple sequence alignment (MAFFT), masking of highly variable alignment positions, building an unrooted tree, and rooting that tree (FastTree). The final `rooted-tree.qza` is what you'll use downstream.

------------------------------------------------------------------------

## Step 7: Alpha and Beta Diversity

*Alpha rarefaction :* checks whether you've sequenced deeply enough to capture the full diversity within each sample:

``` bash
qiime diversity alpha-rarefaction \
    --i-table $WORK/asv-table.qza \
    --i-phylogeny $WORK/rooted-tree.qza \
    --p-min-depth 10 \
    --p-max-depth 2500 \
    --m-metadata-file $DAY2/metadata.tsv \
    --o-visualization $WORK/alpha-rarefaction.qzv
```

**Key parameter:**

- `--p-min-depth` and `--p-max-depth` set the range of sequencing depths tested. The resulting curves should "level off" if sequencing depth was sufficient.

? What would happen if we use --p-max-depth 5000?

👀 **Download the file** `alpha-rarefaction.qzv` **and upload it to [view.qiime2.org](https://view.qiime2.org/)**.

🔍 What to look for: is the curve reaching a plateau?

*Core diversity metrics* : calculates a standard suite of alpha diversity (within-sample) and beta diversity (between-sample) metrics all at once:

``` bash
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny $WORK/rooted-tree.qza \
    --i-table $WORK/asv-table.qza \
    --p-sampling-depth 2500 \
    --m-metadata-file $DAY2/metadata.tsv \
    --output-dir $WORK/core-metrics-results
```

**Key parameter:**

- `--p-sampling-depth 2500` rarefies (subsamples) every sample down to 2,500 reads so diversity comparisons across samples are fair. Results are saved into a new folder, `core-metrics-results`, rather than individual named files.

**Outputs:** Four `*_emperor.qzv` files. Each is an interactive *Emperor* plot. A 3D PCoA (Principal Coordinates Analysis) plot where each dot is a sample, and samples that are closer together are more similar to each other in microbial community composition.

The four plots differ in *how* "similarity" is calculated:

|   | Presence/absence only | Abundance-weighted |
|---------------------|--------------------------|-------------------------|
| **Uses phylogenetic tree** | `unweighted_unifrac_emperor.qzv` | `weighted_unifrac_emperor.qzv` |
| **No phylogenetic tree** | `jaccard_emperor.qzv` | `bray_curtis_emperor.qzv` |

👀 **Download the files** `*_emperor.qzv` **and upload each one of them to [view.qiime2.org](https://view.qiime2.org/)**.

In each Emperor plot, you can color points by any column in your metadata to visually check whether samples cluster by that variable.

------------------------------------------------------------------------

## Step 8: Taxonomic Classification

This step assigns a taxonomic identity (mainly up to genus) to each ASV by comparing it against a reference database. We'll compare two commonly used reference databases: [**SILVA**](https://www.arb-silva.de/) and [**GTDB**](https://gtdb.ecogenomic.org/).

> **Note:** Always double-check classifier versions. Mismatched versions can cause errors. These databases are for Qiime2-2026.07 version.

### SILVA classifier for 16S rRNA V4 region.

``` bash
qiime feature-classifier classify-sklearn \
    --i-classifier $DAY2/silva-138-2-v4-classifier.qza \
    --i-reads $WORK/asv-seqs.qza \
    --o-classification $WORK/silva-taxonomy.qza
```
``` bash
qiime metadata tabulate \
    --m-input-file $WORK/silva-taxonomy.qza \
    --o-visualization $WORK/silva-taxonomy.qzv
```
``` bash
qiime taxa barplot \
    --i-table $WORK/asv-table.qza \
    --i-taxonomy $WORK/silva-taxonomy.qza \
    --m-metadata-file $DAY2/metadata.tsv \
    --o-visualization $WORK/silva-taxonomy-plots.qzv
```

### GTDB classifier for V4 region

``` bash
qiime feature-classifier classify-sklearn \
    --i-classifier $DAY2/gtdb-r232.0-2026.7.0-classifier.qza \
    --i-reads $WORK/asv-seqs.qza \
    --o-classification $WORK/gtdb-taxonomy.qza
```
``` bash
qiime metadata tabulate \
    --m-input-file $WORK/gtdb-taxonomy.qza \
    --o-visualization $WORK/gtdb-taxonomy.qzv
```
``` bash
qiime taxa barplot \
    --i-table $WORK/asv-table.qza \
    --i-taxonomy $WORK/gtdb-taxonomy.qza \
    --m-metadata-file $DAY2/metadata.tsv \
    --o-visualization $WORK/gtdb-taxonomy-plots.qzv
```

**Each classifier block follows the same three-step pattern:**

1\. `classify-sklearn` : assigns taxonomy to each ASV using a pre-trained naive Bayes classifier.

2\. `metadata tabulate` : creates a browsable table of the taxonomy assignments.

3\. `taxa barplot` : creates an interactive stacked bar chart showing the relative taxonomic composition of each sample, viewable at multiple taxonomic levels (phylum, class, genus, etc.).

👀 **Download the files** `*-taxonomy-plots.qzv` **and upload each one to [view.qiime2.org](https://view.qiime2.org/)**.

🔍 What to look for:

- Select a Taxonomic level (Level 6 is for genus).

- Select how to sort the samples (e.g., index, time, dpw)

- You can download the figure file (SVG) or the taxa counts (CSV).

Comparing the SILVA and GTDB results side by side is a good exercise in how reference database choice can affect taxonomic classification outcomes.

## Step 9: ANCOM-BC2: *An*alysis of *Co*mpositions of *M*icrobiomes with *B*ias *C*orrection *2*

Everything so far has described *how similar samples are* (diversity) or *what's present* (taxonomy). ANCOM-BC2 answers a different question: **which specific features (ASVs/taxa) differ significantly in abundance between groups**, for example, between two time points.

> **What ANCOM-BC2 does, conceptually:** Microbiome sequencing data is *compositional*, read counts are relative to the total reads sequenced for that sample, not absolute abundances in the environment. This means standard statistical tests (built for absolute counts) can give misleading results if applied directly. ANCOM-BC2 corrects for this compositional bias and for uneven sampling depth between samples before testing each feature for a significant association with your variable of interest.
>
> ⚠️ **Before you use it:** ANCOM-BC2 is a plugin in QIIME 2, but like any statistical method, it comes with assumptions (e.g., about how much of the community is expected to change) and limitations. Don't treat it as a black box; differential abundance testing is a genuinely hard problem, and a significant q-value is best treated as a hypothesis worth following up with new samples, not as a conclusion on its own. It's worth reading the [composition plugin reference](https://amplicon-docs.qiime2.org/en/stable/references/plugins/composition/) to understand what ANCOM-BC2 can and can't tell you before concluding your results.

``` bash
qiime composition ancombc2 \
  --i-table $WORK/asv-table.qza \
  --m-metadata-file $DAY2/metadata.tsv \
  --p-fixed-effects-formula 'time' \
  --o-ancombc2-output $WORK/ancombc2-output.qza
```
``` bash
qiime composition da-barplot \
  --i-data $WORK/ancombc2-output.qza \
  --i-taxonomy $WORK/gtdb-taxonomy.qza \
  --o-visualization $WORK/ancombc2-barplot-gtdb.qzv
```
``` bash  
qiime composition da-barplot \
  --i-data $WORK/ancombc2-output.qza \
  --i-taxonomy $WORK/silva-taxonomy.qza \
  --o-visualization $WORK/ancombc2-barplot-silva.qzv
```

**What's happening in each command:**

- `qiime composition ancombc2` : runs the differential abundance test itself. `--p-fixed-effects-formula 'time'` tells ANCOM-BC2 which metadata column to test features against. Here it's asking: does this feature's abundance change significantly with `time`?

- `qiime composition da-barplot` (run for each database used): turns the raw ANCOM-BC2 results into a readable bar plot of differentially abundant (**da** = differentially abundant) features, once labeled with GTDB or SILVA taxonomy.

👀 **Download the files** `ancombc2-barplot-*.qzv` **and upload each one to [view.qiime2.org](https://view.qiime2.org/)**.

🔍 What to look for:

- *Define the model variables* : confirm which metadata column(s) were tested (as set by `--p-fixed-effects-formula`), and use the barplot's filtering options to subset which features/taxonomic levels are displayed.

- *Check the q- and p-values* to assess significance. The *p-value* is the significance of a single feature's test before correction; the *q-value* is the p-value adjusted for multiple testing (since you're testing many features at once, some "significant" p-values are expected by chance alone — the q-value accounts for this). Generally, rely on the *q-value* (commonly using a threshold like q \< 0.05) rather than the raw *p-value* when deciding what counts as significant.

- *Identify the significant differentially abundant genera* : which taxa pass your significance threshold, and in which direction (increasing or decreasing) with respect to your variable of interest.f

## 🎉 Congratulations!

You've completed the full QIIME 2 workflow. From raw paired-end FASTQ files to differential abundance testing. Along the way you imported and quality-checked your reads, denoised them into ASVs with DADA2, built a phylogenetic tree, explored alpha and beta diversity, assigned taxonomy with two different reference databases, and used ANCOM-BC2 to pinpoint which taxa actually differ between groups.

That's the same core pipeline used in real published microbiome studies; nice work getting through it today!
