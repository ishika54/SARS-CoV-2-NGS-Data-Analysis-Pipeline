# SARS-CoV-2 NGS Data Analysis Pipeline

> End-to-end next-generation sequencing (NGS) workflow for viral genome
> characterization, variant calling, and phylogenetic analysis using
> publicly available SARS-CoV-2 data.

---

## Project Overview

This project demonstrates a complete bioinformatics pipeline for analyzing
SARS-CoV-2 sequencing data obtained from the NCBI Sequence Read Archive (SRA).
Starting from raw paired-end FASTQ files, the workflow covers quality assessment,
genome alignment, variant detection, and evolutionary analysis — all executed on
a standard laptop using WSL (Windows Subsystem for Linux).

**Dataset:** `SRR36276613` — PCR tiled amplicon sequencing of SARS-CoV-2
(paired-end, ~92k spots, 26M bases)
**Reference Genome:** Wuhan-Hu-1 (`wuhCor1.fa`, RefSeq: NC_045512.2)

---

## Aim

To characterize genetic variants and evolutionary relationships in SARS-CoV-2
sequencing data by implementing a reproducible NGS analysis pipeline — from raw
read quality control through to phylogenetic tree inference.

---

## Objectives

- Retrieve and prepare NGS data from public repositories (NCBI SRA)
- Assess raw read quality and determine the need for preprocessing
- Align reads to the SARS-CoV-2 reference genome and evaluate mapping efficiency
- Perform post-alignment processing: duplicate marking and read group annotation
- Call genetic variants (SNPs and INDELs) and interpret variant statistics
- Validate variant calls through visual inspection
- Construct and interpret a phylogenetic tree across multiple SARS-CoV-2 variants

---

## Workflow Summary
```
Raw FASTQ (SRA)
      │
      ▼
Quality Control (FastQC + MultiQC)
      │
      ▼
Genome Alignment (BWA-MEM → SAM/BAM)
      │
      ▼
Post-Processing (Sort → MarkDuplicates → AddReadGroups → optional BQSR)
      │
      ▼
Variant Calling (GATK HaplotypeCaller → VCF)
      │
      ▼
Variant Validation & Statistics (IGV + bcftools)
      │
      ▼
Multiple Sequence Alignment (MAFFT)
      │
      ▼
Phylogenetic Tree Construction (FastTree)
```

---

## Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| SRA Toolkit (`prefetch`, `fastq-dump`) | 2.11.3 | Data retrieval from NCBI SRA |
| FastQC | — | Per-read quality assessment |
| MultiQC | — | Aggregated QC reporting |
| BWA | 0.7.17-r1188 | Read alignment (BWA-MEM algorithm) |
| SAMtools | 1.13-4 | BAM sorting, indexing, flagstat |
| GATK | 4.5.0.0 | MarkDuplicates, BQSR, HaplotypeCaller |
| bcftools | — | VCF filtering and querying |
| IGV | — | Alignment and variant visualization |
| MAFFT | 7.490 | Multiple sequence alignment (FFT-NS-2) |
| FastTree | 2.1.11 | Maximum-likelihood phylogenetic tree (GTR model) |

---

## Key Results

### Quality Control
- Both read files passed all critical FastQC metrics: per-base quality,
  adapter content, and N content
- WARN/FAIL flags on duplication and GC content are expected for amplicon-based
  viral sequencing — no trimming was required

### Alignment
- **99.42% of reads** mapped to the SARS-CoV-2 reference genome
- **99.35% properly paired** reads; singletons < 0.1%
- High-depth, genome-wide coverage confirmed by IGV visualisation

### Variant Calling
- Variants called across the ~29.9 kb genome using GATK HaplotypeCaller
- High-confidence SNPs identified: AF ≥ 0.90, depth > 200×, QUAL > 9000
- Example fixed mutation — position 18163, REF=A → ALT=G (AF=1.00, DP=586, GQ=99)
- VCF filtered at QUAL > 50, DP > 20 using bcftools

### Phylogenetic Analysis
- MSA performed on 5 SARS-CoV-2 complete genomes (Alpha, Delta, Omicron,
  reference/early isolates) using MAFFT FFT-NS-2
- ML tree inferred with FastTree (GTR+CAT model)
- Omicron showed the highest evolutionary divergence (longest branch length);
  Delta isolates clustered tightly, reflecting shared recent ancestry

---

## Biological Insights

The analysis confirms the expected evolutionary trajectory of SARS-CoV-2:

- **Delta variants** cluster tightly, indicating shared mutations and a recent
  common ancestor
- **Omicron** is the most genetically divergent taxon, consistent with its
  extensive spike protein mutations
- **Alpha** occupies an intermediate phylogenetic position, reflecting earlier
  divergence from ancestral strains
- **BQSR was intentionally omitted** — no validated known-sites VCF exists for
  SARS-CoV-2, and applying BQSR without known sites risks misclassifying true
  viral mutations as sequencing error; this is standard practice in SARS-CoV-2
  pipelines (ARTIC, iVar, LoFreq)
- The >99% mapping rate and genome-wide coverage depth confirm high sequencing
  data quality and appropriate reference selection

---

## Future Improvements

- Replace HaplotypeCaller with **LoFreq** or **iVar** for haploid-aware variant
  calling
- Expand phylogenetic analysis to include recent lineages (XBB, JN.1, KP.x)
- Automate the full pipeline with **Snakemake** or **Nextflow**
- Integrate **Pangolin** for lineage classification of consensus genomes
- Generate consensus FASTA for downstream strain characterisation

---

## Repository Structure
```
├── raw/          # Raw FASTQ files (SRR36276613_1/2.fastq)
├── ref/          # Reference genome (wuhCor1.fa + BWA/GATK indices)
├── align/        # Sorted, deduplicated, recalibrated BAM files
├── qc/           # FastQC and MultiQC HTML reports
├── variants/     # Raw and filtered VCF files
└── phylo/        # MAFFT MSA (.fasta) and FastTree output (.nwk)
```

---

*Analysis performed on Ubuntu 22.04 (WSL2) · Data: NCBI SRA SRR36276613*
