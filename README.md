# SNP Pair QC Pipeline

![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen)
![Version](https://img.shields.io/badge/version-dev1.0-blue)

A Nextflow pipeline for SNP-based quality control to verify sample concordance between paired tissue and cfDNA samples through variant allele frequency (VAF) correlation analysis.

## Overview

This pipeline validates that tissue and cfDNA samples originate from the same individual by comparing variant allele frequencies (VAF) at known SNP sites. It's designed for whole exome sequencing (WES) data and uses Pearson correlation to assess sample matching, helping identify sample swaps, contamination, or labeling errors.

### Key Features

- **Automated sample matching**: Validates tissue-cfDNA sample pairs using SNP correlation
- **High-throughput processing**: Parallel processing with configurable CPU allocation
- **Quality filtering**: Applies mapping quality (MAPQ) and base quality (BQ) thresholds
- **Comprehensive reporting**: Generates detailed QC reports with correlation statistics
- **nf-core modules**: Leverages community-tested workflow components
- **HPC-ready**: Configured for SLURM execution with Singularity containers

## Quick Start

### Prerequisites

- **Nextflow**: ≥ 21.04.0
- **Singularity**: ≥ 3.7.0 (if using containers)
- **Java**: ≥ 11

### Installation

```bash
# Clone the repository
git clone https://github.com/vanreality/snp_pair_qc.git
cd snp_pair_qc

# Verify Nextflow installation
nextflow -version
```

### Basic Usage

```bash
nextflow run snp_pair_qc.nf \
  --input_samplesheet samplesheet.csv \
  --fasta reference.fa \
  --fasta_index reference.fa.fai \
  --known_sites_tsv known_snps.tsv \
  --outdir results/ \
  -profile singularity
```

## Input Preparation

### Samplesheet Format

Create a CSV samplesheet with the following columns:

```text
sample,fastq1,fastq2,pileup
SAMPLE_001,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,/path/to/cfDNA_pileup.tsv.gz
SAMPLE_002,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,/path/to/cfDNA_pileup.tsv.gz
```

**Column descriptions:**

- `sample`: Unique sample identifier
- `fastq1`: Path to tissue forward reads (R1)
- `fastq2`: Path to tissue reverse reads (R2)
- `pileup`: Path to pre-generated cfDNA pileup file

### Known Sites File

The known SNP sites file should be in TSV format (VCF-like) with these columns:

```text
chr    pos      id    ref    alt    qual    filter    info
chr1   12345    .     A      G      .       .         AF=0.234;...
```

The pipeline extracts allele frequency (AF) from the INFO field.

## Pipeline Workflow

```text
Input: Tissue FASTQ + cfDNA Pileup
    ↓
1. BWA-MEM2 Index (reference genome)
    ↓
2. BWA-MEM2 Alignment (tissue FASTQ → BAM)
    ↓
3. Picard MarkDuplicates (remove PCR duplicates)
    ↓
4. SNP Pileup Generation (tissue BAM → pileup)
    ↓
5. VAF Correlation Analysis (tissue pileup ↔ cfDNA pileup)
    ↓
Output: QC Reports + Pipeline Info
```

### Workflow Steps

1. **Genome Indexing**: Creates BWA-MEM2 index from reference genome
2. **Read Alignment**: Maps tissue FASTQ reads to reference genome
3. **Duplicate Marking**: Identifies and removes PCR/optical duplicates
4. **Pileup Generation**: Extracts allele counts at known SNP sites from tissue BAM
5. **Correlation Analysis**: Compares tissue and cfDNA VAF values and computes QC metrics

## Quality Control Criteria

The pipeline assesses sample matching using these criteria:

- **Minimum SNP count**: ≥300 common SNP sites (default, configurable)
- **Correlation range**: Pearson correlation between 0.6 and 0.8
- **P-value**: Statistical significance of correlation

### QC Status

- ✓ **PASS**: Both criteria met - samples likely from same individual
- ✗ **FAIL**: One or both criteria not met - potential sample mismatch

## Output Structure

```text
results/
├── bwamem2_index/
│   └── reference_index/
├── bwamem2_mem/
│   └── sample.bam
├── snp_pileup/
│   └── sample_pileup.tsv.gz
├── snp_correlation/
│   └── sample_report.txt        # QC report with correlation analysis
└── pipeline_info/
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.html
```

### Output Files

- **`*_pileup.tsv.gz`**: Compressed TSV with allele counts at SNP positions
- **`*_report.txt`**: Detailed QC report including:
  - Input parameters
  - Filtering statistics
  - Pearson correlation coefficient and p-value
  - QC pass/fail status with reasoning

## Configuration

### Pipeline Parameters

Key parameters can be set via command line or config files:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_samplesheet` | Input CSV samplesheet | Required |
| `--fasta` | Reference genome FASTA | Required |
| `--fasta_index` | Reference genome index (.fai) | Required |
| `--known_sites_tsv` | Known SNP sites TSV | Required |
| `--outdir` | Output directory | Required |

### Profile Configuration

The pipeline includes a Singularity profile for HPC environments:

```bash
# Use Singularity containers
nextflow run snp_pair_qc.nf -profile singularity ...

# Custom config
nextflow run snp_pair_qc.nf -c custom.config ...
```

### Resource Configuration

Default resource allocations (configurable in `nextflow.config`):

| Process | CPUs | Memory | Time | Queue |
|---------|------|--------|------|-------|
| BWA-MEM2 Index | 16 | 32 GB | 24 h | cn-long |
| BWA-MEM2 Alignment | 16 | 32 GB | 24 h | cn-long |
| MarkDuplicates | 8 | 16 GB | 24 h | cn-long |
| SNP Pileup | 8 | 16 GB | 24 h | cn-long |
| SNP Correlation | 8 | 16 GB | 24 h | cn-long |

## Python Tools

The pipeline includes two standalone Python scripts that can be used independently:

### snp_pileup.py

Generates pileup data from WES BAM files at known SNP sites.

```bash
python bin/snp_pileup.py \
  --input-bam tissue.bam \
  --known-sites snps.tsv \
  --output sample_prefix \
  --min-mapq 20 \
  --min-bq 13 \
  --ncpus 8
```

#### Options

- `--input-bam`: Input BAM file
- `--known-sites`: Known SNP sites TSV file
- `--output`: Output file prefix
- `--min-mapq`: Minimum mapping quality (default: 20)
- `--min-bq`: Minimum base quality (default: 13)
- `--ncpus`: Parallel processes (default: 8)

### snp_correlation.py

Calculates Pearson correlation between tissue and cfDNA VAF values.

```bash
python bin/snp_correlation.py \
  --tissue-pileup tissue_pileup.tsv.gz \
  --cfDNA-pileup cfDNA_pileup.tsv.gz \
  --min-snp-count 300 \
  --min-depth 10 \
  --output sample_prefix
```

#### snp_correlation.py Options

- `--tissue-pileup`: Tissue pileup file
- `--cfDNA-pileup`: cfDNA pileup file
- `--min-snp-count`: Minimum merged SNPs for QC pass (default: 300)
- `--min-depth`: Minimum sequencing depth (default: 10)
- `--output`: Output file prefix

### Python Dependencies

```text
pandas
numpy
scipy
click
pysam
rich
```

## Troubleshooting

### Common Issues

#### Issue: "No common SNP sites found"

- Check that both pileup files use the same reference genome build
- Verify that cfDNA pileup covers the target SNP regions
- Ensure sufficient sequencing depth in both samples

#### Issue: "Low correlation (<0.6)"

- Potential sample swap or contamination
- Verify sample identifiers in samplesheet
- Check tissue and cfDNA sample collection dates
- Review library preparation protocols

#### Issue: "Insufficient SNP count"

- Increase sequencing depth
- Check that known_sites_tsv covers your target regions
- Reduce `--min-depth` threshold (caution: may reduce accuracy)

## Pipeline Information

### Manifest

- **Pipeline Name**: snp_pair_qc
- **Author**: vanreality
- **Version**: dev1.0
- **Description**: A pipeline for SNP pair quality control

### Report Generation

The pipeline automatically generates execution reports:

- **Timeline**: Visual timeline of process execution
- **Report**: Resource usage and execution statistics
- **Trace**: Detailed task-level execution trace
- **DAG**: Directed acyclic graph of workflow

## Support

For issues, questions, or feature requests:

1. Check the [troubleshooting section](#troubleshooting)
2. Review Nextflow logs in `work/` directory
3. Open an issue on the project repository

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319.
- **BWA-MEM2**: Vasimuddin, M., et al. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. *IPDPS*.
- **Picard**: [Broad Institute Picard Tools](http://broadinstitute.github.io/picard/)
- **nf-core modules**: Ewels, P.A., et al. (2020). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology*, 38, 276–278.

## License

This project is available under standard open-source licensing terms (see LICENSE file if available).

---

**Developed by**: vanreality
**Last Updated**: 2025
