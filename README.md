# MUSA: Macaque Unified Set of Alleles

MUSA (Macaque Unified Set of Alleles) is a repository for the alleles curated from Rhesus macaque.

## Repository Structure

| Folder         | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `MUSA/`        | Tables of the curated alleles in the study.                                 |
| `SeedSet/`     | Initial seed set of alleles used to align the genomic and AIRR-seq data.    |
| `scripts/`     | Utility scripts for processing and analysis.                                |
| `data/`        | Additional data required for the analysis and figures.                      |
| `figures/`     | Figures presented in the manuscript.                                        |

## Main Files

| File | Description |
|-------|-------------|
| `MUSA/MUSA_2025-02-05.csv` | Latest version of the MUSA set for all seven Ig loci, containing information on curated and externally sourced alleles, with added details on RSS and Leader variation, including summary information. |

### `MUSA_2025-02-05.csv` Columns

| Column Name | Description |
|-------------|-------------|
| `allele` | The allele name. |
| `iglabel` | The unique coding sequence identifier. |
| `gene_type` | The gene type (e.g., IGHV). |
| `chain` | IG chain (e.g., IGH, IGL, IGK). |
| `seq` | Sequence of the core coding region. |
| `seq_gapped` | IMGT-gapped sequence of the core coding region (V-alleles only). |
| `functional` | Functional status of the allele. |
| `{segment}_heptamer` | Heptamer sequence of the RSS for the segment (e.g., `v_heptamer`). |
| `{segment}_nonamer` | Nonamer sequence of the RSS for the segment (e.g., `v_nonamer`). |
| `spacer_{location}` | Spacer sequence of the RSS for the segment (e.g., `spacer_3`). |
| `l_part1` | Leader part 1 sequence (V-genes only). |
| `l_part2` | Leader part 2 sequence (V-genes only). |
| `spacer_{location}_size` | Length of the Spacer sequence for the segment (e.g., `spacer_3_size`). |
| `sample_count_genomic` | Number of subjects where the allele is detected in genomic data. |
| `samples_genomic` | List of subjects where the allele is detected in genomic data. |
| `sample_count_AIRRseq` | Number of subjects where the allele is detected in AIRR-seq data. |
| `samples_AIRRseq` | List of subjects where the allele is detected in AIRR-seq data. |
| `sure_subject_count` | Number of subjects where the allele appeared with a single RSS variation. |
| `sure_subjects` | List of subjects with a single RSS variation for the allele. |
| `non_sure_subject_count` | Number of subjects where the allele appeared with multiple RSS variations. |
| `non_sure_subjects` | List of subjects with multiple RSS variations for the allele. |
| `pubid`, `genebank`, `kimdb`, `rhgldb+`, `imgt`, `trios`, `vrc`, `guo`, `novel` | Source names of the allele. |
| `notes`, `novel_notes`, `digger_notes`, `rss_notes` | Additional notes on the allele. |

## Usage

To get started, clone this repository:

```bash
git clone https://github.com/ayeletperes/musa.git
cd musa
