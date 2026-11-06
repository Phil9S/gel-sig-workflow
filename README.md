# gel-sig-workflow

`main.nf` is a Nextflow pipeline that processes the tumour/normal BAMs in `config/samples.tsv` to generate copy number profiles:

1. Battenberg clonal and subclonal copy number calls.
2. Extraction of Battenberg purity and ploidy for each sample.
3. Evaluation of the Battenberg profiles using somatic SNV VAF distributions on clonal segments (supplied VCFs, or Strelka2 calls).
4. Downsampled shallow WGS absolute copy number at the Battenberg purity/ploidy, using stage 2 of [swgs-absolutecn v2.0.0](https://github.com/Phil9S/swgs-absolutecn/releases/tag/v2.0.0) (precomputed ploidy/purity, run from the `phil9s/swgs-absolutecn:v2.0.0` container).
5. ASCAT at Affymetrix SNP6.0 loci (alleleCounter) at the Battenberg purity/ploidy.

Before the first run, set the reference paths in `conf/genomes.config` and build the Battenberg image (`containers/battenberg/`).

## Overview

Tumour/normal WGS BAMs → Battenberg copy number → purity/ploidy → three downstream uses
of that fit:

```
config/samples.tsv
      │
      ▼
 BATTENBERG ──► BATTENBERG_PURITY_PLOIDY ──► purity, ploidy, psi per sample
      │                    │
      │                    ├──► SWGS_ABSOLUTECN  downsample tumour BAM → QDNAseqmod → absolute CN
      │                    │                     (swgs-absolutecn v2.0.0 stage 2, precomputed ploidy/purity)
      │                    │
      │                    └──► ASCAT_SNP6       alleleCounter at SNP6.0 loci → ASCAT with
      │                                          rho/psi fixed to Battenberg (+ unconstrained fit)
      ▼
 SNV_VAF_QC   somatic SNVs (supplied VCF or Strelka2) on clonal Battenberg segments
              → expected vs observed VAF / CCF → PASS | REVIEW | INSUFFICIENT_SNVS
      │
      ▼
 CN_SUMMARY   one row per sample: summary/cn_pipeline_summary.tsv
```

The Nextflow pipeline lives alongside the Snakemake signature workflow (`workflow/`)
and reads the same `config/samples.tsv`.

## Requirements

- Nextflow ≥ 24.04 and Java 17+
- Singularity/Apptainer (HPC) or Docker
- Reference data (paths set in `conf/genomes.config`, see below)

### Containers

| Step | Default image | Notes |
|---|---|---|
| Battenberg | `phil9s/battenberg:3.0.0` | **build and push first** from `containers/battenberg/`: pixi environment (R 4.4.3, alleleCounter 4.3.0, IMPUTE2 2.3.2, OpenJDK 11) plus Battenberg v3.0.0, ASCAT v3.2.0 and copynumber from GitHub |
| ASCAT, alleleCounter | `phil9s/ascat_v3_1_2:latest` | the image used by the Snakemake workflow; contains `/ref/<build>/` G1000 references |
| sWGS | `phil9s/swgs-absolutecn:v2.0.0` | swgs-absolutecn v2.0.0 release: R 4.3.3, rswgsabsolutecn v1.1.1, QDNAseqmod, samtools; stage-2 scripts in `/app/workflow/scripts` |
| Strelka2 | `quay.io/biocontainers/strelka:2.9.10--h9ee0642_1` | |
| R utilities | `rocker/tidyverse:4.3.2` (to be replaced) | purity extraction, VAF evaluation, summary; needs R with optparse, vcfR, dplyr and ggplot2 |

```bash
docker build -t phil9s/battenberg:3.0.0 containers/battenberg
docker push phil9s/battenberg:3.0.0
```

The Battenberg image is built the same way as the swgs-absolutecn v2.0.0 image. A build stage
installs the environment locked in `containers/battenberg/pixi.lock` and then installs the GitHub R
packages into it. The slim `debian:12-slim` production stage contains only that environment, at
`/app/.pixi/envs/default`. After changing `pixi.toml`, run `pixi lock` in `containers/battenberg/`
and commit the updated `pixi.lock`; the build fails if the two do not match. Keep OpenJDK at 11:
Battenberg starts Beagle with `-XX:+UseParallelOldGC`, which Java 16 and later reject.
IMPUTE2 is only used with `--battenberg_use_beagle false`. Its binary is a statically linked 2015
build that needs the kernel's legacy vsyscall page, so it crashes (exit 139) on hosts booted with
`vsyscall=none`, such as WSL2. The stock bioconda IMPUTE2 container has the same limitation.

Override any image with `--container_battenberg`, `--container_ascat`, `--container_allelecounter`,
`--container_swgs`, `--container_strelka`, `--container_r`.

On systems without internet access (e.g. the GEL Research Environment) pull the images on a
connected machine and copy them into `--singularity_cache_dir`, using the file names Nextflow
expects (`<registry>-<image>-<tag>.img`, e.g. `phil9s-battenberg-3.0.0.img`):

```bash
singularity pull phil9s-battenberg-3.0.0.img docker://phil9s/battenberg:3.0.0
singularity pull phil9s-ascat_v3_1_2-latest.img docker://phil9s/ascat_v3_1_2:latest
singularity pull phil9s-swgs-absolutecn-v2.0.0.img docker://phil9s/swgs-absolutecn:v2.0.0
singularity pull quay.io-biocontainers-strelka-2.9.10--h9ee0642_1.img docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1
singularity pull rocker-tidyverse-4.3.2.img docker://rocker/tidyverse:4.3.2
```

### Reference data (`conf/genomes.config`)

Each sample's `build` column (`hg19`/`hg38`) selects a block in `params.genomes`. Edit the
file, or override values with `-c my_refs.config`:

```groovy
params {
    genomes {
        hg19 {
            battenberg_ref_dir = '/refs/battenberg/hg19'
            snp6_positions     = '/refs/ascat/GC_AffySNP6_102015.txt'
            fasta              = '/refs/GRCh37/genome.fa'
        }
    }
}
```

| Key | Used by | Description |
|---|---|---|
| `chr_prefix` | alleleCounter | `true` if BAM contigs are `chr1…` (default hg19 `false`, hg38 `true`) |
| `battenberg_ref_dir` | Battenberg | unpacked Battenberg bundle ([GRCh37](https://ora.ox.ac.uk/objects/uuid:2c1fec09-a504-49ab-9ce9-3f17bac531bc), [GRCh38 v3.0.0](https://doi.org/10.48420/30406441)) |
| `battenberg_*` | Battenberg | paths **relative to** `battenberg_ref_dir` (absolute paths also accepted). Defaults follow `inst/example/battenberg_wgs.R` (v3.0.0); check them against your unpacked bundle |
| `ascat_g1000_dir` | ASCAT | ASCAT 1000 Genomes WGS loci/alleles/GC/RT ([Zenodo](https://zenodo.org/records/14008443)); `null` uses `/ref/<build>/` inside the ASCAT container |
| `snp6_positions` | ASCAT | tab-separated SNP6.0 positions with chromosome/position columns, e.g. ASCAT's `GC_AffySNP6_102015.txt` (hg19). For hg38, supply lifted-over positions |
| `fasta`, `fasta_fai` | Strelka2 | only needed when SNVs are called by the pipeline |
| `strelka_call_regions` | Strelka2 | optional bgzipped BED (+ `.tbi`) restricting calling to the primary assembly |

`impute_info.txt` inside the Battenberg bundle may hold paths from the machine it was
built on. The wrapper rewrites them to the staged bundle (it matches the trailing path
components), so the bundle can be moved freely.

## Input: `config/samples.tsv`

| Column | Required | Description |
|---|---|---|
| `sample` | yes | unique sample ID (output prefix) |
| `tumour_bam`, `normal_bam` | yes | indexed BAMs (`.bam.bai` or `.bai`) |
| `tumour_name`, `normal_name` | yes | Battenberg/ASCAT sample names; also used to find the tumour column in supplied VCFs |
| `sex` | yes | `XX`/`XY` (`female`/`male` accepted) |
| `build` | yes | key in `params.genomes` (`hg19`, `hg38`) |
| `cancer` | yes | cancer type, used for output folders |
| `snv_vcf` | no | somatic SNV VCF (Strelka2, Mutect2 or any VCF with `AD` or `AF`+`DP`). Empty/`NA` → Strelka2 is run |
| `swgs_smooth` | no | swgs-absolutecn `smooth` flag: `TRUE` re-segments a noisy sWGS profile (default `--swgs_smooth false`) |

Other columns (e.g. `ploidy`/`purity` used by the Snakemake ASCAT rule) are ignored.

## Running

```bash
# SLURM
nextflow run main.nf -profile slurm,singularity \
    --queue general --long_queue long --cluster_account my_account \
    --singularity_cache_dir /path/to/sif_cache --container_bind /data,/refs

# LSF
nextflow run main.nf -profile lsf,singularity \
    --queue short --long_queue long --cluster_account re_gecip_cancer_pan

# PBS Torque / PBS Pro
nextflow run main.nf -profile pbs,singularity --queue batch --cluster_account my_project
nextflow run main.nf -profile pbspro,singularity --queue workq

# GEL Research Environment (LSF + Singularity, paths from config/config.yaml)
nextflow run main.nf -profile gel

# Local machine
nextflow run main.nf -profile local,docker
```

Add `-resume` to continue a run. Running the Nextflow head job inside a long, low-resource
scheduler job (or `screen`/`tmux` on a login node, where permitted) is recommended.

### Scheduler parameters

| Parameter | Default | Description |
|---|---|---|
| `--queue` | `null` | default queue (LSF/PBS) or partition (SLURM) |
| `--long_queue` | `null` | queue for tasks requesting more than `--long_queue_min_hours` (default 4) |
| `--cluster_account` | `null` | LSF `-P`, PBS `-A`, SLURM `--account` |
| `--slurm_qos` | `null` | SLURM `--qos` |
| `--cluster_options` | `null` | extra raw submission options |
| `--queue_size` | `100` | maximum concurrently submitted jobs |
| `--submit_rate_limit` | `50/1min` | job submission rate |
| `--lsf_per_job_mem_limit` / `--lsf_per_task_reserve` | `false` | LSF memory semantics, if not readable from `lsf.conf` |
| `--max_cpus` / `--max_memory` / `--max_time` | `32` / `256.GB` / `168.h` | ceilings applied to every task |
| `--singularity_cache_dir`, `--container_bind` | `null` | image cache and extra bind paths (comma separated) |

Per-process CPU/memory/time are in `conf/base.config`. Memory and time scale with the
retry attempt (up to 2 retries on scheduler kills). Set `--max_*` on the command line or
in `-params-file`. In a custom `-c` config, set `process.resourceLimits` directly, because
`conf/base.config` resolves the limits before custom configs are read.

### Pipeline parameters

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `config/samples.tsv` | sample sheet |
| `--outdir` | `results/cn_pipeline` | output directory |
| `--skip_vaf_evaluation`, `--skip_swgs`, `--skip_ascat` | `false` | disable a step |
| `--call_snvs` | `true` | run Strelka2 for samples without `snv_vcf` |
| `--ignore_failed_samples` | `true` | a sample whose task still fails after retries is dropped (reported in the trace) instead of stopping the cohort; cohort-level steps always stop the run |
| `--battenberg_*` | as `battenberg_wgs.R` v3.0.0 | e.g. `--battenberg_segmentation_gamma 20`, `--battenberg_max_rho 1.02`, `--battenberg_use_beagle true` |
| `--vaf_min_depth` / `--vaf_min_alt_reads` | `20` / `3` | SNV filters |
| `--vaf_ccf_tolerance` / `--vaf_max_exceed_fraction` | `0.1` / `0.05` | PASS thresholds |
| `--swgs_bin` | `30` | swgs-absolutecn `bins`: one bin size (kb) per run, from 1, 5, 15, 30, 50, 100, 500, 1000 |
| `--swgs_project` | `gel_cn` | swgs-absolutecn `project_name`, prefix of the sWGS outputs |
| `--swgs_use_seed` / `--swgs_seed` | `true` / `9999` | `use_seed` / `seed_val` for CBS segmentation |
| `--swgs_paired_end` | `true` | `pairedEnd`: set `false` for single-end reads |
| `--swgs_auto_smooth` / `--swgs_smooth_threshold` | `false` / `300` | `autoSmooth` / `smoothThreshold`: re-segment profiles with more segments than the threshold |
| `--swgs_publish_bams` | `false` | also publish the downsampled BAMs |
| `--ascat_penalty` / `--ascat_gamma` | `70` / `1` | ASCAT segmentation penalty / platform gamma |
| `--ascat_run_unconstrained` | `true` | also report ASCAT's own purity/ploidy |

## Outputs (`--outdir`)

```
battenberg/<cancer>/<sample>/        <tumour>_copynumber.txt (clonal + subclonal, solution A)
                                     <tumour>_copynumber_extended.txt (solutions A-F)
                                     <tumour>_purity_ploidy.txt, <tumour>_rho_and_psi.txt, *.png, log
purity_ploidy/<cancer>/              <sample>.battenberg_purity_ploidy.tsv
snvs/strelka/<cancer>/<sample>/      Strelka2 somatic SNV/indel VCFs (when called here)
vaf_evaluation/<cancer>/<sample>/    .metrics.tsv, .states.tsv, .snvs.tsv.gz, .pdf
sWGS_fitting/<project>_<bin>kb/      swgs-absolutecn v2.0.0 stage 2 layout:
  absolute_PRE_down_sampling/          <project>_<bin>kb_combined_QC_predownsample.tsv (precomputed fits)
  absolute_POST_down_sampling/         <project>_<bin>kb_ds_abs_fits.tsv, _ds_absCopyNumber.rds,
                                       _ds_absCopyNumber_segTable.tsv (cohort)
    <sample>/                          <project>_<sample>_fitted_QC.tsv
    <sample>/relative_cn_rds/          <project>_<sample>_<bin>kb_relSmoothedCN.rds
    <sample>/abs_cn_rds/               <project>_<sample>_<bin>kb_ds_abs_fits.tsv, _ds_absCopyNumber.rds,
                                       _ds_absCopyNumber_segTable.tsv, plots/<sample>.png
    downsampled_bams/                  <sample>.bam(.bai) with --swgs_publish_bams
ascat_snp6/reference/ascat_snp6_<build>/  SNP6-restricted loci, alleles, GC, RT
ascat_snp6/<cancer>/<sample>/        segments, segments_raw, metrics, RData, plots
summary/                             battenberg_purity_ploidy.tsv, vaf_evaluation_metrics.tsv,
                                     ascat_snp6_metrics.tsv, cn_pipeline_summary.tsv
pipeline_info/                       execution report, timeline, trace, DAG
```

## Method notes

### Purity and ploidy

`BATTENBERG_PURITY_PLOIDY` reads `<tumour>_purity_ploidy.txt` and falls back to the
`FRAC_GENOME` row of `<tumour>_rho_and_psi.txt` (older Battenberg versions):

- `purity`: rho of the FRAC_GENOME solution
- `ploidy`: segment-length-weighted tumour ploidy of the final profile (subclonal states weighted by their cell fraction)
- `psi`: tumour ploidy parameter of the fit (`psit`), the quantity ASCAT expects in `psi_manual`

sWGS absolute copy number uses `purity` and `ploidy`. ASCAT uses `purity` (`rho_manual`) and
`psi` (`psi_manual`). Samples without a usable Battenberg fit are skipped downstream with a
warning.

### VAF evaluation of clonal segments

For a clonal SNV on `m` copies in a segment with total copy number `N = nMaj + nMin`, the
expected VAF is `m·p / (p·N + 2(1−p))`, for `m = 1 … nMaj`. Using PASS SNVs (depth ≥ 20, ≥ 3 alt
reads, autosomes) on clonal Battenberg segments (`frac1_A = 1`), the evaluation reports:

- **Per CN state** (`.states.tsv`): expected and observed VAF peaks (KDE), the rightmost
  observed peak, its nearest multiplicity, and the purity implied by that peak,
  `p = 2·VAF / (m − VAF·(N − 2))`.
- **Per sample** (`.metrics.tsv`):
  - `ccf_clonal_peak` is the clonal CCF cluster and should be ≈ 1.
  - `frac_exceed_max` is the fraction of SNVs significantly above `VAF(nMaj)`; this happens when purity or copy number is too low.
  - `snv_implied_purity` and `purity_delta` give the SNV-based purity, independent of the Battenberg value.
  - `frac_multiplicity_gt1`: very low values with many `nMaj ≥ 2` segments can indicate a spurious genome-doubled solution.
- **`qc_flag`**:
  - `PASS`: `|ccf_clonal_peak − 1| ≤ 0.1` and `frac_exceed_max ≤ 5%`.
  - `REVIEW`: otherwise, with the reason given in `qc_reason`.
  - `INSUFFICIENT_SNVS`: fewer than 100 SNVs on clonal segments.
- **Plots** (`.pdf`):
  1. VAF histograms per CN state, with the expected peaks marked.
  2. Pooled CCF.
  3. Genome-wide VAF, with the expected clonal VAF for each segment.

On simulated data at purity 0.6, a correct fit gives `PASS` (CCF peak 0.99, SNV-implied
purity 0.598). Fits at purity 0.85 and 0.4 are flagged `REVIEW` (CCF peaks 0.65 and 1.39).

### Downsampled sWGS (swgs-absolutecn v2.0.0)

This step runs stage 2 of the [swgs-absolutecn v2.0.0 release](https://github.com/Phil9S/swgs-absolutecn/releases/tag/v2.0.0)
with precomputed ploidy/purity, here the Battenberg `purity` and `ploidy`. The release's own
stage-2 scripts are run in their command-line mode from the `phil9s/swgs-absolutecn:v2.0.0`
image, where they sit in `/app/workflow/scripts` (identical to the v2.0.0 tag). Their R code
is in rswgsabsolutecn v1.1.1.

| swgs-absolutecn v2.0.0 | Nextflow process |
|---|---|
| `processPrecomputed.R` + `splitInput.R` | `SWGS_DOWNSAMPLE_QDNASEQ`: `bin/swgs_precomputed_qc.R`, the per-sample equivalent (same `rswgsabsolutecn::bamCheck`, `getDownsampleDepth` and QC columns) |
| `downsampleBams.R` | `SWGS_DOWNSAMPLE_QDNASEQ` |
| `qdnaseqMod.R` | `SWGS_DOWNSAMPLE_QDNASEQ` |
| `qdnaseqRelToAbs.R` | `SWGS_REL_TO_ABS` |
| `combineResults.R` | `SWGS_COMBINE` |

1. **Precomputed fit.** `downsample_depth = rswgsabsolutecn::getDownsampleDepth(ploidy, purity, nbins)`, where `nbins = 0.932 × usable bins` for the bin size.
2. **Downsampling.** The BAM is downsampled with `samtools view --subsample <downsample_depth / primary mapped reads> --subsample-seed 99`. Ratios above 0.96 keep the full BAM.
3. **Relative copy number.** QDNAseqmod binning (`pairedEnds`), filtering, GC/mappability correction, outlier smoothing, seeded CBS segmentation, then `rswgsabsolutecn::smoothSample` (the per-sample `smooth` flag, or `autoSmooth`).
4. **Absolute copy number.** `rswgsabsolutecn::gridStats` at the fixed ploidy/purity. Fit metrics are reported: segments, clonality, RMSE and segment variance (`.post` columns), plus TP53 copy number and expected allele fraction.
5. **Cohort outputs.** `combineResults.R` writes the cohort `.rds`, segment table and fit table.

Differences from running swgs-absolutecn directly:
- Downsampling and QDNAseq run in one task. A sample kept at full depth is a symlink to the input BAM, which is only guaranteed to be mounted in that task.
- rswgsabsolutecn requires `0 < purity ≤ 1`. Battenberg purities above 1 (`--battenberg_max_rho` defaults to 1.02) are capped at 1 for this step.
- `PATIENT_ID` is the sample ID, and `TP53freq` is `NA` because `samples.tsv` has no observed TP53 allele frequency.

### ASCAT at SNP6.0 loci

1. **SNP6 reference (once per build).** The SNP6.0 positions are intersected with ASCAT's 1000 Genomes WGS loci/alleles. GC and replication-timing files are subset to the same loci and keyed `chr_pos`, so `ascat.correctLogR` matches them. A per-chromosome retention summary is written next to the reference.
2. **alleleCounter, once per unique BAM.** A normal shared by several tumours is counted once. Loci names follow the BAM contig style (`chr_prefix`).
3. **ASCAT.** `ascat.getBAFsAndLogRs` → `correctLogR` → `aspcf` → `runAscat`, with `rho_manual = purity` and `psi_manual = psi` from Battenberg. ASCAT's unconstrained fit is also run, and its purity/ploidy/goodness of fit are reported for comparison (`--ascat_run_unconstrained`).

## Testing

`-stub-run` exercises the full DAG with placeholder outputs (each process has a `stub:`
block). The BAMs, their indexes and the configured reference paths (`battenberg_ref_dir`,
`snp6_positions`, `fasta`) must exist, but they can be empty placeholder files:

```bash
nextflow run main.nf -profile local -stub-run --samples test_samples.tsv --outdir stub_test
```

## Authors

* Philip Smith (@phil9s)
