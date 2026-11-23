# Containers

Two images are built from this directory; all other tools are Biocontainers pulled by
Nextflow (see `conf/containers.config`).

| Image | Recipe | Contents |
|---|---|---|
| `rlang.sif` | `rlang/rlang.def` | R 4.6.1, Bioconductor 3.23, CRAN snapshot 2026-11-20, mumerge@35c8afb, bedtools, python3, UCSC bedClip/bedToBigBed v479 |
| `chipseq_tools.sif` | `chipseq_tools/chipseq_tools.def` + `environment.yml` | bowtie2, samtools, sambamba, bedtools, macs2, epic2, UCSC bedClip/bedToBigBed/bedGraphToBigWig, diffreps-rs 0.1.0 |

## Build on SCC

```bash
ssh scc-i01                      # build node
mkdir -p /scratch/$USER && cd /scratch/$USER
cp -r /path/to/chipseq_nextflow/containers .

cd containers/rlang
singularity build --fakeroot rlang.sif rlang.def                    # ~10-15 min
cd ../chipseq_tools
singularity build --fakeroot chipseq_tools.sif chipseq_tools.def    # ~5 min

# %test runs automatically at the end of each build; re-run any time with:
singularity test ../rlang/rlang.sif
singularity test chipseq_tools.sif

mkdir -p /projectnb/wax-es/routines/singularity/chipseq
mv ../rlang/rlang.sif chipseq_tools.sif /projectnb/wax-es/routines/singularity/chipseq/
```

Check what is inside an image:

```bash
singularity exec rlang.sif cat /opt/r_packages.tsv
singularity exec chipseq_tools.sif cat /opt/environment.lock.yml
```

## Which process uses which image

| Image | Processes |
|---|---|
| `rlang_sif` | all processes running R scripts (incl. `xlsx2track`, `sample_specific_tracks`) |
| `tools_sif` | `bowtie2_align`, `bam_count`, `macs2_callpeak`, `epic2_callpeak`, `create_bigwig_files`, `bed2track` (x2), `total_union`, `fragments_union_overlap_stats`, `diffreps` |
| `containers.trimgalore` | `trim_adapters` |
| `containers.fastqc` | `fastqc` |
| `containers.multiqc` | `multiqc` |
| `containers.picard` | `collect_metrics` |
| `containers.subread` | `calc_coverage`, `COUNT_SEX_MARKERS` |
| `containers.deeptools` | `combine_bigwig_tracks` |
| `containers.poppler` | `combine_mn2_dr_pdfs`, `combine_noxy`, `combine_allchr`, `combine_histograms`, `combine_pdf` |

## Updating

- R packages: change the CRAN snapshot date (and/or the Bioconductor base tag) in `rlang.def`, rebuild.
- Tools: change versions in `environment.yml`; for a new diffreps-rs release change
  `DIFFREPS_VERSION` and `DIFFREPS_SHA256` in `chipseq_tools.def`.
- Biocontainers: change the tag in `conf/containers.config`.
