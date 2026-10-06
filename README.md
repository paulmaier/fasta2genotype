# fasta2genotype

`fasta2genotype.py` converts the per-sample haplotype FASTA file written by the
[Stacks](https://catchenlab.life.illinois.edu/stacks/) `populations` program into input
files for population-genetic software. It can apply quality filters along the way.

It reads Stacks 2 (`populations.samples.fa`) and Stacks 1 (`batch_N.fa`) files, and it
writes these formats:

| `-F` name      | Program / format                                      | Data used                  |
|----------------|-------------------------------------------------------|----------------------------|
| `migrate`      | [migrate-n](https://peterbeerli.com/migrate-html5/)   | DNA sequences              |
| `arlequin`     | [Arlequin](http://cmpg.unibe.ch/software/arlequin35/) | DNA sequences              |
| `diyabc`       | DIYABC (sequence loci, genepop-like format)           | DNA sequences              |
| `gphocs`       | [G-PhoCS](http://compgen.cshl.edu/GPhoCS/)            | DNA sequences (IUPAC)      |
| `phylip`       | PHYLIP alignment (RAxML, IQ-TREE, SVDquartets, ...)   | SNPs or full sequences     |
| `treemix`      | [TreeMix](https://bitbucket.org/nygcresearch/treemix) | biallelic SNPs             |
| `lfmm`         | LFMM / [LEA](https://bioconductor.org/packages/LEA/)  | biallelic SNPs (0/1/2)     |
| `diyabc-snp`   | [DIYABC-RF](https://diyabc.github.io/) SNP format     | biallelic SNPs (0/1/2)     |
| `structure`    | [Structure](https://web.stanford.edu/group/pritchardlab/structure.html) | haplotypes as integers |
| `genepop`      | [Genepop](https://genepop.curtin.edu.au/)             | haplotypes as integers     |
| `arlequin-hap` | Arlequin (STANDARD data type)                         | haplotypes as integers     |
| `genalex`      | [GenAlEx](https://biology-assets.anu.edu.au/GenAlEx/) | haplotypes as integers     |
| `bayescan`     | [BayeScan](http://cmpg.unibe.ch/software/BayeScan/)   | haplotype counts           |
| `sambada`      | [samBada](https://github.com/Sylvie/sambada)          | haplotype presence/absence |
| `allele-freq`  | Table of haplotype frequencies per population         | haplotype frequencies      |

In the "haplotypes as integers" formats, each unique sequence at a locus is treated as
one allele. Short RAD loci then act like multi-allelic markers, similar to microsatellites.

## Contents

- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Input files](#input-files)
  - [Haplotype FASTA](#haplotype-fasta--f-required)
  - [Population map](#population-map--p-required)
  - [Whitelist / blacklist](#whitelist--blacklist--w--b-optional)
  - [VCF](#vcf---vcf-optional)
- [Filters](#filters)
  - [Cut sites and adapters](#cut-sites-and-adapters)
  - [Missing bases (N)](#missing-bases-n)
- [Format notes](#format-notes)
  - [PHYLIP options](#phylip-options)
- [Old-style interactive mode](#old-style-interactive-mode)
- [Changes from version 1.x](#changes-from-version-1x)
- [Testing](#testing)
- [Citation](#citation)
- [License](#license)

## Requirements

Python 3.9 or newer. There are no other dependencies: numpy and scipy are no longer needed.

Install it from PyPI, which adds a `fasta2genotype` command:

```sh
pip install fasta2genotype
fasta2genotype --help
```

Or run the script from a copy of this repository:

```sh
git clone https://github.com/paulmaier/fasta2genotype.git
python3 fasta2genotype/fasta2genotype.py --help
```

`python3 fasta2genotype.py` only finds the script in the current directory; Python does
not search your `PATH`. To run a downloaded copy from anywhere, add the `fasta2genotype`
folder to your `PATH` and run the script by name, `fasta2genotype.py --help`
(macOS/Linux; the script is executable).

The examples below use `python3 fasta2genotype.py`, run from the repository folder. With
the installed command, type `fasta2genotype` instead.

## Quick start

With Stacks 2, run `populations` with `--fasta-samples` (and `--vcf` if you want coverage
filtering). Then run:

```sh
python3 fasta2genotype.py \
    -f populations.samples.fa \
    -p popmap.tsv \
    -o results/myproject \
    -F structure genepop phylip \
    --remove-monomorphic --min-locus-freq 0.8 --min-ind-freq 0.5
```

This writes `results/myproject.str`, `results/myproject.gen` and `results/myproject.phy`,
plus small companion tables such as `myproject.str.pops.tsv`. Use `-F all` to write
every format. Progress messages go to stderr; `-q` limits them to warnings.

Try it on the bundled examples:

```sh
python3 fasta2genotype.py -f examples/stacks_v2/populations.samples.fa \
    -p examples/stacks_v2/popmap.tsv -o example_out/demo -F all --remove-monomorphic
```

## Input files

### Haplotype FASTA (`-f`, required)

Stacks writes one record per sample, locus and gene copy:

```
# Stacks version 2.66; ...
>CLocus_198980_Sample_48_Locus_198980_Allele_0 [MVZ-142736; Bufo_bufo_chr01, 1104705, +]
TGCAGGGAGCCCTGTG...
>CLocus_198980_Sample_48_Locus_198980_Allele_1 [MVZ-142736; Bufo_bufo_chr01, 1104705, +]
TGCAGGGAGCCCTGTG...
```

* Stacks 2: use `populations.samples.fa`. Samples are matched to the population map by the
  name in brackets. Homozygotes appear as two identical alleles.
* Stacks 1 (`batch_N.fa`, e.g. v1.12–1.48): there are no bracketed names, so samples are
  matched by the number after `Sample_`. Homozygotes usually have only `Allele_0`.
* Only alleles 0 and 1 are used, because diploids are assumed. Sequences with a higher
  allele number are ignored and counted in a warning.
* Gzipped files (`.gz`) are read directly.

### Population map (`-p`, required)

Either of two formats is detected automatically. Override with `--popmap-format`.

**Stacks popmap** (the file given to `populations -M`): two tab-separated columns and
no header. Sample names are used in the output files.

```
MVZ-142736	ANBO
BLSU22-0001	ANBO
DEVA23-0001	ANBOxCA
```

**Legacy 3-column file** (fasta2genotype 1.x): a header row, then
`SampleID`, `IndividualID` and `PopulationID`. `SampleID` is matched to the FASTA headers
and `IndividualID` is used in the output. This is the format to use with Stacks 1 files.

```
SampleID	IndividualID	PopulationID
561	S11-0074	M101
577	S11-0057	M101
```

Individuals and populations are written in the order they appear in the population map.
Samples in the FASTA file but not in the map are skipped.

### Whitelist / blacklist (`-w`, `-b`, optional)

One catalog locus ID per line, which is the number after `CLocus_`. Extra columns, such as
the SNP column in a Stacks 2 whitelist, are ignored, and whole loci are kept or removed.
When a whitelist is given, loci are written in whitelist order; otherwise they are
written in numeric order.

### VCF (`--vcf`, optional)

Needed only for `--min-coverage`. Use the SNP VCF from the same `populations` run
(`populations.snps.vcf` in Stacks 2, `batch_N.vcf` in Stacks 1). The read depth (`DP`) of
each sample is combined across the SNPs of each locus using `--coverage-stat`, which is
`mean` by default. The alternatives are `min` and `max`.

The VCF only lists variable loci. Loci that are absent from it (usually monomorphic ones)
are therefore removed when coverage filtering is used.

## Filters

Every filter is off by default. They run in the order listed below.

| Option                     | Effect |
|----------------------------|--------|
| `--min-coverage N`         | Set genotypes with read depth < N to missing (needs `--vcf`). |
| `--remove-monomorphic`     | Remove loci with a single haplotype. Checked again after all other filters. |
| `--het-cutoff H`           | Remove likely paralogs: loci whose observed heterozygosity is ≥ H **and** above the Hardy–Weinberg expectation, with a significant chi-square test of HWE genotype proportions (`--hwe-alpha`, default 0.05). Loci with a heterozygote *deficit* are never removed. |
| `--min-allele-freq F`      | Set genotypes carrying a haplotype with overall frequency < F to missing. |
| `--min-allele-pops F`      | Set genotypes carrying a haplotype found in < F of populations to missing (e.g. 0.25 = 4 of 16 populations). |
| `--min-pop-locus-freq F`   | Within each population, remove a locus if < F of that population's individuals are genotyped. |
| `--min-locus-freq F`       | Remove loci genotyped in < F of all individuals. |
| `--min-ind-freq F`         | Remove individuals genotyped at < F of the remaining loci. |

After filtering, empty loci, empty individuals and empty populations are removed. A
summary of what each step removed is printed.

### Cut sites and adapters

`--clip-5prime` and `--clip-3prime` remove a sequence from the start or end of every
haplotype that carries it. If you give several sequences, only the first match is
removed. For example, with double-digest data where read 1 starts with `TGCAGG` and read 2
with `CGG`:

```sh
--clip-5prime TGCAGG CGG
```

Check your FASTA first: sequences can be reversed or reverse-complemented. A locus whose
haplotypes end up with different lengths is removed with a warning.

### Missing bases (N)

Stacks 2 writes `N` at SNP positions it could not call in a sample. Such haplotypes would
look like new alleles. By default:

* a genotype with an `N` at a variable (SNP) column is set to missing;
* a column that is invariant apart from `N`s is set to `N` in every haplotype, so identical
  haplotypes stay identical;
* columns that are entirely `N` (e.g. the spacer between paired-end reads) are left as is.

In one Stacks 2.66 dataset this affected about 2.5% of genotypes. It also halved the
apparent number of haplotypes per locus. `--keep-ambiguous` turns the handling off and
treats `N` as an ordinary character, as version 1.x did.

## Format notes

* **SNP formats** (`treemix`, `lfmm`, `diyabc-snp`, and `phylip --phylip-sites snps`):
  SNPs are found by comparing the haplotypes that remain after filtering. The Stacks VCF is
  not used. Only biallelic SNPs are written to `treemix`, `lfmm` and `diyabc-snp`. Use
  `--one-snp` to keep only the first biallelic SNP per locus, for unlinked markers.
* **lfmm** writes `PREFIX.lfmm` (individuals × SNPs, minor-allele counts 0/1/2, 9 =
  missing), plus `PREFIX.lfmm.snps.tsv` (locus, 0-based column, major/minor base) and
  `PREFIX.lfmm.inds.tsv`. The `.lfmm` file can be used directly with LEA's `lfmm2()`;
  convert it with `lfmm2geno()` for `snmf()`.
* **treemix** output is gzipped (`PREFIX.treemix.gz`) because TreeMix reads gzipped input.
  Counts are gene copies of the alphabetically first and second base.
* **structure**: two rows per individual, missing data coded `-9`, and populations coded
  as integers in the second column (`PREFIX.str.pops.tsv` maps them to names). Use
  `MARKERNAMES=1`, `POPDATA=1`, `ONEROWPERIND=0` and `MISSING=-9`.
* **genepop**: `--genepop-digits 6` (default; up to 999 alleles per locus) or `4` (up to
  99). `PREFIX.gen.pops.tsv` lists the population order.
* **bayescan**: `PREFIX.bayescan.loci.tsv` and `PREFIX.bayescan.pops.tsv` map the numbers
  used in the file back to locus and population names.
* **migrate** and **arlequin**: missing sequences are written as `?`. migrate-n allows
  names of at most 10 characters (9 plus the `a`/`b` gene-copy suffix). If any name is
  longer, samples are renamed `I001`, `I002`, ... and `PREFIX.migrate.names.tsv` records the
  mapping.
* **gphocs** lists only the individuals genotyped at each locus. Heterozygous sites use
  IUPAC codes.
* **sambada**: two rows per individual (`IDa`, `IDb`), one per gene copy. Each row has a
  0/1 column per haplotype and `NaN` for missing data.
* Titles (`--title`, default: the output prefix) are used by Arlequin, DIYABC, Genepop and
  GenAlEx.

### PHYLIP options

| Option                   | Choices | Meaning |
|--------------------------|---------|---------|
| `--phylip-sites`         | `snps` (default), `full` | concatenate variable sites only, or complete sequences |
| `--phylip-mode`          | `haploid`, `diploid` (default), `population` | one row per gene copy, per individual, or per population; the last two use IUPAC ambiguity codes |
| `--phylip-loci`          | `all` (default), `pi`, `fixed` | keep only sites (with `snps`) or loci (with `full`) that have a fixed difference between rows: shared by at least two rows (`pi`, parsimony-informative) or anywhere (`fixed`) |
| `--phylip-breakpoints`   | | put `!` between loci |
| `--phylip-locus-header`  | | add a tab-separated row of locus names under the header |

Names are padded to at least 10 characters. Longer names are kept, in relaxed PHYLIP
style, which RAxML, IQ-TREE and most current software accept.

## Old-style interactive mode

The 1.x interface still works. Give the five positional arguments, using `NA` for an
unused whitelist or VCF, and answer the questions:

```sh
python3 fasta2genotype.py batch_1.fa whitelist.txt popmap_legacy.tsv NA MyOutput
```

The questions are the same as before, so scripts that pipe answers into the program keep
working. Output files keep their old names: `MyOutput.out`, `MyOutput_pops.out` and
`MyOutput_loci.out`.

## Changes from version 1.x

Version 2.0 is a rewrite. It keeps the formats and options of 1.x, fixes many errors, and
is much faster. On a Stacks 2 dataset of 631 loci and 877 samples, writing a Structure
file took 3 seconds instead of 2 minutes. Most notably, sample sizes, allele counts and population
denominators now follow the documented behavior. See [CHANGELOG.md](CHANGELOG.md) for
every change that can alter results. The last 1.x release is available in the git history
(tag `v1.10`).

## Testing

```sh
python3 -m unittest discover -s tests -v
```

The tests include a check that the Stacks 1 and Stacks 2 encodings of the example data
(`examples/stacks_v1`, `examples/stacks_v2`) give byte-identical output in every format.

## Citation

If you use fasta2genotype, please cite:

> Maier P.A., Vandergast A.G., Ostoja S.M., Aguilar A., Bohonak A.J. (2019). Pleistocene
> glacial cycles drove lineage diversification and fusion in the Yosemite toad
> (*Anaxyrus canorus*). *Evolution*, 73(12), 2476–2496.
> [https://doi.org/10.1111/evo.13868](https://www.doi.org/10.1111/evo.13868)
> ([PDF](https://paulmaierresearch.com/Maier_2019_Evolution.pdf))

Please also cite Stacks and the programs you use the output with.

## License

MIT; see [LICENSE](LICENSE).
