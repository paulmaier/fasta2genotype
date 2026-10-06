# Changelog

## 2.0.0 (2026-10-05)

A rewrite for Python 3 and current Stacks releases. The old interactive interface still
works (see the README). The changes below can alter results compared with 1.x, so read
them before comparing new output with old analyses.

### New

- Reads Stacks 2 output (`populations.samples.fa`, the two-column popmap, and
  `populations.snps.vcf`) as well as Stacks 1 files. The format is detected automatically.
- Full command-line interface (`-f`, `-p`, `-o`, `-F`, and so on). Several formats can be
  written in one run (`-F structure genepop`, or `-F all`).
- Blacklists (`-b`), gzipped input, and loci written in whitelist order when a whitelist
  is given.
- New outputs: `diyabc-snp` (DIYABC-RF SNP format). `lfmm` now writes LEA's `.lfmm`
  genotype matrix.
- Missing bases (`N`) in Stacks 2 haplotypes are handled; see "Missing bases" in the
  README. Use `--keep-ambiguous` for the old behavior.
- `--hwe-alpha` (previously fixed at 0.05), `--coverage-stat`, and `--one-snp` for LFMM
  and DIYABC-RF as well as TreeMix.
- No dependencies: numpy and scipy are no longer needed.
- Installable from PyPI (`pip install fasta2genotype`), which adds a `fasta2genotype`
  command.
- Much faster. Genotypes are stored as integer indices of unique haplotypes rather than in
  nested dictionaries that were searched repeatedly.

### Changed behavior and bug fixes

Filters:

- **Allele-frequency filter**: homozygotes now count as two gene copies, and the frequency
  is relative to the gene copies genotyped at the locus. Previously homozygotes counted
  once and the denominator was twice the total number of individuals.
- **Allele-population filter**: now compares the fraction of populations that carry the
  haplotype with the threshold. 1.x divided by the gene-copy count of one arbitrary
  population instead.
- **Heterozygosity filter**: applied as documented. A locus is removed if its overall
  observed heterozygosity is ≥ the cutoff, exceeds the HWE expectation, *and* fails the HWE
  test. 1.x flagged a locus when any single heterozygous genotype was more common than
  expected. That also removed loci with an overall heterozygote *deficit*, such as those
  showing a Wahlund effect.
- **Individual missing-data filter**: 1.x stopped with an error (`remove_inds`); it now
  works. Removed individuals are also excluded from sample sizes and gene-copy counts.
- Loci, individuals and populations left with no data after filtering are removed. 1.x
  removed individuals only when their data were lost while reading the input (for
  example, to coverage filtering), and otherwise wrote them as missing.
- The monomorphic-locus check runs again after the other filters.
- FASTA samples that are not in the population map are ignored everywhere. 1.x still used
  them in the filter calculations.

Output formats:

- **TreeMix**: homozygotes now count as two gene copies (1.x counted them once). Loci
  whose first SNP is not biallelic are no longer skipped with `--one-snp`. Output is
  gzipped, as TreeMix requires.
- **LFMM**: genotypes come from the filtered haplotypes rather than directly from the VCF,
  so the filters now apply and no VCF or whitelist is needed. 1.x wrote a PLINK-style file
  with base letters and did not apply the filters. It also coded the first individual
  differently from the rest.
- **G-PhoCS**: individuals missing at a locus are now omitted and the per-locus sample
  count is correct. 1.x wrote their line without a line break, so it ran into the next
  line.
- **Structure**: populations are integers (with a `.pops.tsv` key) and missing data are
  `-9`, as Structure expects. 1.x wrote population names and used 0 for missing data.
- **samBada**: each gene-copy row now carries exactly one haplotype, as the 1.x manual
  shows. For heterozygotes, 1.x put both haplotypes in the first row and none in the
  second.
- **Arlequin**: the project title is written as `Title="..."`. Sample names no longer get a
  `Pop_` prefix, and `MissingData='?'`.
- **DIYABC**: the header uses the sex-ratio tag `<NM=1.0NF>`.
- **Genepop**: genotypes are separated by spaces, and an error is raised if a locus has too
  many alleles for 4-digit coding.
- **BayeScan**: the number of alleles per locus no longer includes a missing-data
  pseudo-allele.
- **Haplotype integer codes** are numbered by first appearance in output order, so the
  numbers can differ from 1.x. The genotypes are the same.
- **migrate-n**: names longer than 9 characters are replaced with short codes, with a
  lookup table, instead of stopping with an error. PHYLIP accepts long names (relaxed
  PHYLIP).
- **PHYLIP**: the row count in the header matches the rows written. SNP columns are found
  from the haplotypes present after filtering, ignoring `N`. For the `pi` option,
  homozygous individuals in haploid mode count as two gene-copy rows.
- IUPAC codes ignore missing bases (A + N → A rather than N).
- Individuals and populations are written in population-map order (1.x sorted them as
  text, so `85` came after `536`). Loci are sorted numerically (1.x: as text) unless a
  whitelist gives the order.
- Output files have format-specific extensions (`.str`, `.gen`, `.arp`, `.phy`, ...).
  Old-style interactive runs keep the old `.out` names.

### Coverage filtering

- The VCF is parsed by FORMAT field name, so Stacks 2 VCFs (`GT:DP:AD:GQ:GL`) work.
- A locus's depth for a sample is combined over all its SNPs (`--coverage-stat`, default
  `mean`). In Stacks 1 all SNPs of a locus had the same depth, so this matches 1.x there.

## 1.10 (2017)

Last Python 2 release (git tag `v1.10`).
