# Example data

Five RAD loci (plus five monomorphic loci) from 109 Yosemite toads (*Anaxyrus canorus*)
in 12 populations, from Maier et al. (2019).

| Directory    | Files | Format |
|--------------|-------|--------|
| `stacks_v1/` | `batch_1.fa`, `batch_1.vcf`, `popmap_legacy.tsv`, `whitelist.txt` | Stacks 1.19 output and the 3-column population file used by fasta2genotype 1.x |
| `stacks_v2/` | `populations.samples.fa`, `populations.snps.vcf`, `popmap.tsv`, `whitelist.txt` | The same data re-encoded as Stacks 2 `populations` output (sample names in the FASTA headers, homozygotes written as two alleles, two-column popmap) |

Both directories give identical results, which the test suite checks.

```sh
# Stacks 2 style input, every output format
python3 fasta2genotype.py -f examples/stacks_v2/populations.samples.fa \
    -p examples/stacks_v2/popmap.tsv -o example_out/v2 -F all --remove-monomorphic

# Stacks 1 style input with a read-depth filter
python3 fasta2genotype.py -f examples/stacks_v1/batch_1.fa \
    -p examples/stacks_v1/popmap_legacy.tsv --vcf examples/stacks_v1/batch_1.vcf \
    --min-coverage 20 -o example_out/v1 -F structure genepop

# The fasta2genotype 1.x interactive interface
python3 fasta2genotype.py examples/stacks_v1/batch_1.fa examples/stacks_v1/whitelist.txt \
    examples/stacks_v1/popmap_legacy.tsv examples/stacks_v1/batch_1.vcf example_out/legacy
```
