#!/usr/bin/env python3
"""
fasta2genotype.py -- convert Stacks haplotype FASTA files into population-genetic formats.

Reads the per-sample haplotype FASTA written by the Stacks ``populations`` program
(``populations.samples.fa`` in Stacks 2.x, ``batch_N.fa`` in Stacks 1.x), applies
optional quality filters, and writes one or more of: migrate-n, Arlequin, DIYABC,
LFMM, PHYLIP, G-PhoCS, TreeMix, Structure, Genepop, allele frequencies, samBada,
BayeScan and GenAlEx input files.

Run ``python3 fasta2genotype.py --help`` for usage.

(c) Paul Maier. Released under the MIT License.
Requires Python 3.9+ and no third-party packages.
"""

import argparse
import collections
import gzip
import math
import os
import re
import sys
from array import array

__version__ = "2.0.0"

CITATION = (
    "Maier P.A., Vandergast A.G., Ostoja S.M., Aguilar A., Bohonak A.J. (2019).\n"
    "Pleistocene glacial cycles drove lineage diversification and fusion in the\n"
    "Yosemite toad (Anaxyrus canorus). Evolution, 73(12), 2476-2496.\n"
    "https://doi.org/10.1111/evo.13868"
)

ACGT = frozenset("ACGT")

# IUPAC nucleotide codes -> the set of bases they represent. Missing-data
# characters (N, ?, -, .) represent no information.
IUPAC_BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "", "?": "", "-": "", ".": "",
}
IUPAC_SETS = {code: frozenset(bases) for code, bases in IUPAC_BASES.items()}
BASES_TO_IUPAC = {frozenset(b): c for c, b in IUPAC_BASES.items() if c not in "?-."}

HEADER_RE = re.compile(
    r"^>CLocus_([^_\s]+)_Sample_([^_\s]+)_Locus_[^_\s]+_Allele_([^_\s]+)"
    r"(?:\s+\[\s*([^;\]]+?)\s*[;\]])?"
)


class InputError(Exception):
    """A problem with the user's input files or options."""


# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #

QUIET = False


def info(msg):
    if not QUIET:
        print(msg, file=sys.stderr)


def warn(msg):
    print("WARNING: " + msg, file=sys.stderr)


def preview(items, n=5):
    items = list(items)
    text = ", ".join(str(i) for i in items[:n])
    return text + (", ..." if len(items) > n else "")


# --------------------------------------------------------------------------- #
# Small utilities
# --------------------------------------------------------------------------- #

def open_text(path):
    try:
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "r")
    except OSError as e:
        raise InputError("cannot open '%s': %s" % (path, e.strerror or e))


def open_out(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    if path.endswith(".gz"):
        return gzip.open(path, "wt", newline="\n")
    return open(path, "w", newline="\n")


def locus_sort_key(lid):
    return (0, int(lid), "") if lid.isdigit() else (1, 0, lid)


def iupac(chars):
    """IUPAC code for a collection of nucleotide characters (missing data ignored)."""
    bases = set()
    for c in chars:
        bases.update(IUPAC_BASES.get(c, ""))
    return BASES_TO_IUPAC[frozenset(bases)]


def iupac_merge(seqs):
    """Column-wise IUPAC consensus of equal-length sequences."""
    if len(seqs) == 1:
        return seqs[0]
    return "".join(col[0] if len(set(col)) == 1 else iupac(col) for col in zip(*seqs))


def has_fixed_difference(codes):
    """True if any two IUPAC codes share no possible base (a fixed difference)."""
    sets = [IUPAC_SETS[c] for c in codes if IUPAC_SETS.get(c)]
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            if not sets[i] & sets[j]:
                return True
    return False


def chi2_sf(x, df):
    """Upper-tail probability of the chi-square distribution (replaces scipy)."""
    if df <= 0:
        return 1.0
    if x <= 0:
        return 1.0
    a, x = df / 2.0, x / 2.0
    log_prefix = -x + a * math.log(x) - math.lgamma(a)
    if x < a + 1:  # series expansion of the lower incomplete gamma
        term = total = 1.0 / a
        ap = a
        for _ in range(10000):
            ap += 1
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return max(0.0, 1.0 - total * math.exp(log_prefix))
    # continued fraction for the upper incomplete gamma (modified Lentz)
    tiny = 1e-300
    b = x + 1 - a
    c = 1 / tiny
    d = 1 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = d if abs(d) > tiny else tiny
        c = b + an / c
        c = c if abs(c) > tiny else tiny
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return math.exp(log_prefix) * h


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #

class Sample:
    __slots__ = ("name", "pop", "keys")

    def __init__(self, name, pop, keys):
        self.name = name  # name written to output files
        self.pop = pop
        self.keys = keys  # identifiers that may appear in the FASTA headers


class Locus:
    """One catalog locus.

    ``alleles`` holds the distinct haplotype sequences. ``g0``/``g1`` hold, for
    every sample index, the indices of the sample's two gene copies (equal for a
    homozygote, -1 for missing data).
    """

    __slots__ = ("id", "length", "alleles", "g0", "g1")

    def __init__(self, lid, length, alleles, g0, g1):
        self.id = lid
        self.length = length
        self.alleles = alleles
        self.g0 = g0
        self.g1 = g1

    def set_missing(self, s):
        self.g0[s] = -1
        self.g1[s] = -1

    def n_genotyped(self, samples):
        g0 = self.g0
        return sum(1 for s in samples if g0[s] >= 0)

    def allele_counts(self, samples):
        """Gene-copy counts per allele index among ``samples``."""
        counts = collections.Counter()
        g0, g1 = self.g0, self.g1
        for s in samples:
            a = g0[s]
            if a >= 0:
                counts[a] += 1
                counts[g1[s]] += 1
        return counts


class Dataset:
    def __init__(self, samples, loci, locus_order):
        self.samples = samples
        self.loci = loci
        self.locus_order = locus_order
        self.removed_samples = set()
        self.refresh()

    def refresh(self):
        """Recompute population membership, output order and locus order."""
        members = collections.OrderedDict()
        for i, smp in enumerate(self.samples):
            if i not in self.removed_samples:
                members.setdefault(smp.pop, []).append(i)
        self.members = members
        self.pops = list(members)
        self.order = [s for p in self.pops for s in members[p]]
        self.locus_order = [lid for lid in self.locus_order if lid in self.loci]

    def remove_samples(self, idxs):
        idxs = set(idxs)
        for loc in self.loci.values():
            for s in idxs:
                loc.set_missing(s)
        self.removed_samples |= idxs
        before = self.pops
        self.refresh()
        gone = [p for p in before if p not in self.members]
        if gone:
            info("  %d populations have no individuals left and were removed: %s"
                 % (len(gone), preview(gone)))

    def drop_loci(self, lids):
        for lid in lids:
            del self.loci[lid]
        self.refresh()

    def prune_empty_loci(self):
        empty = [lid for lid, loc in self.loci.items() if loc.n_genotyped(self.order) == 0]
        if empty:
            self.drop_loci(empty)
        return len(empty)

    def used_alleles(self, loc):
        used = set()
        g0, g1 = loc.g0, loc.g1
        for s in self.order:
            a = g0[s]
            if a >= 0:
                used.add(a)
                used.add(g1[s])
        return used


# --------------------------------------------------------------------------- #
# Input
# --------------------------------------------------------------------------- #

def read_popmap(path, fmt="auto"):
    """Read a population map.

    ``stacks`` format (Stacks 1.x/2.x ``-M`` popmap): ``sample<TAB>population``, no header.
    ``legacy`` format (fasta2genotype 1.x): a header row, then
    ``SampleID<TAB>IndividualID<TAB>PopulationID``.
    """
    rows = []
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            rows.append((lineno, [f.strip() for f in fields]))
    if not rows:
        raise InputError("population map '%s' is empty" % path)

    if fmt == "auto":
        fmt = "legacy" if all(len(r[1]) >= 3 for r in rows) else "stacks"
    samples = []
    if fmt == "legacy":
        header = rows.pop(0)[1]
        info("Population map: legacy 3-column format (header row '%s' skipped)"
             % "\t".join(header[:3]))
        for lineno, f in rows:
            if len(f) < 3 or not all(f[:3]):
                raise InputError("%s line %d: expected SampleID, IndividualID and "
                                 "PopulationID columns" % (path, lineno))
            samples.append(Sample(f[1], f[2], (f[0],)))
    else:
        info("Population map: Stacks format (sample, population)")
        for lineno, f in rows:
            if len(f) < 2 or not f[0] or not f[1]:
                raise InputError("%s line %d: expected sample and population columns"
                                 % (path, lineno))
            samples.append(Sample(f[0], f[1], (f[0],)))

    for ids in ([s.keys[0] for s in samples], [s.name for s in samples]):
        dups = [k for k, c in collections.Counter(ids).items() if c > 1]
        if dups:
            raise InputError("population map '%s' lists these samples more than once: %s"
                             % (path, preview(dups)))
    return samples


def read_locus_list(path):
    """Read a white/blacklist: one catalog locus ID per line (extra columns ignored)."""
    loci = []
    snp_columns = False
    with open_text(path) as fh:
        for line in fh:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            loci.append(fields[0])
            snp_columns = snp_columns or len(fields) > 1
    if snp_columns:
        warn("'%s' has more than one column; only the first (locus ID) is used and "
             "whole loci are kept or removed" % path)
    return list(dict.fromkeys(loci))


def iter_fasta(fh):
    header, parts = None, []
    for line in fh:
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts)
            header, parts = line.rstrip(), []
        elif line.startswith("#"):
            continue
        else:
            parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


def read_fasta(path, samples, whitelist=None, blacklist=None, clip5=(), clip3=(),
               keep_ambiguous=False):
    """Parse a Stacks haplotype FASTA into a Dataset."""
    n = len(samples)
    keymap = {}
    for i, smp in enumerate(samples):
        for key in smp.keys:
            keymap[key] = i
    white = set(whitelist) if whitelist else None
    black = set(blacklist) if blacklist else set()

    raw = {}  # locus -> [seq -> raw index, g0 array, g1 array]
    resolve_cache = {}
    unknown_samples = set()
    found = set()
    three_alleles = duplicated = records = 0
    stacks_version = None

    with open_text(path) as fh:
        first = fh.readline()
        m = re.match(r"#\s*Stacks version\s+(\S+?);?\s", first)
        if m:
            stacks_version = m.group(1)
        elif first.startswith(">"):
            fh = _prepend(first, fh)
        for header, seq in iter_fasta(fh):
            m = HEADER_RE.match(header)
            if not m:
                raise InputError("unrecognised FASTA header in '%s': %s\nExpected Stacks "
                                 "format '>CLocus_#_Sample_#_Locus_#_Allele_#'"
                                 % (path, header[:120]))
            records += 1
            lid, snum, allele, sname = m.groups()
            if (white is not None and lid not in white) or lid in black:
                continue
            ck = (snum, sname)
            s = resolve_cache.get(ck)
            if s is None:
                s = keymap.get(sname, keymap.get(snum, -1)) if sname else keymap.get(snum, -1)
                resolve_cache[ck] = s
            if s < 0:
                unknown_samples.add(sname or snum)
                continue
            found.add(s)
            if allele not in ("0", "1"):
                if not allele.isdigit():
                    raise InputError("non-numeric allele in FASTA header: %s" % header[:120])
                three_alleles += 1
                continue

            seq = seq.upper()
            for c in clip5:
                if seq.startswith(c):
                    seq = seq[len(c):]
                    break
            for c in clip3:
                if seq.endswith(c):
                    seq = seq[:len(seq) - len(c)]
                    break

            entry = raw.get(lid)
            if entry is None:
                entry = raw[lid] = [{}, array("i", [-1]) * n, array("i", [-1]) * n]
            index, g0, g1 = entry
            a = index.get(seq)
            if a is None:
                a = index[seq] = len(index)
            g = g0 if allele == "0" else g1
            if g[s] >= 0:
                duplicated += 1
            else:
                g[s] = a

    if records == 0:
        raise InputError("no sequences found in '%s'" % path)
    info("Read %d sequences from '%s' (%s)" % (
        records, path, "Stacks v" + stacks_version if stacks_version else "Stacks v1 format"))
    if unknown_samples:
        warn("%d FASTA sample(s) are not in the population map and were skipped: %s"
             % (len(unknown_samples), preview(sorted(unknown_samples))))
    if not found:
        raise InputError("none of the FASTA samples match the population map. Stacks 2 "
                         "FASTA files are matched by sample name, Stacks 1 files by the "
                         "number in 'Sample_#' (use a legacy 3-column population map).")
    absent = [samples[i].name for i in range(n) if i not in found]
    if absent:
        warn("%d population-map sample(s) have no sequences and were dropped: %s"
             % (len(absent), preview(absent)))
    if duplicated:
        info("  %d duplicate gene copies (same allele number twice) ignored" % duplicated)
    if three_alleles:
        warn("%d sequences with allele numbers above 1 ignored; only alleles 0 and 1 are "
             "used (diploid data assumed)" % three_alleles)

    loci = {}
    inconsistent = ambiguous = 0
    for lid, (index, g0, g1) in raw.items():
        rawseqs = list(index)
        length = len(rawseqs[0])
        if any(len(q) != length for q in rawseqs):
            inconsistent += 1
            continue
        if keep_ambiguous:
            alleles, remap = rawseqs, list(range(len(rawseqs)))
        else:
            alleles, remap = _resolve_ambiguity(rawseqs)
        for s in range(n):
            a, b = g0[s], g1[s]
            if a < 0 and b < 0:
                continue
            if a < 0:
                a = b
            elif b < 0:
                b = a
            a, b = remap[a], remap[b]
            if a < 0 or b < 0:
                ambiguous += 1
                a = b = -1
            g0[s], g1[s] = a, b
        loci[lid] = Locus(lid, length, alleles, g0, g1)
    if inconsistent:
        warn("%d loci with sequences of differing lengths were removed (check the "
             "--clip options)" % inconsistent)
    if ambiguous:
        info("  %d genotypes with uncalled (N) bases at SNP positions set to missing "
             "(use --keep-ambiguous to keep them)" % ambiguous)

    if whitelist:
        order = [lid for lid in whitelist if lid in loci]
        missing = len(whitelist) - len(order)
        if missing:
            info("  %d whitelisted loci were not found in the FASTA file" % missing)
    else:
        order = sorted(loci, key=locus_sort_key)
    return Dataset(samples, loci, order)


def _prepend(first, fh):
    yield first
    for line in fh:
        yield line


def _resolve_ambiguity(rawseqs):
    """Canonicalise haplotypes containing missing bases (N).

    Returns (alleles, remap) where remap[raw_index] is the canonical allele index,
    or -1 when the haplotype has a missing base at a variable (SNP) column and its
    identity is therefore unknown. Columns that are invariant apart from missing
    bases are set to N in every haplotype so they cannot split identical alleles.
    """
    if all(set(q) <= ACGT for q in rawseqs):
        return rawseqs, list(range(len(rawseqs)))
    variable, mask = [], []
    for i, col in enumerate(zip(*rawseqs)):
        chars = set(col)
        if chars <= ACGT:
            if len(chars) > 1:
                variable.append(i)
        elif len(chars & ACGT) >= 2:
            variable.append(i)
        elif chars & ACGT:
            mask.append(i)
    alleles, index, remap = [], {}, []
    for q in rawseqs:
        if any(q[i] not in ACGT for i in variable):
            remap.append(-1)
            continue
        if mask:
            q = list(q)
            for i in mask:
                q[i] = "N"
            q = "".join(q)
        a = index.get(q)
        if a is None:
            a = index[q] = len(alleles)
            alleles.append(q)
        remap.append(a)
    return alleles, remap


def read_vcf_depths(path, ds, stat="mean"):
    """Per-sample read depth for each locus, aggregated over the locus' SNP records.

    Returns {locus_id: [depth per sample index]} (-1 where unknown).
    """
    n = len(ds.samples)
    by_name = {}
    for i, smp in enumerate(ds.samples):
        for key in smp.keys + (smp.name,):
            by_name.setdefault(key, i)
    cols = None
    acc = {}
    dp_index_cache = {}
    records = 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                names = line.rstrip("\r\n").split("\t")[9:]
                cols = [(j, by_name[nm]) for j, nm in enumerate(names) if nm in by_name]
                matched = {i for _, i in cols}
                unmatched = [ds.samples[i].name for i in ds.order if i not in matched]
                if not cols:
                    raise InputError("no VCF sample names match the population map")
                if unmatched:
                    warn("%d sample(s) are not in the VCF; the coverage filter will set all "
                         "their genotypes to missing: %s" % (len(unmatched), preview(unmatched)))
                continue
            if cols is None:
                raise InputError("'%s' has no #CHROM header line; is it a VCF file?" % path)
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                continue
            lid = fields[2].split(":")[0]
            if lid not in ds.loci:
                continue
            fmt = fields[8]
            dpi = dp_index_cache.get(fmt)
            if dpi is None:
                keys = fmt.split(":")
                if "DP" not in keys:
                    raise InputError("VCF record for locus %s has no DP (read depth) field. "
                                     "Use populations.snps.vcf rather than "
                                     "populations.haps.vcf." % lid)
                dpi = dp_index_cache[fmt] = keys.index("DP")
            records += 1
            a = acc.get(lid)
            if a is None:
                a = acc[lid] = ([0] * n, [0] * n, [None] * n, [None] * n)
            tot, cnt, lo, hi = a
            for j, s in cols:
                parts = fields[9 + j].split(":")
                v = parts[dpi] if len(parts) > dpi else "."
                dp = int(v) if v.isdigit() else 0
                tot[s] += dp
                cnt[s] += 1
                lo[s] = dp if lo[s] is None or dp < lo[s] else lo[s]
                hi[s] = dp if hi[s] is None or dp > hi[s] else hi[s]
    info("Read %d VCF records from '%s' covering %d loci" % (records, path, len(acc)))
    depths = {}
    for lid, (tot, cnt, lo, hi) in acc.items():
        if stat == "mean":
            vals = [t / c if c else -1 for t, c in zip(tot, cnt)]
        elif stat == "min":
            vals = [-1 if v is None else v for v in lo]
        else:
            vals = [-1 if v is None else v for v in hi]
        depths[lid] = vals
    return depths


# --------------------------------------------------------------------------- #
# Filters
# --------------------------------------------------------------------------- #

def filter_coverage(ds, depths, min_cov):
    info("Removing genotypes with read depth below %s..." % min_cov)
    no_vcf = [lid for lid in ds.loci if lid not in depths]
    if no_vcf:
        info("  %d loci have no VCF records (e.g. monomorphic loci) and were removed"
             % len(no_vcf))
        ds.drop_loci(no_vcf)
    removed = 0
    for lid, loc in ds.loci.items():
        d = depths[lid]
        for s in ds.order:
            if loc.g0[s] >= 0 and d[s] < min_cov:
                loc.set_missing(s)
                removed += 1
    info("  %d genotypes set to missing" % removed)


def filter_monomorphic(ds):
    mono = [lid for lid, loc in ds.loci.items() if len(loc.allele_counts(ds.order)) < 2]
    if mono:
        ds.drop_loci(mono)
    return len(mono)


def filter_heterozygosity(ds, cutoff, alpha):
    """Remove loci whose observed heterozygosity is >= cutoff, exceeds HWE expectation,
    and departs significantly from Hardy-Weinberg proportions (likely paralogs)."""
    info("Removing loci with excess heterozygosity (Ho >= %s, HWE p < %s)..." % (cutoff, alpha))
    flagged = []
    for lid, loc in ds.loci.items():
        g0, g1 = loc.g0, loc.g1
        counts = collections.Counter()
        genos = collections.Counter()
        n = het = 0
        for s in ds.order:
            a = g0[s]
            if a < 0:
                continue
            b = g1[s]
            n += 1
            counts[a] += 1
            counts[b] += 1
            if a != b:
                het += 1
            genos[(a, b) if a <= b else (b, a)] += 1
        if n == 0 or len(counts) < 2 or het / n < cutoff:
            continue
        freqs = {a: c / (2.0 * n) for a, c in counts.items()}
        if het <= n * (1 - sum(p * p for p in freqs.values())):
            continue
        alleles = sorted(freqs)
        chi = 0.0
        for i, a in enumerate(alleles):
            for b in alleles[i:]:
                exp = n * freqs[a] ** 2 if a == b else 2 * n * freqs[a] * freqs[b]
                chi += (genos.get((a, b), 0) - exp) ** 2 / exp
        k = len(alleles)
        if chi2_sf(chi, k * (k - 1) // 2) < alpha:
            flagged.append(lid)
    ds.drop_loci(flagged)
    info("  %d loci removed" % len(flagged))


def filter_alleles(ds, min_freq, min_pops):
    """Set genotypes carrying rare alleles to missing."""
    if min_freq:
        info("Removing alleles with overall frequency below %s..." % min_freq)
    if min_pops:
        info("Removing alleles found in fewer than %s of populations..." % min_pops)
    n_pops = len(ds.pops)
    removed_alleles = removed_genos = 0
    for loc in ds.loci.values():
        counts = loc.allele_counts(ds.order)
        total = sum(counts.values())
        if not total:
            continue
        flagged = set()
        if min_freq:
            flagged.update(a for a, c in counts.items() if c / total < min_freq)
        if min_pops:
            pops_with = collections.Counter()
            for pop in ds.pops:
                pops_with.update(loc.allele_counts(ds.members[pop]).keys())
            flagged.update(a for a in counts if pops_with[a] / n_pops < min_pops)
        if flagged:
            removed_alleles += len(flagged)
            g0, g1 = loc.g0, loc.g1
            for s in ds.order:
                if g0[s] in flagged or g1[s] in flagged:
                    loc.set_missing(s)
                    removed_genos += 1
    info("  %d alleles removed (%d genotypes set to missing)" % (removed_alleles, removed_genos))


def filter_missing(ds, pop_thresh, locus_thresh, ind_thresh):
    if pop_thresh:
        info("Removing loci from populations where fewer than %s of individuals are "
             "genotyped..." % pop_thresh)
        removed = 0
        for loc in ds.loci.values():
            for pop in ds.pops:
                members = ds.members[pop]
                g = loc.n_genotyped(members)
                if g and g / len(members) < pop_thresh:
                    for s in members:
                        loc.set_missing(s)
                    removed += 1
        info("  %d locus x population combinations removed" % removed)
        ds.prune_empty_loci()

    if locus_thresh:
        info("Removing loci genotyped in fewer than %s of individuals..." % locus_thresh)
        n = len(ds.order)
        low = [lid for lid, loc in ds.loci.items() if loc.n_genotyped(ds.order) / n < locus_thresh]
        ds.drop_loci(low)
        info("  %d loci removed" % len(low))

    if ind_thresh:
        info("Removing individuals genotyped at fewer than %s of loci..." % ind_thresh)
        ds.prune_empty_loci()
        if not ds.loci:
            raise InputError("all loci were removed by the filters; check the thresholds")
        have = collections.Counter()
        for loc in ds.loci.values():
            g0 = loc.g0
            for s in ds.order:
                if g0[s] >= 0:
                    have[s] += 1
        nloci = len(ds.loci)
        low = [s for s in ds.order if have[s] / nloci < ind_thresh]
        if low:
            ds.remove_samples(low)
        info("  %d individuals removed%s" % (
            len(low), (": " + preview(ds.samples[s].name for s in low)) if low else ""))


def run_filters(ds, args, depths=None):
    if depths is not None:
        filter_coverage(ds, depths, args.min_coverage)
    if args.remove_monomorphic:
        info("Removing monomorphic loci...")
        info("  %d loci removed" % filter_monomorphic(ds))
    if args.het_cutoff is not None:
        filter_heterozygosity(ds, args.het_cutoff, args.hwe_alpha)
    if args.min_allele_freq or args.min_allele_pops:
        filter_alleles(ds, args.min_allele_freq, args.min_allele_pops)
    if args.min_pop_locus_freq or args.min_locus_freq or args.min_ind_freq:
        filter_missing(ds, args.min_pop_locus_freq, args.min_locus_freq, args.min_ind_freq)
    empty = ds.prune_empty_loci()
    if empty:
        info("Removed %d loci left with no genotypes" % empty)
    if args.remove_monomorphic:
        mono = filter_monomorphic(ds)
        if mono:
            info("Removed %d loci that became monomorphic after filtering" % mono)
    have = set()
    for loc in ds.loci.values():
        have.update(s for s in ds.order if loc.g0[s] >= 0)
    empty = [s for s in ds.order if s not in have]
    if empty:
        info("Removed %d individuals left with no genotypes: %s"
             % (len(empty), preview(ds.samples[s].name for s in empty)))
        ds.remove_samples(empty)
    if not ds.loci:
        raise InputError("no loci remain after filtering; check the input files and thresholds")
    if not ds.order:
        raise InputError("no individuals remain after filtering")


# --------------------------------------------------------------------------- #
# Shared helpers for writers
# --------------------------------------------------------------------------- #

class OutputNamer:
    LEGACY_SUFFIX = {"pops": "_pops.out", "loci": "_loci.out", "snps": ".snp",
                     "inds": "_inds.out", "names": "_names.out"}

    def __init__(self, prefix, legacy=False):
        self.prefix = prefix
        self.legacy = legacy

    def main(self, ext):
        return self.prefix + (".out" if self.legacy else ext)

    def aux(self, ext, kind):
        if self.legacy:
            return self.prefix + self.LEGACY_SUFFIX[kind]
        base = re.sub(r"(\.gz|\.txt|\.tsv)+$", "", ext)
        return self.prefix + base + "." + kind + ".tsv"


def write_table(path, header, rows):
    with open_out(path) as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")
    return path


def short_names(ds, maxlen, namer, ext, program):
    """Sample names, replaced by short codes if any exceeds ``maxlen`` characters."""
    names = {s: ds.samples[s].name for s in ds.order}
    if max(len(v) for v in names.values()) <= maxlen:
        return names, []
    width = len(str(len(ds.order)))
    short = {s: "I%0*d" % (width, i + 1) for i, s in enumerate(ds.order)}
    path = write_table(namer.aux(ext, "names"), ["code", "sample", "population"],
                       [(short[s], names[s], ds.samples[s].pop) for s in ds.order])
    warn("%s allows at most %d-character names; samples were renamed %s, %s, ... "
         "(see %s)" % (program, maxlen, short[ds.order[0]], "I%0*d" % (width, 2), path))
    return short, [path]


def haplotype_codes(ds):
    """Integer code (1, 2, ...) per allele, numbered by first appearance in output order."""
    codes = {}
    for lid in ds.locus_order:
        loc = ds.loci[lid]
        g0, g1 = loc.g0, loc.g1
        m = {}
        for s in ds.order:
            a = g0[s]
            if a >= 0:
                if a not in m:
                    m[a] = len(m) + 1
                b = g1[s]
                if b not in m:
                    m[b] = len(m) + 1
        codes[lid] = m
    return codes


def variable_columns(loc, used):
    seqs = [loc.alleles[a] for a in sorted(used)]
    cols = []
    if len(seqs) < 2:
        return cols
    for i, col in enumerate(zip(*seqs)):
        if len(set(col) & ACGT) >= 2:
            cols.append(i)
    return cols


def biallelic_snps(ds, one_snp=False):
    """Yield (locus, column, base_by_allele, (base1, base2)) for biallelic SNPs.

    base1 < base2 alphabetically; base_by_allele maps allele index -> base (None if
    not an unambiguous A/C/G/T).
    """
    for lid in ds.locus_order:
        loc = ds.loci[lid]
        used = ds.used_alleles(loc)
        for col in variable_columns(loc, used):
            bases = {a: loc.alleles[a][col] for a in used}
            observed = set(bases.values()) & ACGT
            if len(observed) != 2:
                continue
            by_allele = {a: (b if b in ACGT else None) for a, b in bases.items()}
            yield loc, col, by_allele, tuple(sorted(observed))
            if one_snp:
                break


def snp_genotype_counts(ds, loc, by_allele, alt):
    """Per-sample count (0/1/2) of base ``alt`` in output order; None if missing."""
    out = []
    g0, g1 = loc.g0, loc.g1
    for s in ds.order:
        a = g0[s]
        if a < 0:
            out.append(None)
            continue
        x, y = by_allele[a], by_allele[g1[s]]
        out.append(None if x is None or y is None else (x == alt) + (y == alt))
    return out


def minor_base(genos_by_base):
    """Pick the minor allele of a biallelic SNP from {base: gene-copy count}."""
    (b1, c1), (b2, c2) = sorted(genos_by_base.items())
    return b1 if c1 < c2 else b2


def title_of(opts):
    return opts.title or os.path.basename(opts.out)


# --------------------------------------------------------------------------- #
# Writers: sequence formats
# --------------------------------------------------------------------------- #

def write_migrate(ds, namer, opts):
    ext = ".migrate"
    path = namer.main(ext)
    names, extra = short_names(ds, 9, namer, ext, "migrate-n")
    loci = [ds.loci[lid] for lid in ds.locus_order]
    with open_out(path) as f:
        first = "%d\t%d" % (len(ds.pops), len(loci))
        f.write(first + ("\t" + opts.title if opts.title else "") + "\n")
        f.write("\t".join(str(loc.length) for loc in loci) + "\n")
        for pop in ds.pops:
            members = ds.members[pop]
            f.write("%d\tPop_%s\n" % (2 * len(members), pop))
            for loc in loci:
                missing = "?" * loc.length
                alleles, g0, g1 = loc.alleles, loc.g0, loc.g1
                for s in members:
                    a = g0[s]
                    sa, sb = (missing, missing) if a < 0 else (alleles[a], alleles[g1[s]])
                    nm = names[s]
                    f.write("%-10s%s\n%-10s%s\n" % (nm + "a", sa, nm + "b", sb))
    return [path] + extra


def _arlequin(ds, path, opts, datatype, cell):
    """cell(loc, allele_index or None) -> string for one gene copy."""
    loci = [ds.loci[lid] for lid in ds.locus_order]
    with open_out(path) as f:
        f.write("[Profile]\n\n")
        f.write('\tTitle="%s"\n' % title_of(opts))
        f.write("\tNbSamples=%d\n" % len(ds.pops))
        f.write("\tGenotypicData=1\n\tGameticPhase=0\n")
        f.write("\tDataType=%s\n\tLocusSeparator=TAB\n\tMissingData='?'\n\n" % datatype)
        f.write("[Data]\n\n\t[[Samples]]\n\n")
        for pop in ds.pops:
            members = ds.members[pop]
            f.write('\t\tSampleName="%s"\n\t\tSampleSize=%d\n\t\tSampleData={\n'
                    % (pop, len(members)))
            for s in members:
                row_a, row_b = [], []
                for loc in loci:
                    a = loc.g0[s]
                    if a < 0:
                        row_a.append(cell(loc, None))
                        row_b.append(cell(loc, None))
                    else:
                        row_a.append(cell(loc, a))
                        row_b.append(cell(loc, loc.g1[s]))
                f.write("%s\t1\t%s\n" % (ds.samples[s].name, "\t".join(row_a)))
                f.write("\t\t%s\n" % "\t".join(row_b))
            f.write("}\n\n")
    return [path]


def write_arlequin(ds, namer, opts):
    def cell(loc, a):
        return "?" * loc.length if a is None else loc.alleles[a]
    return _arlequin(ds, namer.main(".arp"), opts, "DNA", cell)


def write_diyabc(ds, namer, opts):
    path = namer.main(".diyabc")
    loci = [ds.loci[lid] for lid in ds.locus_order]
    with open_out(path) as f:
        f.write("%s <NM=1.0NF>\n" % title_of(opts))
        for loc in loci:
            f.write("%s\t<A>\n" % loc.id)
        for pop in ds.pops:
            f.write("Pop\n")
            for s in ds.members[pop]:
                cells = []
                for loc in loci:
                    a = loc.g0[s]
                    if a < 0:
                        cells.append("<[][]>")
                    else:
                        cells.append("<[%s][%s]>" % (loc.alleles[a], loc.alleles[loc.g1[s]]))
                f.write("%s\t,\t%s\n" % (ds.samples[s].name, "\t".join(cells)))
    return [path]


def write_gphocs(ds, namer, opts):
    path = namer.main(".gphocs")
    with open_out(path) as f:
        f.write("%d\n\n" % len(ds.locus_order))
        for lid in ds.locus_order:
            loc = ds.loci[lid]
            alleles, g0, g1 = loc.alleles, loc.g0, loc.g1
            cache = {}
            rows = []
            for s in ds.order:
                a = g0[s]
                if a < 0:
                    continue
                b = g1[s]
                key = (a, b)
                seq = cache.get(key)
                if seq is None:
                    seq = cache[key] = iupac_merge([alleles[a], alleles[b]])
                rows.append("%s\t%s\n" % (ds.samples[s].name, seq))
            f.write("%s\t%d\t%d\n" % (lid, len(rows), loc.length))
            f.writelines(rows)
            f.write("\n")
    return [path]


def write_phylip(ds, namer, opts):
    path = namer.main(".phy")
    sites, mode, keep = opts.phylip_sites, opts.phylip_mode, opts.phylip_loci

    if mode == "haploid":
        units = [(ds.samples[s].name + c, s, i) for s in ds.order for i, c in enumerate("ab")]
    elif mode == "diploid":
        units = [(ds.samples[s].name, s, None) for s in ds.order]
    else:
        units = [(pop, ds.members[pop], None) for pop in ds.pops]

    blocks = []  # (locus id, width, [sequence or None per unit])
    for lid in ds.locus_order:
        loc = ds.loci[lid]
        used = ds.used_alleles(loc)
        if sites == "snps":
            cols = variable_columns(loc, used)
            proj = {a: "".join(loc.alleles[a][c] for c in cols) for a in used}
            width = len(cols)
        else:
            proj = {a: loc.alleles[a] for a in used}
            width = loc.length
        g0, g1 = loc.g0, loc.g1
        cache = {}
        seqs = []
        for _, who, copy in units:
            # key = the allele indices summarised by this row at this locus
            if mode == "population":
                key = frozenset(x for s in who if g0[s] >= 0 for x in (g0[s], g1[s]))
            elif g0[who] < 0:
                key = None
            elif mode == "haploid":
                key = (g0[who] if copy == 0 else g1[who],)
            else:
                key = (g0[who], g1[who])
            if not key:
                seqs.append(None)
                continue
            seq = cache.get(key)
            if seq is None:
                seq = cache[key] = iupac_merge([proj[x] for x in sorted(set(key))])
            seqs.append(seq)

        if keep != "all":
            tally = collections.Counter(q for q in seqs if q is not None)
            kept = []
            for c in range(width):
                states = collections.Counter()
                for q, mult in tally.items():
                    if IUPAC_SETS.get(q[c]):
                        states[q[c]] += mult
                if keep == "pi":
                    states = [x for x, m in states.items() if m >= 2]
                if has_fixed_difference(states):
                    kept.append(c)
            if not kept:
                continue
            if sites == "snps":
                seqs = [None if q is None else "".join(q[c] for c in kept) for q in seqs]
                width = len(kept)
        if width:
            blocks.append((lid, width, seqs))

    if not blocks:
        raise InputError("no loci qualify for the PHYLIP alignment with these options")
    if keep != "all":
        info("  %d of %d loci kept as %s" % (
            len(blocks), len(ds.locus_order),
            "phylogenetically informative" if keep == "pi" else "having fixed differences"))
    namewidth = max(10, max(len(u[0]) for u in units) + 1)
    bang = "!" if opts.phylip_breakpoints else ""
    with open_out(path) as f:
        f.write("%d %d\n" % (len(units), sum(b[1] for b in blocks)))
        if opts.phylip_locus_header:
            f.write("\t" + "\t".join(b[0] for b in blocks) + "\n")
        for i, (name, _, _) in enumerate(units):
            parts = [name.ljust(namewidth)]
            for _, width, seqs in blocks:
                parts.append(seqs[i] if seqs[i] is not None else "N" * width)
                parts.append(bang)
            f.write("".join(parts) + "\n")
    return [path]


# --------------------------------------------------------------------------- #
# Writers: SNP formats
# --------------------------------------------------------------------------- #

def _snp_matrix(ds, opts):
    """Biallelic SNPs coded as minor-allele counts: [(snp_info, [count or None per sample])]."""
    out = []
    for loc, col, by_allele, (b1, b2) in biallelic_snps(ds, opts.one_snp):
        counts = snp_genotype_counts(ds, loc, by_allele, b2)
        n2 = sum(c for c in counts if c is not None)
        n1 = 2 * sum(1 for c in counts if c is not None) - n2
        minor = minor_base({b1: n1, b2: n2})
        major = b1 if minor == b2 else b2
        if minor == b1:
            counts = [None if c is None else 2 - c for c in counts]
        out.append(((loc.id, col, major, minor), counts))
    if not out:
        raise InputError("no biallelic SNPs remain after filtering")
    info("  %d biallelic SNPs%s" % (len(out), " (one per locus)" if opts.one_snp else ""))
    return out


def write_lfmm(ds, namer, opts):
    ext = ".lfmm"
    path = namer.main(ext)
    snps = _snp_matrix(ds, opts)
    with open_out(path) as f:
        for i in range(len(ds.order)):
            f.write(" ".join("9" if c[i] is None else str(c[i]) for _, c in snps) + "\n")
    p1 = write_table(namer.aux(ext, "snps"), ["snp", "locus", "column", "major", "minor"],
                     [("%s_%d" % (lid, col), lid, col, maj, mnr)
                      for (lid, col, maj, mnr), _ in snps])
    p2 = write_table(namer.aux(ext, "inds"), ["individual", "population"],
                     [(ds.samples[s].name, ds.samples[s].pop) for s in ds.order])
    return [path, p1, p2]


def write_diyabc_snp(ds, namer, opts):
    path = namer.main(".diyabc.snp")
    snps = _snp_matrix(ds, opts)
    with open_out(path) as f:
        f.write("%s <NM=1.0NF> <MAF=hudson>\n" % title_of(opts))
        f.write("IND\tSEX\tPOP\t" + "\t".join("A" for _ in snps) + "\n")
        for i, s in enumerate(ds.order):
            smp = ds.samples[s]
            f.write("%s\t9\t%s\t%s\n" % (smp.name, smp.pop, "\t".join(
                "9" if c[i] is None else str(c[i]) for _, c in snps)))
    return [path]


def write_treemix(ds, namer, opts):
    ext = ".treemix.gz"
    path = namer.main(ext)
    n = 0
    with open_out(path) as f:
        f.write(" ".join(ds.pops) + "\n")
        for loc, col, by_allele, (b1, b2) in biallelic_snps(ds, opts.one_snp):
            g0, g1 = loc.g0, loc.g1
            cells = []
            for pop in ds.pops:
                c1 = c2 = 0
                for s in ds.members[pop]:
                    a = g0[s]
                    if a < 0:
                        continue
                    x, y = by_allele[a], by_allele[g1[s]]
                    if x is None or y is None:
                        continue
                    c1 += (x == b1) + (y == b1)
                    c2 += (x == b2) + (y == b2)
                cells.append("%d,%d" % (c1, c2))
            f.write(" ".join(cells) + "\n")
            n += 1
    if not n:
        raise InputError("no biallelic SNPs remain after filtering")
    info("  %d biallelic SNPs%s" % (n, " (one per locus)" if opts.one_snp else ""))
    return [path]


# --------------------------------------------------------------------------- #
# Writers: haplotype (allele-integer) formats
# --------------------------------------------------------------------------- #

def _pop_table(ds, namer, ext):
    return write_table(namer.aux(ext, "pops"), ["number", "population"],
                       [(i + 1, p) for i, p in enumerate(ds.pops)])


def _loci_table(ds, namer, ext):
    return write_table(namer.aux(ext, "loci"), ["number", "locus"],
                       [(i + 1, lid) for i, lid in enumerate(ds.locus_order)])


def _coded(ds, codes, s):
    """[(code_a, code_b) or None] for sample s across loci in output order."""
    out = []
    for lid in ds.locus_order:
        loc = ds.loci[lid]
        a = loc.g0[s]
        out.append(None if a < 0 else (codes[lid][a], codes[lid][loc.g1[s]]))
    return out


def write_structure(ds, namer, opts):
    ext = ".str"
    path = namer.main(ext)
    codes = haplotype_codes(ds)
    popnum = {p: i + 1 for i, p in enumerate(ds.pops)}
    with open_out(path) as f:
        f.write("\t\t" + "\t".join(ds.locus_order) + "\n")
        for s in ds.order:
            smp = ds.samples[s]
            gts = _coded(ds, codes, s)
            for copy in (0, 1):
                f.write("%s\t%d\t%s\n" % (smp.name, popnum[smp.pop], "\t".join(
                    "-9" if g is None else str(g[copy]) for g in gts)))
    return [path, _pop_table(ds, namer, ext)]


def write_genepop(ds, namer, opts):
    ext = ".gen"
    path = namer.main(ext)
    codes = haplotype_codes(ds)
    width = 2 if opts.genepop_digits == 4 else 3
    most = max((len(m) for m in codes.values()), default=0)
    if most >= 10 ** width:
        raise InputError("a locus has %d alleles, too many for %d-digit Genepop format; "
                         "use --genepop-digits 6" % (most, opts.genepop_digits))
    missing = "0" * (2 * width)
    with open_out(path) as f:
        f.write(title_of(opts) + "\n")
        for lid in ds.locus_order:
            f.write(lid + "\n")
        for pop in ds.pops:
            f.write("Pop\n")
            for s in ds.members[pop]:
                cells = [missing if g is None else "%0*d%0*d" % (width, g[0], width, g[1])
                         for g in _coded(ds, codes, s)]
                f.write("%s ,  %s\n" % (ds.samples[s].name, " ".join(cells)))
    return [path, _pop_table(ds, namer, ext)]


def _pop_allele_counts(ds, codes, pop):
    """{locus: Counter(code -> gene copies)} for one population."""
    out = {}
    for lid in ds.locus_order:
        loc = ds.loci[lid]
        m = codes[lid]
        out[lid] = collections.Counter({m[a]: c for a, c in loc.allele_counts(ds.members[pop]).items()})
    return out


def write_allele_freq(ds, namer, opts):
    path = namer.main(".allelefreq.tsv")
    codes = haplotype_codes(ds)
    header = [""] + ["%s_%d" % (lid, c) for lid in ds.locus_order
                     for c in range(1, len(codes[lid]) + 1)]
    with open_out(path) as f:
        f.write("\t".join(header) + "\n")
        for pop in ds.pops:
            counts = _pop_allele_counts(ds, codes, pop)
            row = [pop]
            for lid in ds.locus_order:
                total = sum(counts[lid].values())
                for c in range(1, len(codes[lid]) + 1):
                    row.append("%.5f" % (counts[lid][c] / total if total else 0.0))
            f.write("\t".join(row) + "\n")
    return [path]


def write_sambada(ds, namer, opts):
    path = namer.main(".sambada.txt")
    codes = haplotype_codes(ds)
    nalleles = [len(codes[lid]) for lid in ds.locus_order]
    with open_out(path) as f:
        f.write("\t".join(["ID"] + ["%s_%d" % (lid, c) for lid, k in zip(ds.locus_order, nalleles)
                                    for c in range(1, k + 1)]) + "\n")
        for s in ds.order:
            gts = _coded(ds, codes, s)
            for copy, suffix in ((0, "a"), (1, "b")):
                row = [ds.samples[s].name + suffix]
                for g, k in zip(gts, nalleles):
                    if g is None:
                        row.extend(["NaN"] * k)
                    else:
                        row.extend("1" if c == g[copy] else "0" for c in range(1, k + 1))
                f.write("\t".join(row) + "\n")
    return [path]


def write_bayescan(ds, namer, opts):
    ext = ".bayescan.txt"
    path = namer.main(ext)
    codes = haplotype_codes(ds)
    with open_out(path) as f:
        f.write("[loci]=%d\n\n[populations]=%d\n" % (len(ds.locus_order), len(ds.pops)))
        for i, pop in enumerate(ds.pops):
            f.write("\n[pop]=%d\n" % (i + 1))
            counts = _pop_allele_counts(ds, codes, pop)
            for j, lid in enumerate(ds.locus_order):
                k = len(codes[lid])
                c = counts[lid]
                f.write("%d\t%d\t%d\t%s\n" % (j + 1, sum(c.values()), k,
                                              "\t".join(str(c[x]) for x in range(1, k + 1))))
    return [path, _loci_table(ds, namer, ext), _pop_table(ds, namer, ext)]


def write_arlequin_hap(ds, namer, opts):
    codes = haplotype_codes(ds)

    def cell(loc, a):
        return "?" if a is None else str(codes[loc.id][a])
    return _arlequin(ds, namer.main(".hap.arp"), opts, "STANDARD", cell)


def write_genalex(ds, namer, opts):
    path = namer.main(".genalex.txt")
    codes = haplotype_codes(ds)
    with open_out(path) as f:
        f.write("\t".join([str(len(ds.locus_order)), str(len(ds.order)), str(len(ds.pops))]
                          + [str(len(ds.members[p])) for p in ds.pops]) + "\n")
        f.write("\t".join([title_of(opts), "", ""] + ds.pops) + "\n")
        f.write("IndID\tPopID\t" + "\t\t".join(ds.locus_order) + "\n")
        for s in ds.order:
            smp = ds.samples[s]
            cells = []
            for g in _coded(ds, codes, s):
                cells.extend(("0", "0") if g is None else (str(g[0]), str(g[1])))
            f.write("%s\t%s\t%s\n" % (smp.name, smp.pop, "\t".join(cells)))
    return [path]


FORMATS = collections.OrderedDict([
    ("migrate", (write_migrate, "migrate-n DNA sequence file")),
    ("arlequin", (write_arlequin, "Arlequin project, DNA sequences (.arp)")),
    ("diyabc", (write_diyabc, "DIYABC sequence (genepop-like) file")),
    ("diyabc-snp", (write_diyabc_snp, "DIYABC-RF SNP file, one row per individual")),
    ("lfmm", (write_lfmm, "LFMM/LEA genotype matrix (0/1/2, 9 = missing)")),
    ("phylip", (write_phylip, "PHYLIP alignment (see --phylip-* options)")),
    ("gphocs", (write_gphocs, "G-PhoCS sequence file")),
    ("treemix", (write_treemix, "TreeMix allele counts (gzipped)")),
    ("structure", (write_structure, "Structure, haplotypes coded as integers")),
    ("genepop", (write_genepop, "Genepop, haplotypes coded as integers")),
    ("allele-freq", (write_allele_freq, "haplotype frequencies per population")),
    ("sambada", (write_sambada, "samBada presence/absence per gene copy")),
    ("bayescan", (write_bayescan, "BayeScan haplotype counts")),
    ("arlequin-hap", (write_arlequin_hap, "Arlequin project, haplotypes coded as integers")),
    ("genalex", (write_genalex, "GenAlEx, haplotypes coded as integers")),
])


# --------------------------------------------------------------------------- #
# Command line
# --------------------------------------------------------------------------- #

def fraction(text):
    try:
        v = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError("'%s' is not a number" % text)
    if not 0 <= v <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return v


def dna(text):
    if not re.match(r"^[ACGTacgt]+$", text):
        raise argparse.ArgumentTypeError("'%s' is not a DNA sequence (A/C/G/T)" % text)
    return text.upper()


def nonneg_int(text):
    try:
        v = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("'%s' is not an integer" % text)
    if v < 0:
        raise argparse.ArgumentTypeError("must be 0 or greater")
    return v


EXAMPLES = """\
examples:
  # Stacks 2: Structure and Genepop files, basic filtering
  %(prog)s -f populations.samples.fa -p popmap.tsv -o results/toads \\
      -F structure genepop --remove-monomorphic --min-locus-freq 0.8

  # PHYLIP alignment of parsimony-informative SNPs, one row per population
  %(prog)s -f populations.samples.fa -p popmap.tsv -o toads -F phylip \\
      --phylip-sites snps --phylip-mode population --phylip-loci pi

  # Stacks 1 input with a read-depth filter, clipping the SbfI cut site
  %(prog)s -f batch_1.fa -p popmap_legacy.tsv --vcf batch_1.vcf --min-coverage 10 \\
      --clip-5prime TGCAGG -o toads -F migrate

  # Old-style interactive mode (fasta2genotype 1.x arguments; use NA to omit a file)
  %(prog)s batch_1.fa whitelist.txt popmap_legacy.tsv NA toads

Please cite:
""" + CITATION


def build_parser():
    p = argparse.ArgumentParser(
        description="Convert a Stacks haplotype FASTA file (populations.samples.fa) into "
                    "input files for population-genetic software.",
        epilog=EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("legacy", nargs="*", help=argparse.SUPPRESS)
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    p.add_argument("-q", "--quiet", action="store_true", help="only print warnings and errors")

    g = p.add_argument_group("input")
    g.add_argument("-f", "--fasta", metavar="FILE",
                   help="Stacks haplotype FASTA: populations.samples.fa (Stacks 2) or "
                        "batch_N.fa (Stacks 1); may be gzipped")
    g.add_argument("-p", "--popmap", metavar="FILE",
                   help="population map: 'sample<TAB>pop' (Stacks) or the 3-column "
                        "fasta2genotype 1.x format with a header row")
    g.add_argument("--popmap-format", choices=["auto", "stacks", "legacy"], default="auto",
                   help="population map format (default: auto = legacy if every row has "
                        "3+ columns)")
    g.add_argument("-w", "--whitelist", metavar="FILE",
                   help="keep only these catalog loci (one ID per line); also sets the "
                        "output order of loci")
    g.add_argument("-b", "--blacklist", metavar="FILE", help="remove these catalog loci")
    g.add_argument("--vcf", metavar="FILE",
                   help="VCF from the same populations run, for --min-coverage "
                        "(populations.snps.vcf in Stacks 2)")

    g = p.add_argument_group("output")
    g.add_argument("-o", "--out", metavar="PREFIX", help="output file prefix (may include a directory)")
    g.add_argument("-F", "--format", nargs="+", choices=list(FORMATS) + ["all"], metavar="FORMAT",
                   help="one or more output formats, or 'all': " + ", ".join(FORMATS))
    g.add_argument("--title", help="project title for Arlequin, DIYABC, Genepop, GenAlEx and "
                                   "migrate-n (default: output prefix)")
    g.add_argument("--one-snp", action="store_true",
                   help="treemix/lfmm/diyabc-snp: use only the first biallelic SNP per locus")
    g.add_argument("--genepop-digits", type=int, choices=[4, 6], default=6,
                   help="Genepop genotype width (default: 6)")
    g.add_argument("--phylip-sites", choices=["snps", "full"], default="snps",
                   help="variable sites only, or full sequences (default: snps)")
    g.add_argument("--phylip-mode", choices=["haploid", "diploid", "population"],
                   default="diploid",
                   help="one row per gene copy, per individual (IUPAC codes), or per "
                        "population (IUPAC codes) (default: diploid)")
    g.add_argument("--phylip-loci", choices=["all", "pi", "fixed"], default="all",
                   help="keep all sites/loci, only parsimony-informative ones (fixed "
                        "difference shared by 2+ rows), or those with any fixed difference "
                        "(default: all)")
    g.add_argument("--phylip-breakpoints", action="store_true",
                   help="mark boundaries between loci with '!'")
    g.add_argument("--phylip-locus-header", action="store_true",
                   help="add a row of locus names after the PHYLIP header")

    g = p.add_argument_group("sequence handling")
    g.add_argument("--clip-5prime", nargs="+", type=dna, default=[], metavar="SEQ",
                   help="remove the first of these sequences found at the 5' end (e.g. a "
                        "restriction site such as TGCAGG)")
    g.add_argument("--clip-3prime", nargs="+", type=dna, default=[], metavar="SEQ",
                   help="remove the first of these sequences found at the 3' end")
    g.add_argument("--keep-ambiguous", action="store_true",
                   help="keep haplotypes with N at SNP positions as distinct alleles "
                        "(default: set those genotypes to missing)")

    g = p.add_argument_group("filters (applied in this order)")
    g.add_argument("--min-coverage", type=nonneg_int, default=0, metavar="N",
                   help="set genotypes with read depth below N to missing (needs --vcf); "
                        "loci absent from the VCF are removed")
    g.add_argument("--coverage-stat", choices=["mean", "min", "max"], default="mean",
                   help="how to combine the depths of a locus' SNPs (default: mean)")
    g.add_argument("--remove-monomorphic", action="store_true",
                   help="remove loci with a single haplotype (checked before and after the "
                        "other filters)")
    g.add_argument("--het-cutoff", type=fraction, metavar="H",
                   help="remove loci with observed heterozygosity >= H that exceed "
                        "Hardy-Weinberg expectations significantly (possible paralogs)")
    g.add_argument("--hwe-alpha", type=fraction, default=0.05, metavar="P",
                   help="significance level for --het-cutoff (default: 0.05)")
    g.add_argument("--min-allele-freq", type=fraction, default=0, metavar="F",
                   help="set genotypes carrying haplotypes with frequency < F to missing")
    g.add_argument("--min-allele-pops", type=fraction, default=0, metavar="F",
                   help="set genotypes carrying haplotypes found in < F of populations to "
                        "missing")
    g.add_argument("--min-pop-locus-freq", type=fraction, default=0, metavar="F",
                   help="remove a locus from any population where < F of individuals are "
                        "genotyped")
    g.add_argument("--min-locus-freq", type=fraction, default=0, metavar="F",
                   help="remove loci genotyped in < F of all individuals")
    g.add_argument("--min-ind-freq", type=fraction, default=0, metavar="F",
                   help="remove individuals genotyped at < F of loci")
    return p


def legacy_prompts(p, args):
    """Reproduce the interactive questions of fasta2genotype 1.x."""
    fasta, whitelist, popmap, vcf, out = args.legacy
    is_na = lambda x: x.lower() == "na"

    def ask(prompt, conv, ok):
        while True:
            try:
                text = input(prompt)
            except EOFError:
                raise InputError("input ended while answering the interactive questions")
            try:
                v = conv(text.strip())
            except ValueError:
                v = None
            if v is not None and ok(v):
                return v
            print("     ** Warning: Not a valid option. **")

    def choose(prompt, n):
        return ask(prompt, int, lambda v: 1 <= v <= n)

    def frac(prompt):
        return ask(prompt, float, lambda v: 0 <= v <= 1)

    def seqs(prompt):
        return ask(prompt, lambda t: t.upper().split(),
                   lambda v: all(re.match("^[ACGT]+$", x) for x in v))

    fmt = {1: "migrate", 2: "arlequin", 3: "diyabc", 4: "lfmm", 5: "phylip", 6: "gphocs",
           7: "treemix"}.get(choose(
               "Output type? [1] Migrate [2] Arlequin [3] DIYABC [4] LFMM [5] Phylip "
               "[6] G-Phocs [7] Treemix [8] Haplotype: ", 8))
    if fmt in ("arlequin", "diyabc"):
        args.title = input("Title of project? : ")
    if fmt == "phylip":
        args.phylip_sites = ["snps", "full"][choose(
            "Use SNPs or full sequences for alignment? [1] SNPs [2] Full Sequences : ", 2) - 1]
        args.phylip_mode = ["haploid", "diploid", "population"][choose(
            "Type of sequences for alignment? [1] Haploid [2] Diploid [3] Population: ", 3) - 1]
        args.phylip_loci = ["pi", "fixed", "all"][choose(
            "Keep only phylogenetically informative (PI) loci, fixed loci, or all loci? "
            "[1] PI [2] Fixed [3] All: ", 3) - 1]
        args.phylip_breakpoints = choose(
            "Flag break points between loci with '!' symbol? [1] Yes [2] No : ", 2) == 1
        args.phylip_locus_header = choose(
            "Insert locus name headers in first row? [1] Yes [2] No : ", 2) == 1
    if fmt == "treemix":
        args.one_snp = choose("How many SNPs to keep per locus? [1] Only one [2] All : ", 2) == 1
    if fmt is None:
        fmt = ["structure", "genepop", "allele-freq", "sambada", "bayescan", "arlequin-hap",
               "genalex"][choose(
                   "Specific output type? [1] Structure [2] Genepop [3] AlleleFreqency "
                   "[4] SamBada [5] Bayescan [6] Arlequin [7] GenAlEx : ", 7) - 1]
        if fmt == "genepop":
            args.genepop_digits = [4, 6][choose(
                "Genepop in four [1] or six [2] digit format? ", 2) - 1]
        if fmt in ("genepop", "arlequin-hap"):
            args.title = input("Title of project? : ")
    args.format = [fmt]

    if not is_na(whitelist) and choose("Loci to use? [1] Whitelist [2] All: ", 2) == 1:
        args.whitelist = whitelist
    if choose("Remove restriction enzyme or adapter sequences? These may bias data. "
              "[1] Yes [2] No: ", 2) == 1:
        args.clip_5prime = seqs("Beginning (5') sequence(s) to remove? (If multiple use "
                                "spaces, if none leave blank): ")
        args.clip_3prime = seqs("Ending (3') sequence(s) to remove? (If multiple use "
                                "spaces, if none leave blank): ")
    if not is_na(vcf):
        args.min_coverage = ask("Coverage Cutoff (number reads for locus)? Use '0' to "
                                "ignore coverage: ", int, lambda v: v >= 0)
        if args.min_coverage:
            args.vcf = vcf
    args.remove_monomorphic = choose("Remove monomorphic loci? [1] Yes [2] No: ", 2) == 1
    if choose("Remove loci with excess heterozygosity? This can remove paralogs. "
              "[1] Yes [2] No: ", 2) == 1:
        args.het_cutoff = frac("Maximum heterozygosity cutoff for removing loci out of "
                               "Hardy-Weinberg? ")
    if choose("Filter for allele frequency? False alleles might bias data. [1] Yes [2] No: ",
              2) == 1:
        args.min_allele_freq = frac("Allele frequency threshold for removal across all "
                                    "individuals? Use '0' to ignore this: ")
        args.min_allele_pops = frac("Frequency of populations containing allele for removal "
                                    "across all individuals? Use '0' to ignore this: ")
    if choose("Filter for missing genotypes? These might bias data. [1] Yes [2] No: ", 2) == 1:
        args.min_locus_freq = frac("Locus frequency threshold for locus removal across all "
                                   "individuals? Use '0' to ignore this: ")
        args.min_pop_locus_freq = frac("Population frequency threshold for locus removal "
                                       "across each population? Use '0' to ignore this: ")
        args.min_ind_freq = frac("Individual frequency threshold for individual removal "
                                 "across all loci? Use '0' to ignore this: ")
    args.fasta, args.popmap, args.out = fasta, popmap, out
    return args


def parse_args(argv=None):
    p = build_parser()
    args = p.parse_args(argv)
    args.legacy_names = False
    if args.legacy:
        if len(args.legacy) != 5 or args.fasta or args.popmap or args.out or args.format:
            p.error("use either the 5 old-style positional arguments "
                    "([fasta] [whitelist] [population file] [VCF] [output name]) "
                    "or the -f/-p/-o/-F options, not both")
        args = legacy_prompts(p, args)
        args.legacy_names = True
    else:
        missing = [flag for flag, v in (("-f/--fasta", args.fasta), ("-p/--popmap", args.popmap),
                                        ("-o/--out", args.out), ("-F/--format", args.format))
                   if not v]
        if missing:
            p.error("the following arguments are required: " + ", ".join(missing))
    if args.min_coverage and not args.vcf:
        p.error("--min-coverage needs --vcf")
    if args.vcf and not args.min_coverage:
        warn("--vcf is only used with --min-coverage; ignoring it")
    args.format = list(FORMATS) if "all" in args.format else list(dict.fromkeys(args.format))
    return args


def main(argv=None):
    global QUIET
    args = parse_args(argv)
    QUIET = args.quiet
    info("fasta2genotype %s -- please cite:\n%s\n" % (__version__, CITATION))
    try:
        samples = read_popmap(args.popmap, args.popmap_format)
        whitelist = read_locus_list(args.whitelist) if args.whitelist else None
        blacklist = read_locus_list(args.blacklist) if args.blacklist else None
        ds = read_fasta(args.fasta, samples, whitelist, blacklist, args.clip_5prime,
                        args.clip_3prime, args.keep_ambiguous)
        info("  %d individuals in %d populations, %d loci" % (
            len(ds.order), len(ds.pops), len(ds.loci)))
        depths = None
        if args.min_coverage:
            depths = read_vcf_depths(args.vcf, ds, args.coverage_stat)
        run_filters(ds, args, depths)
        info("After filtering: %d individuals in %d populations, %d loci" % (
            len(ds.order), len(ds.pops), len(ds.loci)))

        namer = OutputNamer(args.out, legacy=args.legacy_names)
        written = []
        for fmt in args.format:
            info("Writing %s..." % fmt)
            written.extend(FORMATS[fmt][0](ds, namer, args))
    except InputError as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 1
    info("Done. Wrote:\n  " + "\n  ".join(written))
    return 0


if __name__ == "__main__":
    sys.exit(main())
