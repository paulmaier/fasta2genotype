"""Tests for fasta2genotype.py.  Run from the repository root with:

    python3 -m unittest discover -s tests -v
"""

import builtins
import contextlib
import filecmp
import gzip
import io
import math
import os
import sys
import tempfile
import unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import fasta2genotype as f2g  # noqa: E402

EX1 = os.path.join(ROOT, "examples", "stacks_v1")
EX2 = os.path.join(ROOT, "examples", "stacks_v2")

A = "ACGTACGTAC"   # three haplotypes of one 10-bp locus
B = "ACGTTCGTAC"   # differs from A at column 4
C = "ACGTACGTAG"   # differs from A at column 9


def run(argv):
    """Run main() quietly; return (exit code, stderr text)."""
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        code = f2g.main(argv)
    return code, err.getvalue()


class TempDirTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        f2g.QUIET = True

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def read(self, name):
        with open(os.path.join(self.tmp, name)) as fh:
            return fh.read()

    def v2_fasta(self, genotypes, name="in.fa"):
        """genotypes: {locus: {sample: (seq0, seq1)}} -> Stacks 2 style FASTA."""
        lines = ["# Stacks version 2.66; test"]
        for locus, by_sample in genotypes.items():
            for i, (sample, seqs) in enumerate(by_sample.items()):
                for allele, seq in enumerate(seqs):
                    lines.append(">CLocus_%s_Sample_%d_Locus_%s_Allele_%d [%s; chr1, 100, +]"
                                 % (locus, i + 1, locus, allele, sample))
                    lines.append(seq)
        return self.write(name, "\n".join(lines) + "\n")

    def popmap(self, pops, name="popmap.tsv"):
        return self.write(name, "".join("%s\t%s\n" % kv for kv in pops.items()))

    def dataset(self, genotypes, pops, **kw):
        samples = f2g.read_popmap(self.popmap(pops))
        return f2g.read_fasta(self.v2_fasta(genotypes), samples, **kw)


class TestUtilities(unittest.TestCase):
    def test_chi2_sf_closed_forms(self):
        for x in (0.1, 1.0, 3.84, 10.0, 40.0):
            self.assertAlmostEqual(f2g.chi2_sf(x, 1), math.erfc(math.sqrt(x / 2)), places=10)
            self.assertAlmostEqual(f2g.chi2_sf(x, 2), math.exp(-x / 2), places=10)
            self.assertAlmostEqual(f2g.chi2_sf(x, 4), math.exp(-x / 2) * (1 + x / 2), places=10)
        self.assertEqual(f2g.chi2_sf(0, 3), 1.0)
        self.assertEqual(f2g.chi2_sf(5, 0), 1.0)

    def test_iupac(self):
        self.assertEqual(f2g.iupac("AG"), "R")
        self.assertEqual(f2g.iupac("ACG"), "V")
        self.assertEqual(f2g.iupac("AN"), "A")
        self.assertEqual(f2g.iupac("NN"), "N")
        self.assertEqual(f2g.iupac("RC"), "V")
        self.assertEqual(f2g.iupac_merge(["ACGT", "ACGA"]), "ACGW")

    def test_fixed_difference(self):
        self.assertTrue(f2g.has_fixed_difference(["A", "G"]))
        self.assertTrue(f2g.has_fixed_difference(["M", "K"]))
        self.assertFalse(f2g.has_fixed_difference(["A", "R"]))
        self.assertFalse(f2g.has_fixed_difference(["A", "N", "A"]))

    def test_locus_sort(self):
        self.assertEqual(sorted(["110", "45", "9", "x"], key=f2g.locus_sort_key),
                         ["9", "45", "110", "x"])


class TestInput(TempDirTest):
    def test_stacks2_fasta(self):
        ds = self.dataset({"1": {"s1": (A, B), "s2": (A, A)}}, {"s1": "p1", "s2": "p1"})
        loc = ds.loci["1"]
        self.assertEqual(loc.length, 10)
        self.assertNotEqual(loc.g0[0], loc.g1[0])   # heterozygote
        self.assertEqual(loc.g0[1], loc.g1[1])      # homozygote written twice
        self.assertEqual(len(loc.alleles), 2)

    def test_stacks1_fasta_and_legacy_popmap(self):
        fa = self.write("v1.fa", ">CLocus_7_Sample_3_Locus_1_Allele_0\n%s\n"
                                 ">CLocus_7_Sample_3_Locus_1_Allele_1\n%s\n"
                                 ">CLocus_7_Sample_4_Locus_9_Allele_0\n%s\n"
                                 ">CLocus_7_Sample_4_Locus_9_Allele_2\n%s\n" % (A, B, A, C))
        pm = self.write("pm.tsv", "SampleID\tIndividualID\tPopulationID\n3\tind3\tP\n4\tind4\tP\n")
        samples = f2g.read_popmap(pm)
        self.assertEqual([s.name for s in samples], ["ind3", "ind4"])
        ds = f2g.read_fasta(fa, samples)
        loc = ds.loci["7"]
        self.assertEqual((loc.g0[1], loc.g1[1]), (0, 0))  # allele 0 only -> homozygote
        self.assertEqual(len(loc.alleles), 2)              # allele 2 ignored

    def test_unknown_header_is_an_error(self):
        fa = self.write("bad.fa", ">locus1\nACGT\n")
        samples = f2g.read_popmap(self.popmap({"s1": "p"}))
        with self.assertRaises(f2g.InputError):
            f2g.read_fasta(fa, samples)

    def test_gzipped_input(self):
        path = self.v2_fasta({"1": {"s1": (A, B)}})
        with open(path, "rb") as src, gzip.open(path + ".gz", "wb") as dst:
            dst.write(src.read())
        ds = f2g.read_fasta(path + ".gz", f2g.read_popmap(self.popmap({"s1": "p"})))
        self.assertEqual(len(ds.loci), 1)

    def test_whitelist_order_and_blacklist(self):
        geno = {lid: {"s1": (A, B)} for lid in ("1", "2", "3", "10")}
        ds = self.dataset(geno, {"s1": "p"}, whitelist=["10", "2", "1"], blacklist=["1"])
        self.assertEqual(ds.locus_order, ["10", "2"])
        ds = self.dataset(geno, {"s1": "p"})
        self.assertEqual(ds.locus_order, ["1", "2", "3", "10"])

    def test_clipping(self):
        ds = self.dataset({"1": {"s1": ("TGCAGG" + A, "TGCAGG" + B)}}, {"s1": "p"},
                          clip5=["CGG", "TGCAGG"], clip3=["AC"])
        self.assertEqual(sorted(ds.loci["1"].alleles), sorted([A[:-2], B[:-2]]))

    def test_ambiguous_snp_sets_genotype_missing(self):
        amb = A[:4] + "N" + A[5:]             # N at the SNP column 4
        geno = {"1": {"s1": (A, B), "s2": (amb, A), "s3": (B, B)}}
        pops = {"s1": "p", "s2": "p", "s3": "p"}
        loc = self.dataset(geno, pops).loci["1"]
        self.assertEqual(loc.g0[1], -1)
        self.assertEqual(len({loc.g0[0], loc.g1[0]}), 2)
        loc = self.dataset(geno, pops, keep_ambiguous=True).loci["1"]
        self.assertGreaterEqual(loc.g0[1], 0)
        self.assertEqual(len(loc.alleles), 3)

    def test_n_at_invariant_column_does_not_split_alleles(self):
        amb = A[:1] + "N" + A[2:]             # column 1 is invariant
        loc = self.dataset({"1": {"s1": (A, B), "s2": (amb, amb)}},
                           {"s1": "p", "s2": "p"}).loci["1"]
        self.assertEqual(len(loc.alleles), 2)
        self.assertEqual(loc.g0[1], loc.g0[0])
        self.assertTrue(all(q[1] == "N" for q in loc.alleles))

    def test_vcf_depths(self):
        ds = self.dataset({"5": {"s1": (A, B), "s2": (A, A)}}, {"s1": "p", "s2": "p"})
        vcf = self.write("in.vcf", "##fileformat=VCFv4.2\n"
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
                         "c\t1\t5:4:+\tA\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:10:5,5\t0/0:30:30,0\n"
                         "c\t2\t5:9:+\tC\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:20:5,5\t./.\n")
        self.assertEqual(f2g.read_vcf_depths(vcf, ds, "mean")["5"], [15, 15])
        self.assertEqual(f2g.read_vcf_depths(vcf, ds, "min")["5"], [10, 0])
        self.assertEqual(f2g.read_vcf_depths(vcf, ds, "max")["5"], [20, 30])


class TestFilters(TempDirTest):
    def args(self, **kw):
        ns = f2g.build_parser().parse_args([])
        for k, v in kw.items():
            setattr(ns, k, v)
        return ns

    def test_monomorphic(self):
        ds = self.dataset({"1": {"s1": (A, A), "s2": (A, A)}, "2": {"s1": (A, B), "s2": (A, A)}},
                          {"s1": "p", "s2": "p"})
        f2g.run_filters(ds, self.args(remove_monomorphic=True))
        self.assertEqual(list(ds.loci), ["2"])

    def test_allele_frequency(self):
        geno = {"1": {"s%d" % i: (A, A) for i in range(9)}}
        geno["1"]["s9"] = (A, C)              # C is 1 of 20 gene copies
        ds = self.dataset(geno, {"s%d" % i: "p" for i in range(10)})
        f2g.filter_alleles(ds, 0.06, 0)
        self.assertEqual(ds.loci["1"].g0[9], -1)
        self.assertEqual(ds.loci["1"].n_genotyped(ds.order), 9)

    def test_allele_population_fraction(self):
        geno = {"1": {"a1": (A, B), "b1": (A, A), "c1": (A, A), "d1": (A, A)}}
        ds = self.dataset(geno, {"a1": "A", "b1": "B", "c1": "C", "d1": "D"})
        f2g.filter_alleles(ds, 0, 0.3)        # B is in 1 of 4 populations
        self.assertEqual(ds.loci["1"].g0[0], -1)

    def test_missing_data(self):
        geno = {"1": {"s1": (A, B), "s2": (A, A), "s3": (A, A), "s4": (B, B)},
                "2": {"s1": (A, B)},
                "3": {"s1": (A, B), "s2": (A, B), "s3": (A, B)}}
        ds = self.dataset(geno, {"s1": "p", "s2": "p", "s3": "q", "s4": "q"})
        f2g.filter_missing(ds, 0, 0.5, 0)
        self.assertEqual(sorted(ds.loci), ["1", "3"])
        f2g.filter_missing(ds, 0.6, 0, 0)     # locus 3 has 1 of 2 in pop q
        self.assertEqual(ds.loci["3"].n_genotyped(ds.members["q"]), 0)
        f2g.filter_missing(ds, 0, 0, 0.75)    # s3 and s4 now have 1 of 2 loci
        self.assertEqual([ds.samples[s].name for s in ds.order], ["s1", "s2"])
        self.assertEqual(ds.pops, ["p"])

    def test_heterozygosity_excess_removed_deficit_kept(self):
        pops = {"s%d" % i: "p" for i in range(40)}
        excess = {"s%d" % i: (A, B) for i in range(40)}            # all heterozygous
        deficit = {"s%d" % i: ((A, A) if i < 20 else (B, B)) for i in range(40)}
        ds = self.dataset({"1": excess, "2": deficit}, pops)
        f2g.filter_heterozygosity(ds, 0.0, 0.05)
        self.assertEqual(list(ds.loci), ["2"])

    def test_individuals_without_data_are_dropped(self):
        ds = self.dataset({"1": {"s1": (A, B), "s2": (A, A)}, "2": {"s2": (A, C)}},
                          {"s1": "p", "s2": "q"})
        f2g.filter_missing(ds, 0, 0.9, 0)     # removes locus 2 (1 of 2 individuals)
        f2g.run_filters(ds, self.args())
        self.assertEqual(len(ds.order), 2)
        ds.loci["1"].set_missing(0)
        f2g.run_filters(ds, self.args())
        self.assertEqual([ds.samples[s].name for s in ds.order], ["s2"])
        self.assertEqual(ds.pops, ["q"])


class TestWriters(TempDirTest):
    GENO = {"1": {"s1": (A, B), "s2": (B, B), "s3": (A, A)},
            "2": {"s1": (A, A), "s3": (A, C)}}
    POPS = {"s1": "west", "s2": "west", "s3": "east"}

    def write_format(self, fmt, *extra):
        out = os.path.join(self.tmp, "out")
        code, err = run(["-q", "-f", self.v2_fasta(self.GENO), "-p", self.popmap(self.POPS),
                         "-o", out, "-F", fmt] + list(extra))
        self.assertEqual(code, 0, err)
        return out

    def test_structure(self):
        self.write_format("structure")
        rows = [l.split("\t") for l in self.read("out.str").splitlines()]
        self.assertEqual(rows[0], ["", "", "1", "2"])
        self.assertEqual(rows[1], ["s1", "1", "1", "1"])
        self.assertEqual(rows[2], ["s1", "1", "2", "1"])
        self.assertEqual(rows[3], ["s2", "1", "2", "-9"])
        self.assertEqual(rows[5][:2], ["s3", "2"])
        self.assertIn("2\teast", self.read("out.str.pops.tsv"))

    def test_genepop(self):
        self.write_format("genepop", "--genepop-digits", "4", "--title", "T")
        lines = self.read("out.gen").splitlines()
        self.assertEqual(lines[:5], ["T", "1", "2", "Pop", "s1 ,  0102 0101"])
        self.assertEqual(lines[5], "s2 ,  0202 0000")

    def test_treemix_counts_homozygotes_twice(self):
        self.write_format("treemix")
        with gzip.open(os.path.join(self.tmp, "out.treemix.gz"), "rt") as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[0], "west east")
        self.assertEqual(lines[1], "1,3 2,0")   # locus 1 column 4: A/T
        self.assertEqual(lines[2], "2,0 1,1")   # locus 2 column 9: C/G

    def test_migrate(self):
        self.write_format("migrate")
        lines = self.read("out.migrate").splitlines()
        self.assertEqual(lines[0], "2\t2")
        self.assertEqual(lines[1], "10\t10")
        self.assertEqual(lines[2], "4\tPop_west")
        self.assertEqual(lines[3], "s1a       " + A)
        self.assertIn("s2b       " + "?" * 10, lines)

    def test_migrate_long_names(self):
        self.GENO = {"1": {"sample_long_name": (A, B), "s2": (A, A)}}
        self.POPS = {"sample_long_name": "p", "s2": "p"}
        self.write_format("migrate")
        self.assertIn("I1a       " + A, self.read("out.migrate"))
        self.assertIn("I1\tsample_long_name\tp", self.read("out.migrate.names.tsv"))

    def test_phylip_modes(self):
        self.write_format("phylip")
        lines = self.read("out.phy").splitlines()
        self.assertEqual(lines[0], "3 2")
        self.assertEqual(lines[1], "s1        WC")
        self.assertEqual(lines[2], "s2        TN")
        self.write_format("phylip", "--phylip-mode", "population", "--phylip-sites", "full")
        lines = self.read("out.phy").splitlines()
        self.assertEqual(lines[0], "2 20")
        self.assertEqual(lines[1], "west      ACGTWCGTAC" + A)
        self.write_format("phylip", "--phylip-mode", "haploid", "--phylip-breakpoints")
        self.assertEqual(self.read("out.phy").splitlines()[1], "s1a       A!C!")

    def test_phylip_informative_sites(self):
        self.write_format("phylip", "--phylip-loci", "pi", "--phylip-mode", "haploid")
        lines = self.read("out.phy").splitlines()
        self.assertEqual(lines[0], "6 1")      # only locus 1 column 4 is informative
        self.assertEqual(lines[1], "s1a       A")
        out = os.path.join(self.tmp, "pi")
        code, err = run(["-q", "-f", self.v2_fasta(self.GENO), "-p", self.popmap(self.POPS),
                         "-o", out, "-F", "phylip", "--phylip-loci", "pi"])
        self.assertEqual(code, 1)              # no site is informative among individuals
        self.assertIn("no loci qualify", err)

    def test_gphocs_lists_only_genotyped(self):
        self.write_format("gphocs")
        text = self.read("out.gphocs")
        self.assertTrue(text.startswith("2\n\n1\t3\t10\n"))
        self.assertIn("\n2\t2\t10\ns1\t%s\ns3\t%s\n" % (A, A[:9] + "S"), text)

    def test_lfmm(self):
        self.write_format("lfmm")
        self.assertEqual(self.read("out.lfmm").splitlines(), ["1 0", "2 9", "0 1"])
        self.assertIn("1_4\t1\t4\tA\tT", self.read("out.lfmm.snps.tsv"))

    def test_bayescan_and_allele_freq(self):
        self.write_format("bayescan", "allele-freq")
        self.assertEqual(self.read("out.bayescan.txt").splitlines()[4:7],
                         ["[pop]=1", "1\t4\t2\t1\t3", "2\t2\t2\t2\t0"])
        freq = self.read("out.allelefreq.tsv").splitlines()
        self.assertEqual(freq[0], "\t1_1\t1_2\t2_1\t2_2")
        self.assertEqual(freq[2], "east\t1.00000\t0.00000\t0.50000\t0.50000")

    def test_sambada_one_allele_per_gene_copy(self):
        self.write_format("sambada")
        lines = [l.split("\t") for l in self.read("out.sambada.txt").splitlines()]
        self.assertEqual(lines[1], ["s1a", "1", "0", "1", "0"])
        self.assertEqual(lines[2], ["s1b", "0", "1", "1", "0"])
        self.assertEqual(lines[3], ["s2a", "0", "1", "NaN", "NaN"])

    def test_arlequin(self):
        self.write_format("arlequin", "arlequin-hap", "--title", "My project")
        text = self.read("out.arp")
        self.assertIn('Title="My project"', text)
        self.assertIn('SampleName="west"\n\t\tSampleSize=2', text)
        self.assertIn("s2\t1\t%s\t%s\n\t\t%s\t%s\n" % (B, "?" * 10, B, "?" * 10), text)
        self.assertIn("s2\t1\t2\t?\n\t\t2\t?\n", self.read("out.hap.arp"))

    def test_genalex_and_diyabc(self):
        self.write_format("genalex", "diyabc", "diyabc-snp", "--title", "T")
        gx = self.read("out.genalex.txt").splitlines()
        self.assertEqual(gx[0], "2\t3\t2\t2\t1")
        self.assertEqual(gx[3], "s1\twest\t1\t2\t1\t1")
        dy = self.read("out.diyabc").splitlines()
        self.assertEqual(dy[:4], ["T <NM=1.0NF>", "1\t<A>", "2\t<A>", "Pop"])
        self.assertEqual(dy[5], "s2\t,\t<[%s][%s]>\t<[][]>" % (B, B))
        snp = self.read("out.diyabc.snp").splitlines()
        self.assertEqual(snp[1], "IND\tSEX\tPOP\tA\tA")
        self.assertEqual(snp[3], "s2\t9\twest\t2\t9")


class TestCommandLine(TempDirTest):
    def test_examples_v1_and_v2_give_identical_output(self):
        a, b = os.path.join(self.tmp, "v1"), os.path.join(self.tmp, "v2")
        opts = ["-q", "-F", "all", "--remove-monomorphic", "--title", "ex"]
        self.assertEqual(run(opts + ["-f", os.path.join(EX1, "batch_1.fa"),
                                     "-p", os.path.join(EX1, "popmap_legacy.tsv"),
                                     "-o", os.path.join(a, "ex")])[0], 0)
        self.assertEqual(run(opts + ["-f", os.path.join(EX2, "populations.samples.fa"),
                                     "-p", os.path.join(EX2, "popmap.tsv"),
                                     "-o", os.path.join(b, "ex")])[0], 0)
        names = sorted(os.listdir(a))
        self.assertEqual(names, sorted(os.listdir(b)))
        self.assertEqual(len(names), 21)
        match, mismatch, errors = filecmp.cmpfiles(a, b, names, shallow=False)
        # the gzip header stores the file name, so compare treemix content separately
        mismatch = [m for m in mismatch if not m.endswith(".gz")]
        self.assertEqual(mismatch + errors, [])

    def test_coverage_filter_v1_and_v2_agree(self):
        outs = []
        for d, fa, pm, vcf in ((EX1, "batch_1.fa", "popmap_legacy.tsv", "batch_1.vcf"),
                               (EX2, "populations.samples.fa", "popmap.tsv",
                                "populations.snps.vcf")):
            out = os.path.join(self.tmp, os.path.basename(d))
            code, err = run(["-q", "-f", os.path.join(d, fa), "-p", os.path.join(d, pm),
                             "--vcf", os.path.join(d, vcf), "--min-coverage", "50",
                             "-o", out, "-F", "migrate"])
            self.assertEqual(code, 0, err)
            outs.append(open(out + ".migrate").read())
        self.assertEqual(outs[0], outs[1])

    def test_legacy_interactive_mode(self):
        answers = iter(["8", "1", "1", "2", "1", "2", "2", "2"])
        out = os.path.join(self.tmp, "legacy")
        with mock.patch.object(builtins, "input", lambda prompt="": next(answers)), \
                contextlib.redirect_stdout(io.StringIO()):
            code, err = run(["-q", os.path.join(EX1, "batch_1.fa"),
                             os.path.join(EX1, "whitelist.txt"),
                             os.path.join(EX1, "popmap_legacy.tsv"), "NA", out])
        self.assertEqual(code, 0, err)
        self.assertTrue(os.path.exists(out + ".out"))
        self.assertTrue(os.path.exists(out + "_pops.out"))
        # loci follow the whitelist order
        self.assertEqual(open(out + ".out").readline().split(), ["110", "112", "134", "45", "91"])

    def test_argument_errors(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                f2g.parse_args(["-f", "x.fa"])
            with self.assertRaises(SystemExit):
                f2g.parse_args(["-f", "x", "-p", "y", "-o", "z", "-F", "structure",
                                "--min-coverage", "5"])
            with self.assertRaises(SystemExit):
                f2g.parse_args(["-f", "x", "-p", "y", "-o", "z", "-F", "structure",
                                "--min-locus-freq", "1.5"])

    def test_missing_file_is_reported(self):
        code, err = run(["-q", "-f", "nope.fa", "-p", "nope.tsv", "-o",
                         os.path.join(self.tmp, "x"), "-F", "structure"])
        self.assertEqual(code, 1)
        self.assertIn("cannot open", err)


if __name__ == "__main__":
    unittest.main()
