#!/usr/bin/env python3
"""Compare Rust fragment BEDs against Biopython, which is test-only."""
import os
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from Bio import Restriction
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic')).resolve()
TABLE = ROOT / 'data/restriction-enzymes.tsv'


def expected_bed(seq, names):
    boundaries = sorted({cut + 1 for name in names
                         for cut in getattr(Restriction, name).search(Seq(seq))
                         if 0 < cut + 1 < len(seq)})
    points = [0, *boundaries, len(seq)]
    return ''.join(f'chr1\t{a}\t{b}\n' for a, b in zip(points, points[1:]))


class DigestComparison(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='rust-hic-digest-')
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def digest(self, seq, enzyme):
        fasta, sizes, output = (self.root / x for x in ('reference.fa', 'chrom.sizes', 'out.bed'))
        fasta.write_text('>chr1\n' + seq + '\n')
        sizes.write_text(f'chr1\t{len(seq)}\n')
        output.unlink(missing_ok=True)
        result = subprocess.run([str(BIN), 'digest', '-r', str(fasta), '--chromsizes', str(sizes),
                                 '-e', enzyme, '-o', str(output)], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return output.read_text()

    def test_all_enzyme_definitions_both_orientations(self):
        concrete = str.maketrans({'R':'A','Y':'C','S':'C','W':'A','K':'G','M':'A',
                                  'B':'C','D':'A','H':'A','V':'A','N':'A','U':'T'})
        rows = (line.split('\t') for line in TABLE.read_text().splitlines()
                if line and not line.startswith('#') and not line.startswith('name\t'))
        for fields in rows:
            name, site = fields[:2]
            with self.subTest(enzyme=name):
                motif = site.translate(concrete)
                reverse = str(Seq(motif).reverse_complement())
                seq = 'T' * 80 + motif + 'A' * 80 + reverse + 'T' * 80
                self.assertEqual(self.digest(seq, name), expected_bed(seq, [name]))

    def test_cocktail_deduplicates_shared_cuts(self):
        seq = 'AAAAAGATCAAAAA'
        self.assertEqual(self.digest(seq, 'MboI-DpnII'), expected_bed(seq, ['MboI', 'DpnII']))

    def test_chromosome_edges_and_ambiguous_bases(self):
        for enzyme, motif in [('MboI', 'GATC'), ('HinfI', 'GANTC'),
                              ('BsaI', 'GAGACC'), ('BcgI', 'CGAAAAAAATGC')]:
            for padding in (0, 1, 2, 5, 15, 30):
                for seq in (motif + 'A' * padding, 'A' * padding + motif,
                            'N' * padding + motif.lower() + 'N' * padding):
                    with self.subTest(enzyme=enzyme, sequence=seq):
                        self.assertEqual(self.digest(seq, enzyme), expected_bed(seq, [enzyme]))

    def test_no_overwrite_or_invalid_enzyme(self):
        seq = 'A' * 20 + 'GATC' + 'A' * 20
        self.digest(seq, 'MboI')
        fasta, sizes, output = (self.root / x for x in ('reference.fa', 'chrom.sizes', 'out.bed'))
        before = output.read_bytes()
        args = [str(BIN), 'digest', '-r', str(fasta), '--chromsizes', str(sizes), '-o', str(output)]
        self.assertNotEqual(subprocess.run([*args, '-e', 'MboI'], capture_output=True).returncode, 0)
        self.assertEqual(output.read_bytes(), before)
        output.unlink()
        self.assertNotEqual(subprocess.run([*args, '-e', 'not-an-enzyme'], capture_output=True).returncode, 0)
        self.assertFalse(output.exists())
        sizes.write_text(f'chr1\t{len(seq)}\nmissing\t10\n')
        self.assertNotEqual(subprocess.run([*args, '-e', 'MboI'], capture_output=True).returncode, 0)
        self.assertFalse(output.exists())

    def test_reuse_generated_bed_in_workflow(self):
        seq = 'A' * 20 + 'GATC' + 'A' * 20
        self.digest(seq, 'MboI')
        fasta, sizes, bed = (self.root / x for x in ('reference.fa', 'chrom.sizes', 'out.bed'))
        pairs = self.root / 'input.pairs'
        pairs.write_text('## pairs format v1.0.0\n#sorted: chr1-chr2-pos1-pos2\n'
                         '#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n'
                         'read1\tchr1\t5\tchr1\t35\t+\t-\tUU\n')
        results = []
        for mode in ('digest', 'bed'):
            output = self.root / ('result-' + mode)
            config = {'assembly': 'test', 'chromsizes': str(sizes), 'input_pairs': [str(pairs)],
                      'output_dir': str(output), 'scratch_dir': str(self.root), 'keep_pairs': True,
                      'resolutions': [10], 'balance_max_iterations': 1}
            if mode == 'digest':
                config['digest'] = {'enzyme': 'MboI', 'reference': str(fasta)}
            else:
                config['fragments'] = str(bed)
            path = self.root / (mode + '.json')
            path.write_text(json.dumps(config))
            result = subprocess.run([str(BIN), 'run', str(path)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(next((output / 'qc').glob('*.report.json')).read_text())
            self.assertEqual(report['restriction']['mode'], 'digest' if mode == 'digest' else 'bed')
            if mode == 'digest':
                self.assertEqual(report['restriction']['enzyme'], 'MboI')
            results.append(next((output / 'pairs').glob('*.pairs')).read_bytes())
        self.assertEqual(*results)


if __name__ == '__main__':
    unittest.main(verbosity=2)
