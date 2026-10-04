"""Every row worked out by hand over a graph small enough to check on paper.

Nothing downstream of the converter would catch a flipped chain's coordinates
and CIGAR direction, or where a chain breaks, so these tests pin each one.
"""

import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gfa_to_pairwise_paf.py')

# GRCh38 chrA walks 1 2 4 5 6 7 8 9 8 10 (59 bp; node 8 twice, at 33 and 49).
# HG01109 ctgA arrives as two W pieces: the first takes the SNP allele 3 for 2,
# skips 5, inserts 11, and reaches 8 the first time; the second starts at
# offset 50 on 9 and reaches 8 the second time. HG01123 ctgB traverses 7 and 5
# backwards, skipping 6 between them. HG00097 ctgC is one shared node.
NODES = [
    ('1', 'ACGTACGTAC'),
    ('2', 'G'),
    ('3', 'T'),
    ('4', 'CCCCC'),
    ('5', 'AAA'),
    ('6', 'GGGGGGGG'),
    ('7', 'TTTTTT'),
    ('8', 'CAGCAG'),
    ('9', 'AAAAAAAAAA'),
    ('10', 'CGCG'),
    ('11', 'TTAA'),
]
WALKS = (
    'W\tGRCh38\t0\tchrA\t0\t59\t>1>2>4>5>6>7>8>9>8>10\n'
    'W\tHG01109\t1\tctgA\t0\t50\t>1>3>4>6>11>7>8>9\n'
    'W\tHG01109\t1\tctgA\t50\t70\t>9>8>10\n'
    'W\tHG01123\t1\tctgB\t0\t23\t>1<7<5>10\n'
    'W\tHG00097\t1\tctgC\t0\t10\t>1\n'
)
PATHS = (
    'P\tGRCh38#0#chrA\t1+,2+,4+,5+,6+,7+,8+,9+,8+,10+\t*\n'
    'P\tHG01109#1#ctgA\t1+,3+,4+,6+,11+,7+,8+,9+,9+,8+,10+\t*\n'
    'P\tHG01123#1#ctgB\t1+,7-,5-,10+\t*\n'
    'P\tHG00097#1#ctgC\t1+\t*\n'
)
HEADER = 'H\tVN:Z:1.1\n' + ''.join(f'S\t{i}\t{seq}\n' for i, seq in NODES) + 'L\t1\t+\t2\t+\t0M\n'
WITH_WALKS = HEADER + WALKS
WITH_PATHS = HEADER + PATHS


def row(*fields):
    return [str(f) for f in fields]


# qname qlen qstart qend strand tname tlen tstart tend matches blocklen mapq cg
CTGA_FIRST_PIECE = row('HG01109#1#ctgA', 70, 0, 50, '+', 'GRCh38#0#chrA', 59, 0, 49, 45, 53, 255, 'cg:Z:10=1X5=3D8=4I22=')
CTGA_SECOND_PIECE = row('HG01109#1#ctgA', 70, 50, 70, '+', 'GRCh38#0#chrA', 59, 39, 59, 20, 20, 255, 'cg:Z:20=')
CTGB_FORWARD = row('HG01123#1#ctgB', 23, 0, 10, '+', 'GRCh38#0#chrA', 59, 0, 10, 10, 10, 255, 'cg:Z:10=')
# 7 then 5 backwards: reference 16..33 with the 8 bp of node 6 deleted, and
# the CIGAR reads along the reference, so node 5 comes first
CTGB_INVERTED = row('HG01123#1#ctgB', 23, 10, 19, '-', 'GRCh38#0#chrA', 59, 16, 33, 9, 17, 255, 'cg:Z:3=8D6=')
CTGB_TAIL = row('HG01123#1#ctgB', 23, 19, 23, '+', 'GRCh38#0#chrA', 59, 55, 59, 4, 4, 255, 'cg:Z:4=')
CTGC = row('HG00097#1#ctgC', 10, 0, 10, '+', 'GRCh38#0#chrA', 59, 0, 10, 10, 10, 255, 'cg:Z:10=')
ALL_ROWS = [CTGC, CTGA_FIRST_PIECE, CTGA_SECOND_PIECE, CTGB_FORWARD, CTGB_INVERTED, CTGB_TAIL]


def convert(gfa, args):
    return subprocess.run(
        [sys.executable, SCRIPT, *args],
        input=gfa,
        capture_output=True,
        text=True,
    )


def paf_rows(stdout):
    rows = [line.split('\t') for line in stdout.split('\n') if line]
    return sorted(rows, key=lambda r: (r[0], int(r[2])))


def identical(q, qlen, qs, qe, strand, ts, te):
    n = qe - qs
    return row(q, qlen, qs, qe, strand, 'GRCh38#0#chrA', 59, ts, te, n, n, 255, f'cg:Z:{n}=')


def with_qlen(fields, qlen):
    return [fields[0], str(qlen), *fields[2:]]


class GfaToPairwisePaf(unittest.TestCase):
    def test_every_row_by_hand(self):
        run = convert(WITH_WALKS, ['--reference', 'GRCh38#0'])
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(paf_rows(run.stdout), ALL_ROWS)
        self.assertRegex(run.stderr, r'HG01109#1: 2 walks, 9 anchors -> 2 chains, 65 bp =, 73 columns')
        self.assertRegex(run.stderr, r'11 nodes, 10 GRCh38#0 steps on 1 walks')

    def test_p_lines_and_bare_reference(self):
        run = convert(WITH_PATHS, ['--reference', 'GRCh38'])
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(paf_rows(run.stdout), ALL_ROWS)

    def test_max_gap_breaks_a_chain(self):
        run = convert(WITH_WALKS, ['--reference', 'GRCh38', '--queries', 'HG01109#1,HG01123#1', '--max-gap', '2'])
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(
            paf_rows(run.stdout),
            [
                row('HG01109#1#ctgA', 70, 0, 16, '+', 'GRCh38#0#chrA', 59, 0, 16, 15, 16, 255, 'cg:Z:10=1X5='),
                identical('HG01109#1#ctgA', 70, 16, 24, '+', 19, 27),
                identical('HG01109#1#ctgA', 70, 28, 50, '+', 27, 49),
                CTGA_SECOND_PIECE,
                CTGB_FORWARD,
                identical('HG01123#1#ctgB', 23, 10, 16, '-', 27, 33),
                identical('HG01123#1#ctgB', 23, 16, 19, '-', 16, 19),
                CTGB_TAIL,
            ],
        )

    def test_no_x_writes_insertion_then_deletion(self):
        run = convert(WITH_WALKS, ['--reference', 'GRCh38', '--queries', 'HG01109#1', '--no-x'])
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(
            paf_rows(run.stdout),
            [
                row('HG01109#1#ctgA', 70, 0, 50, '+', 'GRCh38#0#chrA', 59, 0, 49, 45, 54, 255, 'cg:Z:10=1I1D5=3D8=4I22='),
                CTGA_SECOND_PIECE,
            ],
        )

    def test_min_block_and_queries(self):
        run = convert(WITH_WALKS, ['--reference', 'GRCh38', '--queries', 'HG01123#1', '--min-block', '15'])
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(paf_rows(run.stdout), [CTGB_INVERTED])

    def test_chrom_sizes_per_query(self):
        with tempfile.TemporaryDirectory() as tmp:
            lengths = os.path.join(tmp, 'lengths.txt')
            with open(lengths, 'w') as f:
                f.write('ctgA\t100\n')
            sizes = os.path.join(tmp, 'sizes')
            run = convert(
                WITH_WALKS,
                ['--reference', 'GRCh38', '--chrom-sizes-dir', sizes, '--contig-lengths', lengths],
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            for name, expected in [
                ('HG01109.1.chrom.sizes', 'ctgA\t100\n'),
                ('HG01123.1.chrom.sizes', 'ctgB\t23\n'),
                ('HG00097.1.chrom.sizes', 'ctgC\t10\n'),
            ]:
                with open(os.path.join(sizes, name)) as f:
                    self.assertEqual(f.read(), expected)
        self.assertEqual(
            [r for r in paf_rows(run.stdout) if r[0] == 'HG01109#1#ctgA'],
            [with_qlen(CTGA_FIRST_PIECE, 100), with_qlen(CTGA_SECOND_PIECE, 100)],
        )

    def test_query_walks_ahead_of_the_reference(self):
        walks = [line for line in WALKS.split('\n') if line]
        queries_first = HEADER + ''.join(f'{line}\n' for line in [*[w for w in walks if 'GRCh38' not in w], walks[0]])
        held = convert(queries_first, ['--reference', 'GRCh38'])
        self.assertEqual(held.returncode, 0, held.stderr)
        self.assertEqual(paf_rows(held.stdout), ALL_ROWS)

        # minigraph-cactus writes one chromosome's S, L and W lines after
        # another's, so a second reference contig can follow the query walks;
        # node 11 was private to ctgA when it was aligned, which the guard catches
        late_reference = WITH_WALKS + 'W\tGRCh38\t0\tchrB\t0\t4\t>11\n'
        refused = convert(late_reference, ['--reference', 'GRCh38'])
        self.assertNotEqual(refused.returncode, 0)
        self.assertRegex(refused.stderr, r'GRCh38#0 chrB:0 arrived after a query walk that visits 1 of its nodes')

        held_all = convert(late_reference, ['--reference', 'GRCh38', '--hold-queries'])
        self.assertEqual(held_all.returncode, 0, held_all.stderr)
        self.assertEqual(
            [r for r in paf_rows(held_all.stdout) if r[0] == 'HG01109#1#ctgA'],
            [
                row('HG01109#1#ctgA', 70, 0, 24, '+', 'GRCh38#0#chrA', 59, 0, 27, 23, 27, 255, 'cg:Z:10=1X5=3D8='),
                row('HG01109#1#ctgA', 70, 24, 28, '+', 'GRCh38#0#chrB', 4, 0, 4, 4, 4, 255, 'cg:Z:4='),
                row('HG01109#1#ctgA', 70, 28, 50, '+', 'GRCh38#0#chrA', 59, 27, 49, 22, 22, 255, 'cg:Z:22='),
                CTGA_SECOND_PIECE,
            ],
        )

    def test_wrong_reference_fails(self):
        run = convert(WITH_WALKS, ['--reference', 'CHM13'])
        self.assertNotEqual(run.returncode, 0)
        self.assertRegex(run.stderr, r'no CHM13#0 walk in the input')


if __name__ == '__main__':
    unittest.main()
