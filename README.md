# gfa-to-pairwise-paf

`gfa_to_pairwise_paf.py` reads a pangenome graph's GFA once and writes PAF: one
record per chain of graph nodes a query haplotype shares with a reference path.
It needs no HAL, no MAF and no projection through the reference, and it reads
minigraph-cactus and pggb output alike (W or P paths, gzip accepted).

## Install

One file, Python 3.9+, standard library only:

```bash
curl -fO https://raw.githubusercontent.com/cmdcolin/gfa-to-pairwise-paf/v1.0.0/gfa_to_pairwise_paf.py
```

## Usage

```bash
# HPRC: two haplotypes against GRCh38, with a chrom.sizes per query
pigz -dc hprc-v2.1-mc-grch38.gfa.gz \
  | python3 gfa_to_pairwise_paf.py --reference GRCh38#0 \
      --queries HG01109#1,HG01123#1 --chrom-sizes-dir sizes/ > hprc.paf

# any path can be the reference, so two strains align to each other directly
python3 gfa_to_pairwise_paf.py ecoli.gfa.gz --reference Sakai#0 \
  --queries CFT073#0 > sakai_vs_cft073.paf
```

| Option                     | Effect                                                                                     |
| -------------------------- | ------------------------------------------------------------------------------------------ |
| `--reference <sample#hap>` | Required. A bare sample means `#0`.                                                        |
| `--queries <list>`         | Comma-separated `sample#hap`s; default every other one in the file.                        |
| `--max-gap <bp>`           | Private bp a chain may skip on either side between two anchors (default 10000).            |
| `--min-block <bp>`         | Drop records spanning fewer reference bp.                                                  |
| `--chrom-sizes-dir <dir>`  | Write `<sample>.<hap>.chrom.sizes` per query.                                              |
| `--contig-lengths <file>`  | A chrom.sizes or `.fai` with exact contig lengths, keyed by contig or `sample#hap#contig`. |
| `--no-x`                   | Write private runs as `I` then `D` instead of pairing them as `X`.                         |
| `--hold-queries`           | Align every query walk after the whole file is read; see [Line order](#line-order).        |

`python3 gfa_to_pairwise_paf.py --help` prints the full description.

## How records are built

The converter indexes each reference contig as node → (rank, offset,
orientation), then follows each query walk step by step. A step on a node the
reference visits is an anchor, and a chain is a run of anchors whose reference
ranks advance monotonically: increasing when the query traverses the nodes in
the reference's orientation, decreasing when it traverses them flipped. A chain
ends where the next anchor breaks monotonicity, changes relative orientation,
leaves the contig or sits past `--max-gap`.

Each chain becomes one PAF row: the query interval in forward coordinates
(strand `-` for a flipped chain), the reference interval, mapq 255, and a
`cg:Z:` CIGAR in the reference's forward direction, as minimap2 writes a `-`
row. A shared node is `<len>=`, since two paths through one node carry identical
sequence. The private bp between two anchors are `min(q,r)X`, then the remainder
as `I` (query only) or `D` (reference only).

**An `X` is not a base comparison.** It says the graph put different sequence
between the same two anchors; the converter never realigns the private runs.

A contig's length is the largest W `end` seen for it (a P path's is its walk
length). That is the true length only when the last piece reaches the contig's
end, so pass `--contig-lengths` when it matters. Rows go out once the input
ends, because that length sits in column 2 of every row.

Only the reference's walks and the requested queries' are parsed, and node
sequences are never kept, only their lengths, which keeps a 464-haplotype graph
tractable. stderr reports progress per reference walk, then per query the walks
read, anchors, chains, `=` bp and total columns.

### Paths

A W line is `sample hap contig start end walk`, and a contig may arrive as
several W lines with different starts (minigraph-cactus writes one per
unclipped stretch). A P line's name is read as PanSN `sample#hap#contig` when it
has two `#`; otherwise the whole name is both sample and contig, haplotype 0.

### Line order

S lines may follow W lines, as when minigraph-cactus writes one chromosome's S,
L and W lines after another's. A query walk arriving before any reference walk
waits until the end; one arriving after is aligned at once. That is right unless
a later reference walk visits a node the query walked as private, and the
converter stops with a message if one does. `--hold-queries` holds every query
walk until the file is read, which is always right and costs memory in
proportion to the walks held.

## Tests

```bash
python3 -m unittest -v
```

The fixture is a graph small enough to work every row out by hand: a SNP, an
indel each way, an inversion, a contig in two W pieces, and a node the reference
visits twice.

## See also

[gfa-to-pairwise-paf-rs](https://github.com/cmdcolin/gfa-to-pairwise-paf-rs) is a Rust port with the same flags and byte-identical output, about ten times faster, with prebuilt binaries.

## Origin

Extracted from [jbrowse-components](https://github.com/GMOD/jbrowse-components),
where it builds the pairwise alignments behind the HPRC multiway synteny demo.
The history of `gfa_to_pairwise_paf.py` comes with it.

## License

Apache 2.0, as jbrowse-components.
