# sear2k - short sequence search in genomes

`sear` scans a genome FASTA with short query sequences using `ssearch36`, splitting the genome into overlapping fragments to avoid boundary misses. It merges fragment hits back into genome coordinates, extracts sequence with flanks, and can optionally keep the top-N best hits and align them with MAFFT.

For **queries ≥ 1000 bp**, sear automatically switches to `minimap2` for faster whole-genome alignment (no fragment scanning needed). This is especially useful for LINE reconstruction or other long-element searches.

## Highlights

- Fragmented scanning (2500 bp windows, 2000 bp step) for short queries
- **Automatic minimap2 mode** for queries ≥ 1000 bp with sensitivity cascade (asm5 → asm10 → asm20)
- Strand-aware merging and coordinate reconstruction
- Optional best-hit selection + MAFFT alignment with auto-sizing for long queries
- Safe handling of long headers and cleanup of intermediates
- Optional slop control for output BED/FASTA

## Requirements

Required:
- `ssearch36` (FASTA package) — for queries < 1000 bp
- `minimap2` — for queries ≥ 1000 bp (auto-detected)
- `bedtools`
- `seqkit` — for queries < 1000 bp (fragment splitting)
- `samtools`

Optional (best-hit mode):
- `mafft`

Note: `ssearch36` and `seqkit` are only required in ssearch36 mode (queries < 1000 bp). If you only use sear with long queries and have `minimap2`, they are not needed.

## Install

Make the script executable and ensure dependencies are on your `PATH`:

```bash
chmod +x sear
```

## Usage

```bash
sear [options] <query.fasta> <bank.fasta> [Slf=0.8] [Homology=65] [Flank=50] [KeepSplits=0]
```

Options:
- `-b, --best [N]`: keep top N unique hits by bitscore and run MAFFT (default 100 for ssearch36 mode; auto-sized for minimap2 mode — see below)
- `-m, --minus`: enable minus-bank depletion (multi-round per split; ssearch36 mode only)
- `-s, --stop [N]`: stop after N total hits (default 1000)
- `-k, --keep`: keep split files for reuse
- `--slop L,R`: override flank sizes for output bed/bnk and best-hit flanks
- `--force-ssearch`: disable minimap2 auto-switch, always use ssearch36

Notes:
- `--slop` overrides the positional `Flank` value.
- If your query FASTA does not end with a newline, `sear` inserts one before MAFFT input to avoid header leakage.

## Minimap2 mode

When the query is **≥ 1000 bp**, sear automatically uses `minimap2` instead of ssearch36 fragment-scanning. This is transparent — the output format (BED7) is identical.

**Sensitivity cascade:** sear tries minimap2 presets in order of stringency:
1. `asm5` (≤ 5% divergence)
2. `asm10` (≤ 10% divergence)
3. `asm20` (≤ 20% divergence)

It stops at the first preset that returns any hits. This ensures the most exact matches are preferred while still finding divergent copies.

**`-b` auto-sizing (minimap2 mode only):** When `-b` is used without an explicit N:
- Query < 5000 bp → `-b 10` (10 best hits + MAFFT alignment)
- Query ≥ 5000 bp → `-b 5` (5 best hits + MAFFT alignment)

An explicit value (e.g. `-b 20`) always overrides auto-sizing.

**Override:** `--force-ssearch` disables minimap2 and forces ssearch36 for all query lengths.

## Examples

Basic run (short query, ssearch36):

```bash
sear query.fa genome.fa
```

Long query (auto-switches to minimap2):

```bash
sear long_element.fa genome.fa
```

Best-hit mode:

```bash
sear --best 200 query.fa genome.fa
```

Best-hit mode with auto-sizing (long query):

```bash
sear --best long_element.fa genome.fa
# auto-sets -b 10 for queries <5 kb, -b 5 for ≥5 kb
```

Force ssearch36 for a long query:

```bash
sear --force-ssearch long_element.fa genome.fa
```

Multi-round search with minus bank:

```bash
sear --minus query.fa genome.fa
```

Custom slop for outputs:

```bash
sear --slop 10,20 query.fa genome.fa
```

## Outputs

Let:
- `qname` = query filename without extension
- `taxname` = first 3 characters of the working genome filename

Main outputs:
- `taxname-qname.bed`: merged hit intervals (BED7: chrom, start, end, homology%, length, strand, bitscore)
- `taxname-qname.bnk`: extracted FASTA with flanks (or `--slop` values)

Best-hit mode outputs:
- `taxname-qname.bestN.mafft.fa`: MAFFT alignment of query + top-N hits

Intermediates are cleaned automatically unless `--keep` is used for split reuse.

## BED7 output format

Both ssearch36 and minimap2 modes produce the same 7-column BED:

```
chrom  start  end  homology%  length  strand  bitscore
```

In minimap2 mode, the bitscore column uses the minimap2 `AS` (alignment score) tag.

## License

MIT
