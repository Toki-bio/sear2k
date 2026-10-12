# sear2k - short sequence search in genomes

`sear` scans a genome FASTA with short query sequences using `ssearch36`, splitting the genome into overlapping fragments to avoid boundary misses. It merges fragment hits back into genome coordinates, extracts sequence with flanks, and can optionally keep the top-N best hits and align them with MAFFT.

For **queries ≥ 1000 bp**, sear automatically switches to `minimap2`, with the query as the index and the genome streamed through as overlapping windows. This is what makes it work for repeats such as LINEs (see Minimap2 mode).

## Highlights

- Fragmented scanning (2500 bp windows, 2000 bp step) for short queries
- **Automatic minimap2 mode** for queries ≥ 1000 bp: query indexed, genome as overlapping windows, base-level alignments (`-c`)
- Strand-aware merging and coordinate reconstruction
- Optional best-hit selection + MAFFT alignment with auto-sizing for long queries
- Safe handling of long headers and cleanup of intermediates
- Optional slop control for output BED/FASTA

## Requirements

Required:
- `ssearch36` (FASTA package) — for queries < 1000 bp
- `minimap2` — for queries ≥ 1000 bp (auto-detected)
- `bedtools`
- `seqkit` — genome windows (both modes)
- `samtools`

Optional (best-hit mode):
- `mafft`

Note: `ssearch36` is only required in ssearch36 mode (queries < 1000 bp).

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
- `--mode regular|fast`, `--fast` (or `SEAR_MODE=fast`): fast mode runs round 1 with `fasta36` (ktup 6; HSPs within 50 bp count as one locus for the length rule) and then one `ssearch36 -z 1` round on the minus bank, the fragments with the round-1 loci cut out, so what fasta36 misses ssearch36 still finds. Same loci within 0.1 %, about 2x faster (more on hit-rich families). Needs `fasta36` in PATH.

Notes:
- `--slop` overrides the positional `Flank` value.
- If your query FASTA does not end with a newline, `sear` inserts one before MAFFT input to avoid header leakage.

## Minimap2 mode

When the query is **≥ 1000 bp**, sear uses `minimap2` instead of ssearch36. The output format (BED7) is the same.

**How:** the roles are swapped. The query is minimap2's reference and the genome is streamed through `seqkit sliding` as overlapping windows (the "reads"): step `max(query length, 2 kb)`, width `step + query length + 1 kb`, so every copy up to query length + 1 kb lies whole inside at least one window. Alignments are base-level (`-x map-ont -k 11 -w 3 -c`); hits from overlapping windows are merged, then the usual identity and length filters apply. Nothing is written to disk but the hits.

**Why:** indexing the genome instead (what sear did before 2026-10) loses almost every copy of an interspersed repeat, because its minimizers occur hundreds of times and minimap2 drops high-occurrence seeds. Test on human chr21 with L1.3 (6 kb) against RepeatMasker: the old mode found 4 of 73 full-length L1s, the swapped mode 65 at any length; with sear's default length rule (hit ≥ 0.8 × query) it keeps 52, as many as blastn dc-megablast keeps under the same rule (50). The old mode also ran without `-c`, so identity read ~70 % for 98 % copies and the score column was always 0.

**Columns:** identity = matching bases / alignment block length (PAF `$10/$11`); score = minimap2 `AS:i`.

**Settings:** `map-ont` with `-k 11 -w 3` by default. Set `SEAR_MM2_PRESET` (e.g. `asm20`) or `SEAR_MM2_OPTS` (e.g. `"-k 15 -w 10"`) to change them. The denser seeds raised recall on a 1.5 kb L1 3′ end from 77 % to 89 % of dc-megablast's loci, for about twice the time (a 6 kb query takes ~9 s on human chr21 and ~2 min on a 1.65 Gbp genome with 4 threads).

**`-b` auto-sizing (minimap2 mode only):** When `-b` is used without an explicit N:
- Query < 5000 bp → `-b 10` (10 best hits + MAFFT alignment)
- Query ≥ 5000 bp → `-b 5` (5 best hits + MAFFT alignment)

An explicit value (e.g. `-b 20`) always overrides auto-sizing. Best hits are ranked by `AS`.

**Override:** `--force-ssearch` uses ssearch36 for any query length. ssearch36 scans 2.5 kb windows, so it can only report hits up to 2500 bp: if query length × Slf exceeds that (above 3125 bp at Slf 0.8), sear warns and skips the search instead of returning an empty result after a full scan.

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
