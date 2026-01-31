# sear2k — short sequence search in genomes (ssearch36 + fragment scanning)

**License:** MIT

`sear2k` searches for short query sequences (e.g., SINEs, small repeats, motifs) across genome FASTA files using `ssearch36` (FASTA package). To avoid missing matches on long contigs, it **splits the genome into overlapping fragments**, searches each fragment, and then converts fragment hits back into genome coordinates.

The pipeline is aimed at fast genome-wide screening for short elements where a single-pass search against full-length contigs may be inefficient or miss boundary-spanning alignments.

---

## How it works

1. **Index genome** (`samtools faidx`) if needed
2. **Create overlapping genome fragments** (default: 2500 bp windows with 2000 bp step → 500 bp overlap)
3. **Run `ssearch36`** for the query against each fragment bank part
4. **Filter hits** by minimum aligned length and percent identity
5. **(Optional) Multi-round depletion search**: remove hit intervals from the fragment bank and re-search remaining sequence
6. **Convert fragment coordinates back to genome coordinates**
7. **Merge overlapping hits** by strand
8. **Extract hit sequences** (FASTA) with flanks
9. **(Optional) Best-hit mode**: pick top N hits by bitscore and align them with MAFFT

---

## Features

* Fragment-based scanning with overlap to reduce missed hits
* Adjustable thresholds:

  * minimum aligned length as a fraction of query length (default `0.9`)
  * minimum percent identity (default `65`)
* Multi-round depletion search (default) for more comprehensive discovery
  *(can be disabled with `-m/--minus` for a faster single-round run)*
* Output:

  * merged **BED** with hit stats (includes **bitscore**)
  * extracted **FASTA** sequences with flanks
* Optional **top-N best hits** output and **MAFFT alignment**
* Handles underscores in FASTA sequence names by using a temporary header-safe copy and restoring names in outputs

---

## Requirements

### Required

* `ssearch36` (from the FASTA package)
* `bedtools`
* `seqkit`
* `samtools`

### Optional (only for `-b/--best`)

* `mafft`
* `nproc` (typically from coreutils)

---

## Installation

Make the script executable and ensure dependencies are on your `PATH`:

```bash
chmod +x sear2k
```

---

## Usage

```bash
sear2k [-b|--best [N]] [-m|--minus] <query.fasta> <genome.fasta> \
  [sequence_length_fraction=0.9] [similarity_threshold=65] [flank_size=50] [keep_splits=0]
```

### Options

* `-b`, `--best [N]`
  Select top **N** hits by **bitscore** (default **100**), extract them with **±50 bp flanks**, and align them to the query using MAFFT.

* `-m`, `--minus`
  Skip the minus-bank depletion steps (i.e., **only one round per split bank**). This is faster but may be less thorough.

### Positional arguments

* `query.fasta`
  Query sequence(s).

* `genome.fasta`
  Genome FASTA to search.

* `sequence_length_fraction` *(default: `0.9`)*
  Minimum aligned length threshold relative to query length.

* `similarity_threshold` *(default: `65`)*
  Minimum percent identity.

* `flank_size` *(default: `50`)*
  Flank size added to both sides of extracted sequences in the main `.bnk` output.

* `keep_splits` *(default: `0`)*
  `0` deletes temporary split files; `1` keeps them for reuse.

---

## Examples

### Basic run (default thresholds)

```bash
sear2k query.fa genome.fa
```

### Change identity threshold and flanks

```bash
sear2k query.fa genome.fa 0.9 60 200 0
```

### Keep split files for reuse

```bash
sear2k query.fa genome.fa 0.9 65 50 1
```

### Fast single-round mode (no depletion)

```bash
sear2k --minus query.fa genome.fa
# or
sear2k -m query.fa genome.fa
```

### Best-hit mode (top 100 by default) + MAFFT

```bash
sear2k -b query.fa genome.fa
```

### Best-hit mode with explicit N

```bash
sear2k --best 200 query.fa genome.fa
```

### Combine best-hit mode + single-round mode

```bash
sear2k -m -b 200 query.fa genome.fa
```

---

## Output

Let:

* `qname` = query filename without extension
* `taxname` = first 3 characters of the working genome filename (as defined in the script)

### Main outputs

* `taxname-qname.bed`
  Merged hit intervals (strand-aware). The script reports hit statistics per merged interval and includes **bitscore**.

* `taxname-qname.bnk`
  FASTA sequences extracted from the genome for the merged hits with **±flank_size** bases.

### Best mode outputs (`-b/--best`)

* `taxname-qname.bestN.bed`
  BED of the top-N hits (score = bitscore).

* `taxname-qname.bestN.flank50.bed`
  Top-N hits expanded by **±50 bp** for extraction/alignment.

* `taxname-qname.bestN.fa`
  FASTA sequences for best hits (with ±50 bp flanks).

* `taxname-qname.bestN.mafft.fa`
  MAFFT alignment of `query + best-hit sequences`.

---

## Notes / caveats

* **Fragment parameters**: default windows/overlap are tuned for short elements. For longer queries, increase window size and overlap in the script.
* **High-copy elements**: runtime and intermediate file size can grow quickly.
* **Heuristic variability**: `ssearch36` is heuristic; borderline hits can vary across runs.
* **Underscores in FASTA names**: the script uses a temporary header-safe copy and restores names in output files.

---

## License

MIT
