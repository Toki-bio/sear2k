sear2k searches for short query sequences (e.g., SINEs, small repeats, motifs) across large genome FASTA files using ssearch36 from the FASTA package. To avoid missing matches in long contigs, it splits the genome into overlapping fragments, searches each fragment, then converts fragment hit coordinates back into genome coordinates.

The workflow is designed for high-copy elements and fragmented genomes where running a single ssearch36 against whole contigs can miss hits near internal boundaries or become unwieldy.

How it works

Index genome
Creates a .fai index (via samtools faidx) if not present.

Fragment the genome with overlap
Builds a sliding-window bank (default windows are 2500 bp with 2000 bp step → 500 bp overlap), then splits it into multiple FASTA part files for speed and manageability.

Search each part with ssearch36
Runs ssearch36 for the query against each fragment part.

Filter + merge hits
Filters hits by:

minimum aligned length (fraction of query length)

minimum percent identity
Overlapping hits are merged by coordinate and strand.

Convert fragment coordinates back to genome coordinates
Uses the fragment header encoding (*_sliding:*) to compute the original genome coordinates.

Extract sequences around hits
Uses BED coordinates to extract hit sequences from the genome, with user-defined flanks on both sides.

Optional: top hits + MAFFT alignment (-b/--best)
Selects the top-N hits by bitscore, extracts them (with flanks), and aligns them against the query with MAFFT.

Requirements
Required tools

ssearch36 (FASTA package)

bedtools

seqkit

samtools

Optional (only for -b/--best)

mafft

nproc

Usage
sear2k [-b|--best [N]] <query.fasta> <genome.fasta> [sequence_length_fraction=0.9] [similarity_threshold=65] [flank_size=50] [keep_splits=0]

Arguments

<query.fasta>
Query sequence(s). Typically one element, but multi-FASTA also works (behavior depends on ssearch36 output and your filtering).

<genome.fasta>
Genome FASTA to search.

sequence_length_fraction (default: 0.9)
Minimum aligned-length threshold relative to query length.

similarity_threshold (default: 65)
Minimum percent identity.

flank_size (default: 50)
Flanking bases added to both sides of extracted sequences in the main .bnk output.

keep_splits (default: 0)
0 removes temporary split files; 1 keeps them for reuse.

Best mode

-b / --best [N]
Extracts the top N hits by bitscore (default N=100) and runs MAFFT alignment against the query.

Examples:

# default thresholds and flanks
sear2k query.fa genome.fa

# relaxed identity, larger flanks, keep split files
sear2k query.fa genome.fa 0.9 60 200 1

# top 100 hits by bitscore + MAFFT
sear2k -b query.fa genome.fa

# top 200 hits by bitscore + MAFFT
sear2k --best 200 query.fa genome.fa
