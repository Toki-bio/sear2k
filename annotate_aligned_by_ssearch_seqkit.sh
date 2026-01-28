#!/usr/bin/env bash
# annotate_aligned_by_ssearch_seqkit.sh
# Usage:
#   ./annotate_aligned_by_ssearch_seqkit.sh <aligned_in.clw|.msf|.fa> <consensus_db.fasta> [out_prefix] [threads]
# Example:
#   ./annotate_aligned_by_ssearch_seqkit.sh sample.clw consensus_db.fasta sample.annot 8
#
# Requirements: seqkit, ssearch36 (FASTA), awk, xargs, python3
set -euo pipefail

ALIGN="$1"
DB="$2"
OUT_PREFIX="${3:-${ALIGN%.*}.annot}"
THREADS="${4:-$(nproc)}"

command -v seqkit >/dev/null 2>&1 || { echo "seqkit required but not found" >&2; exit 1; }
command -v ssearch36 >/dev/null 2>&1 || { echo "ssearch36 required but not found" >&2; exit 1; }
command -v xargs >/dev/null 2>&1 || { echo "xargs required but not found" >&2; exit 1; }

WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/annotate_split_XXXX")
trap 'rm -rf "$WORKDIR"' EXIT

# 1) produce aligned FASTA (preserve gaps) for final renaming step
ALIGNED_FA="${WORKDIR}/aligned_with_gaps.fasta"
echo "Converting input alignment -> aligned FASTA (preserve gaps): $ALIGNED_FA"
seqkit seq "$ALIGN" > "$ALIGNED_FA"

# 2) produce ungapped FASTA (seqkit seq -g) and split into single-sequence files
UNALIGNED_FA="${WORKDIR}/unaligned_ungapped.fasta"
SPLITDIR="${WORKDIR}/singles"
mkdir -p "$SPLITDIR"
echo "Producing ungapped FASTA and splitting into single-sequence files -> $SPLITDIR"
seqkit seq -g "$ALIGN" > "$UNALIGNED_FA"
seqkit split -s 1 -O "$SPLITDIR" "$UNALIGNED_FA" >/dev/null

# 3) per-file worker: run ssearch36 and write a .hit file with qid<TAB>subject<TAB>bitscore
PROCESS_SCRIPT="${WORKDIR}/process_query.sh"
cat > "$PROCESS_SCRIPT" <<'SH'
#!/usr/bin/env bash
f="$1"
db="$2"
qid=$(sed -n '1s/^>//p' "$f" | awk '{print $1}')
# run ssearch36 (tabular -m 8); pick best hit by highest numeric value in last column (NF)
ssearch36 -m 8 "$f" "$db" 2>/dev/null \
  | awk -F"\t" -v q="$qid" 'BEGIN{best_s="None";best=0;found=0} {score=$(NF)+0; if(found==0 || score>best){best=score;best_s=$2;found=1}} END{ if(found) print q "\t" best_s "\t" best; else print q "\tNone\t0"}' > "${f}.hit"
SH
chmod +x "$PROCESS_SCRIPT"

# 4) run all workers in parallel
echo "Running ssearch36 in parallel on $(find "$SPLITDIR" -type f -name '*.fasta' | wc -l) single-sequence files ($THREADS threads)"
find "$SPLITDIR" -type f -name "*.fasta" -print0 \
  | xargs -0 -n1 -P "$THREADS" -I {} bash -c '"$0" "{}" "'"$DB"'"' "$PROCESS_SCRIPT"

# 5) collect hits into single TSV (qid \t subject \t bitscore)
BEST_TSV="${OUT_PREFIX}.best_hits.tsv"
echo "Collecting per-file hits -> $BEST_TSV"
find "$SPLITDIR" -type f -name "*.hit" -print0 | xargs -0 cat > "$BEST_TSV" || true
# dedupe safety: keep last occurrence for each qid
awk -F"\t" '{best[$1]=$0} END{for(k in best) print best[k]}' "$BEST_TSV" > "${BEST_TSV}.uniq" && mv "${BEST_TSV}.uniq" "$BEST_TSV"

# 6) annotate aligned FASTA headers by prepending topHit and bitscore
ANNOT_FA="${OUT_PREFIX}.aligned.annot.fasta"
echo "Annotating aligned FASTA headers -> $ANNOT_FA"
python3 - "$ALIGNED_FA" "$BEST_TSV" "$ANNOT_FA" <<'PY'
import sys
if len(sys.argv)!=4:
    print("usage: script aligned_fa best_tsv out_annot_fa", file=sys.stderr); sys.exit(2)
aligned_fa, best_tsv, out_fa = sys.argv[1:]
best = {}
with open(best_tsv) as fh:
    for line in fh:
        line=line.rstrip("\n")
        if not line: continue
        cols=line.split("\t")
        q=cols[0]
        s=cols[1] if len(cols)>1 else "None"
        b=cols[2] if len(cols)>2 else "0"
        best[q]=(s,b)
def clean(x):
    return x.replace(" ", "_").replace("|","_")
with open(aligned_fa) as inf, open(out_fa,"w") as out:
    header=None; seq=[]
    for line in inf:
        if line.startswith(">"):
            if header is not None:
                out.write(header+"\n")
                out.write("".join(seq))
            raw=line[1:].rstrip("\n")
            qid=raw.split()[0]
            s,b = best.get(qid, ("None","0"))
            s_clean = clean(s)
            # prepend annotation and keep original header after pipe
            new_header = f">topHit={s_clean}|bitscore={b}|{raw}"
            header=new_header
            seq=[]
        else:
            seq.append(line)
    if header is not None:
        out.write(header+"\n")
        out.write("".join(seq))
PY

echo "Done.
Annotated aligned FASTA: $ANNOT_FA
Best-hits TSV: $BEST_TSV
Temporary workdir (auto-removed): $WORKDIR
"

exit 0
