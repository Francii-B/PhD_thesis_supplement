#!/usr/bin/env python3
"""
2-gz_plasmid_filter_cigar_id.py : Sequence identity is based on the reported CIGAR (computed as "=" ÷ ("=" + "X" + "I" + "D"))

Example usage:
  python 592-Phylign-Fulgor_sam-to-tsv.py alignments.sam.gz --qlen 55000 --min_cov 0.95 --min_id 0.95 --min_mapq 0 --filter

Input: headerless SAM (only alignment lines, no @SQ/@PG headers)

Output TSV fields:
  qname  genomeID  qcov  aligned_bases  qlen  mean_identity  max_mapq  n_segments  contigs
"""

import argparse, sys, re
from collections import defaultdict
import gzip

def parse_args():
    p = argparse.ArgumentParser(description="Compute query coverage per genome from headerless SAM")
    p.add_argument("infile", help="input SAM file (headerless, can be .gz)")
    p.add_argument("--qlen", type=int, required=True, help="query length (plasmid length)")
    p.add_argument("--min_cov", type=float, default=0.95)
    p.add_argument("--min_id", type=float, default=0.95)
    p.add_argument("--min_mapq", type=int, default=0)
    p.add_argument("--filter", action="store_true", help="output only refs passing thresholds")
    return p.parse_args()

# --- helpers ---

cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")

def parse_cigar(cigar):
    """
    Parse CIGAR string into query intervals and counts.
    Returns: intervals, aligned_bases, match_count, denom_bases
    - aligned_bases: query-aligned bases (M,=,X,I)
    - match_count: exact matches (=)
    - denom_bases: = + X + I + D (for identity denominator)
    """
    q_index = 0
    intervals = []
    aligned = 0
    first_q = None
    last_q = None
    match_count = 0
    denom_bases = 0

    for length, op in cigar_re.findall(cigar):
        length = int(length)
        if op in ("M","=","X","I"):  # consumes query
            if first_q is None:
                first_q = q_index
            last_q = q_index + length
            aligned += length
            q_index += length
        elif op == "S":  # soft clip, consumes query
            q_index += length
        # D consumes reference, not query
        if op == "=":
            match_count += length
            denom_bases += length
        elif op == "X":
            denom_bases += length
        elif op == "I":
            denom_bases += length
        elif op == "D":
            denom_bases += length

    if first_q is not None:
        intervals.append((first_q, last_q))
    return intervals, aligned, match_count, denom_bases

def merge_intervals(intervals):
    """Merge overlapping query intervals"""
    if not intervals:
        return [], 0
    ints = sorted(intervals, key=lambda x: x[0])
    merged = []
    cur_s, cur_e = ints[0]
    for s,e in ints[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    total = sum(e-s for s,e in merged)
    return merged, total

def main():
    args = parse_args()
    qlen = args.qlen

    data_intervals = defaultdict(list)
    data_aligned_bases = defaultdict(int)
    data_match = defaultdict(int)
    data_denom = defaultdict(int)
    data_mapq = defaultdict(int)
    data_segcount = defaultdict(int)
    data_contigs = defaultdict(set)

    opener = gzip.open if args.infile.endswith(".gz") else open
    with opener(args.infile, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue  # skip malformed
            qname = fields[0]
            flag = int(fields[1])
            rname_full = fields[2]
            mapq = int(fields[4])
            cigar = fields[5]

            # skip unmapped
            if rname_full == "*" or cigar == "*":
                continue
            if mapq < args.min_mapq:
                continue
            if flag & 0x100:  # secondary alignment
                continue

            genome_id, contig_id = rname_full.split(".",1) if "." in rname_full else (rname_full, rname_full)
            key = (qname, genome_id)

            intervals, aligned, match_count, denom_bases = parse_cigar(cigar)
            if not intervals:
                continue

            data_intervals[key].extend(intervals)
            data_aligned_bases[key] += aligned
            data_match[key] += match_count
            data_denom[key] += denom_bases
            data_contigs[key].add(contig_id)

            data_mapq[key] = max(data_mapq.get(key, 0), mapq)
            data_segcount[key] += 1

    # --- output ---
    print("\t".join(["qname","genomeID","qcov","aligned_bases","qlen","mean_identity","max_mapq","n_segments","contigs"]))
    for (qname, genome_id), intervals in data_intervals.items():
        merged, total = merge_intervals(intervals)
        aligned_bases = data_aligned_bases[(qname, genome_id)]
        match_count = data_match[(qname, genome_id)]
        denom_bases = data_denom[(qname, genome_id)]
        if denom_bases > 0:
            mean_identity = match_count / denom_bases
            mean_identity_str = f"{mean_identity:.4f}"
        else:
            mean_identity_str = "NA"
        qcov = total / qlen
        max_mapq = data_mapq[(qname, genome_id)]
        nseg = data_segcount[(qname, genome_id)]
        contigs_str = ";".join(sorted(data_contigs[(qname, genome_id)]))
        row = [qname, genome_id, f"{qcov:.4f}", str(aligned_bases), str(qlen),
               mean_identity_str, str(max_mapq), str(nseg), contigs_str]
        if args.filter:
            try:
                mid = float(mean_identity_str) if mean_identity_str!="NA" else 0.0
            except:
                mid = 0.0
            if qcov >= args.min_cov and mid >= args.min_id:
                print("\t".join(row))
        else:
            print("\t".join(row))

if __name__ == "__main__":
    main()

