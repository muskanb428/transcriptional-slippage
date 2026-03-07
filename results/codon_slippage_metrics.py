#!/usr/bin/env python3
"""
codon_slippage_metrics.py

Compute codon heterogeneity and DNA homopolymer suppression metrics for
amino-acid homopolymer regions encoded by CDS annotations.

Fix included:
  Filters out non-protein-coding annotations by restricting to transcripts
  (and/or genes) annotated as protein_coding in the GTF.

Outputs:
  <species>.runs.tsv
  <species>.aa_summary.tsv
  <species>.species_summary.tsv
"""

from __future__ import annotations

import argparse
import collections
import csv
import dataclasses
import math
import os
import random
import re
import statistics
import sys
from typing import Dict, List, Tuple, Optional, Set


# -----------------------------
# Genetic code (standard)
# -----------------------------
CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_TO_CODONS: Dict[str, List[str]] = collections.defaultdict(list)
for cod, aa in CODON_TABLE.items():
    if aa != "*":
        AA_TO_CODONS[aa].append(cod)

DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


# -----------------------------
# FASTA access
# -----------------------------
class GenomeFASTA:
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self._use_pyfaidx = False
        self._fa = None
        self._seqs: Optional[Dict[str, str]] = None

        try:
            import pyfaidx  # type: ignore
            self._fa = pyfaidx.Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
            self._use_pyfaidx = True
        except Exception:
            self._seqs = self._load_into_memory(fasta_path)

    def _load_into_memory(self, fasta_path: str) -> Dict[str, str]:
        seqs: Dict[str, List[str]] = {}
        name = None
        with open(fasta_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    if name in seqs:
                        raise ValueError(f"Duplicate FASTA contig name: {name}")
                    seqs[name] = []
                else:
                    if name is None:
                        raise ValueError("FASTA malformed: sequence before header")
                    seqs[name].append(line.upper())
        return {k: "".join(v) for k, v in seqs.items()}

    def fetch(self, chrom: str, start_1based: int, end_1based: int) -> str:
        """1-based inclusive coordinates."""
        if start_1based > end_1based:
            return ""
        if self._use_pyfaidx:
            return str(self._fa[chrom][start_1based - 1:end_1based])
        assert self._seqs is not None
        seq = self._seqs.get(chrom)
        if seq is None:
            raise KeyError(f"Chrom not found in FASTA: {chrom}")
        return seq[start_1based - 1:end_1based]


# -----------------------------
# GTF parsing with protein-coding filtering
# -----------------------------
def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_field):
        attrs[m.group(1)] = m.group(2)
    return attrs


def get_biotype(attrs: Dict[str, str], level: str) -> Optional[str]:
    """
    level: 'transcript' or 'gene'
    Returns a biotype/type if found, else None.
    """
    keys = []
    if level == "transcript":
        keys = ["transcript_biotype", "transcript_type", "biotype"]
    elif level == "gene":
        keys = ["gene_biotype", "gene_type", "biotype"]
    else:
        keys = ["biotype"]
    for k in keys:
        if k in attrs and attrs[k]:
            return attrs[k]
    return None


def is_protein_coding_biotype(b: Optional[str]) -> bool:
    if b is None:
        return False
    return b.lower() == "protein_coding"


@dataclasses.dataclass
class CDSExon:
    chrom: str
    start: int
    end: int
    strand: str
    frame: Optional[int]  # 0/1/2 or None


@dataclasses.dataclass
class TranscriptCDS:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: List[CDSExon]


def load_protein_coding_ids_from_gtf(gtf_path: str) -> Tuple[Set[str], Set[str]]:
    """
    First pass over GTF to collect protein_coding transcript IDs and gene IDs
    from 'transcript' and 'gene' feature lines.

    Returns: (protein_coding_transcripts, protein_coding_genes)
    """
    pc_tx: Set[str] = set()
    pc_genes: Set[str] = set()

    with open(gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs_s = parts
            if feature not in ("transcript", "gene"):
                continue
            attrs = parse_gtf_attributes(attrs_s)
            gid = attrs.get("gene_id")
            tid = attrs.get("transcript_id")

            if feature == "gene":
                b = get_biotype(attrs, "gene")
                if gid and is_protein_coding_biotype(b):
                    pc_genes.add(gid)
            elif feature == "transcript":
                b = get_biotype(attrs, "transcript")
                # Some GTFs only have gene biotype on transcript lines
                if b is None:
                    b = get_biotype(attrs, "gene")
                if tid and is_protein_coding_biotype(b):
                    pc_tx.add(tid)
                if gid and is_protein_coding_biotype(get_biotype(attrs, "gene")):
                    pc_genes.add(gid)

    return pc_tx, pc_genes


def load_cds_from_gtf_protein_coding_only(gtf_path: str) -> Dict[str, TranscriptCDS]:
    """
    Two-pass load:
      - Pass 1: identify protein_coding transcript/gene IDs from transcript/gene lines
      - Pass 2: load CDS lines only if transcript_id (or gene_id) is protein_coding
               If GTF lacks transcript/gene lines, fall back to CDS biotype if present.
    """
    pc_tx, pc_genes = load_protein_coding_ids_from_gtf(gtf_path)
    have_pc_lists = (len(pc_tx) > 0) or (len(pc_genes) > 0)

    tx: Dict[str, TranscriptCDS] = {}
    skipped_noncoding = 0
    kept = 0
    cds_with_biotype = 0

    with open(gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs_s = parts
            if feature != "CDS":
                continue

            attrs = parse_gtf_attributes(attrs_s)
            tid = attrs.get("transcript_id")
            gid = attrs.get("gene_id", "NA")
            if tid is None:
                continue

            # Determine if this CDS belongs to a protein_coding transcript/gene
            keep = True
            if have_pc_lists:
                keep = (tid in pc_tx) or (gid in pc_genes)
            else:
                # Fallback: use biotype fields on CDS line if transcript/gene features missing
                b_t = get_biotype(attrs, "transcript")
                b_g = get_biotype(attrs, "gene")
                if b_t is not None or b_g is not None:
                    cds_with_biotype += 1
                    keep = is_protein_coding_biotype(b_t) or is_protein_coding_biotype(b_g)
                else:
                    # No way to tell; keep (but warn later)
                    keep = True

            if not keep:
                skipped_noncoding += 1
                continue

            frame_i: Optional[int] = None
            if frame in ("0", "1", "2"):
                frame_i = int(frame)

            exon = CDSExon(chrom=chrom, start=int(start), end=int(end), strand=strand, frame=frame_i)

            if tid not in tx:
                tx[tid] = TranscriptCDS(transcript_id=tid, gene_id=gid, chrom=chrom, strand=strand, exons=[])
            tx[tid].exons.append(exon)
            kept += 1

    # Sort exons in transcript order
    for t in tx.values():
        if t.strand == "+":
            t.exons.sort(key=lambda e: (e.start, e.end))
        else:
            t.exons.sort(key=lambda e: (e.start, e.end), reverse=True)

    # Informative warnings (stderr)
    if have_pc_lists:
        if kept == 0:
            print(
                "WARNING: protein_coding filtering found no CDS to keep. "
                "Check if your GTF uses different biotype keys or lacks transcript/gene lines.",
                file=sys.stderr
            )
        else:
            print(
                f"INFO: Protein-coding filter active. Kept {kept} CDS records; skipped {skipped_noncoding}.",
                file=sys.stderr
            )
    else:
        if cds_with_biotype == 0:
            print(
                "WARNING: No transcript/gene features and no biotype fields on CDS lines. "
                "Unable to filter non-protein-coding models using GTF metadata, so all CDS were kept.",
                file=sys.stderr
            )
        else:
            print(
                f"INFO: Filtered using CDS-line biotypes. Kept {kept} CDS records; skipped {skipped_noncoding}.",
                file=sys.stderr
            )

    return tx


# -----------------------------
# Core helpers
# -----------------------------
def translate_cds(cds_seq: str) -> str:
    cds_seq = cds_seq.upper()
    usable = (len(cds_seq) // 3) * 3
    cds_seq = cds_seq[:usable]
    prot: List[str] = []
    for i in range(0, len(cds_seq), 3):
        cod = cds_seq[i:i + 3]
        aa = CODON_TABLE.get(cod, "X")
        if aa == "*":
            break
        prot.append(aa)
    return "".join(prot)


def longest_homopolymer_nt(seq: str) -> int:
    if not seq:
        return 0
    best = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


def homopolymer_excess(seq: str, threshold: int = 3) -> int:
    if not seq:
        return 0
    total = 0
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
        else:
            if cur > threshold:
                total += (cur - threshold)
            cur = 1
    if cur > threshold:
        total += (cur - threshold)
    return total


def shannon_entropy_from_counts(counts: Dict[str, int]) -> float:
    tot = sum(counts.values())
    if tot == 0:
        return 0.0
    ent = 0.0
    for v in counts.values():
        p = v / tot
        ent -= p * math.log(p, 2)
    return ent


def gini_simpson_from_counts(counts: Dict[str, int]) -> float:
    tot = sum(counts.values())
    if tot == 0:
        return 0.0
    s = 0.0
    for v in counts.values():
        p = v / tot
        s += p * p
    return 1.0 - s


def find_aa_homopolymers(protein: str, min_run: int) -> List[Tuple[int, int, str]]:
    runs: List[Tuple[int, int, str]] = []
    i = 0
    while i < len(protein):
        aa = protein[i]
        j = i + 1
        while j < len(protein) and protein[j] == aa:
            j += 1
        if aa != "X" and (j - i) >= min_run:
            runs.append((i, j, aa))
        i = j
    return runs


def build_spliced_cds(genome: GenomeFASTA, t: TranscriptCDS) -> str:
    pieces = [genome.fetch(e.chrom, e.start, e.end) for e in t.exons]
    cds = "".join(pieces)
    if t.strand == "-":
        cds = revcomp(cds)

    first_frame: Optional[int] = None
    for e in t.exons:
        if e.frame is not None:
            first_frame = e.frame
            break
    if first_frame in (1, 2):
        cds = cds[first_frame:]
    return cds


def estimate_codon_usage(all_cds: List[str]) -> Dict[str, float]:
    counts = collections.Counter()
    total = 0
    for cds in all_cds:
        cds = cds.upper()
        usable = (len(cds) // 3) * 3
        cds = cds[:usable]
        for i in range(0, len(cds), 3):
            cod = cds[i:i + 3]
            aa = CODON_TABLE.get(cod)
            if aa and aa != "*":
                counts[cod] += 1
                total += 1
    if total == 0:
        return {}
    return {c: v / total for c, v in counts.items()}


def min_possible_max_homopolymer_for_aa_run(aa: str, run_len: int) -> int:
    codons = AA_TO_CODONS.get(aa, [])
    if not codons or run_len <= 0:
        return 0

    def end_run_len(cod: str) -> int:
        if cod[2] == cod[1] == cod[0]:
            return 3
        if cod[2] == cod[1]:
            return 2
        return 1

    def prefix_run_len(cod: str) -> int:
        if cod[0] == cod[1] == cod[2]:
            return 3
        if cod[0] == cod[1]:
            return 2
        return 1

    dp: Dict[Tuple[str, int], int] = {}

    for cod in codons:
        within = longest_homopolymer_nt(cod)
        key = (cod[-1], end_run_len(cod))
        dp[key] = min(dp.get(key, 10**9), within)

    for _ in range(1, run_len):
        ndp: Dict[Tuple[str, int], int] = {}
        for (last_base, last_run), gmax in dp.items():
            for cod in codons:
                within = longest_homopolymer_nt(cod)
                pref = prefix_run_len(cod)
                boundary = (last_run + pref) if cod[0] == last_base else pref
                new_gmax = max(gmax, within, boundary)
                key = (cod[-1], end_run_len(cod))
                ndp[key] = min(ndp.get(key, 10**9), new_gmax)
        dp = ndp

    return min(dp.values()) if dp else 0


def expected_random_max_homopolymer_for_aa_run(
    aa: str,
    run_len: int,
    codon_usage: Dict[str, float],
    n_mc: int,
    rng: random.Random
) -> float:
    codons = AA_TO_CODONS.get(aa, [])
    if not codons or run_len <= 0:
        return 0.0

    weights = [codon_usage.get(c, 0.0) for c in codons]
    s = sum(weights)
    if s <= 0:
        weights = [1.0] * len(codons)
        s = float(len(codons))
    weights = [w / s for w in weights]

    acc = 0.0
    for _ in range(n_mc):
        dna = "".join(rng.choices(codons, weights=weights, k=run_len))
        acc += float(longest_homopolymer_nt(dna))
    return acc / float(n_mc)


def suppression_index(obs: float, rnd: float, best: float) -> float:
    denom = rnd - best
    if denom <= 1e-9:
        return 0.0
    s = (rnd - obs) / denom
    return max(-2.0, min(2.0, s))


def run_species(
    genome_fa: str,
    gtf: str,
    species: str,
    outdir: str,
    min_run: int,
    homopoly_excess_thresh: int,
    mc: int,
    seed: int
) -> Tuple[str, str, str]:
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)

    genome = GenomeFASTA(genome_fa)

    # FIX: protein-coding only
    tx = load_cds_from_gtf_protein_coding_only(gtf)
    if not tx:
        raise RuntimeError("No protein-coding CDS entries found (after filtering). Check GTF biotype fields.")

    all_cds: List[str] = []
    cds_by_tx: Dict[str, Tuple[TranscriptCDS, str]] = {}

    for tid, t in tx.items():
        try:
            cds = build_spliced_cds(genome, t)
        except Exception:
            continue
        if len(cds) >= 3:
            all_cds.append(cds)
            cds_by_tx[tid] = (t, cds)

    codon_usage = estimate_codon_usage(all_cds)

    runs_path = os.path.join(outdir, f"{species}.runs.tsv")
    aa_path = os.path.join(outdir, f"{species}.aa_summary.tsv")
    sp_path = os.path.join(outdir, f"{species}.species_summary.tsv")

    with open(runs_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "species", "transcript_id", "gene_id", "chrom", "strand",
            "protein_len",
            "aa", "aa_start", "aa_end", "aa_len",
            "codon_entropy", "codon_gini_simpson", "n_unique_codons",
            "dna_len",
            "obs_max_nt_homopolymer", "obs_homopoly_excess",
            "best_possible_max_nt_homopolymer",
            "random_expected_max_nt_homopolymer",
            "delta_maxhp", "delta_norm",
            "suppression_index",
            "dna_A_frac", "dna_T_frac", "dna_G_frac", "dna_C_frac"
        ])

        for tid, (t, cds) in cds_by_tx.items():
            protein = translate_cds(cds)
            if not protein:
                continue

            usable = (len(cds) // 3) * 3
            codons_list = [cds[i:i + 3].upper() for i in range(0, usable, 3)]

            runs = find_aa_homopolymers(protein, min_run=min_run)
            if not runs:
                continue

            for s, e, aa in runs:
                if e > len(codons_list):
                    continue

                codon_seg = codons_list[s:e]
                dna = "".join(codon_seg).upper()

                counts = collections.Counter(codon_seg)
                ent = shannon_entropy_from_counts(counts)
                gs = gini_simpson_from_counts(counts)
                nuniq = len(counts)

                obs_maxhp = float(longest_homopolymer_nt(dna))
                obs_exc = int(homopolymer_excess(dna, threshold=homopoly_excess_thresh))

                best = float(min_possible_max_homopolymer_for_aa_run(aa, e - s))
                rnd = float(expected_random_max_homopolymer_for_aa_run(aa, e - s, codon_usage, mc, rng))

                delta_maxhp = rnd - obs_maxhp
                denom = (rnd - best)
                delta_norm = 0.0 if denom <= 1e-9 else (rnd - obs_maxhp) / denom

                sup = float(suppression_index(obs_maxhp, rnd, best))

                L = max(1, len(dna))
                a_frac = dna.count("A") / L
                t_frac = dna.count("T") / L
                g_frac = dna.count("G") / L
                c_frac = dna.count("C") / L

                w.writerow([
                    species, tid, t.gene_id, t.chrom, t.strand,
                    len(protein),
                    aa, s, e, (e - s),
                    f"{ent:.6f}", f"{gs:.6f}", nuniq,
                    len(dna),
                    int(obs_maxhp), obs_exc,
                    int(best),
                    f"{rnd:.6f}",
                    f"{delta_maxhp:.6f}", f"{delta_norm:.6f}",
                    f"{sup:.6f}",
                    f"{a_frac:.6f}", f"{t_frac:.6f}", f"{g_frac:.6f}", f"{c_frac:.6f}"
                ])

    # AA summary by scanning runs.tsv
    aa_rows: Dict[str, List[dict]] = collections.defaultdict(list)
    with open(runs_path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            aa = p[idx["aa"]]
            aa_rows[aa].append({
                "aa_len": int(p[idx["aa_len"]]),
                "obs": float(p[idx["obs_max_nt_homopolymer"]]),
                "best": float(p[idx["best_possible_max_nt_homopolymer"]]),
                "rnd": float(p[idx["random_expected_max_nt_homopolymer"]]),
                "sup": float(p[idx["suppression_index"]]),
                "ent": float(p[idx["codon_entropy"]]),
                "nuniq": float(p[idx["n_unique_codons"]]),
                "exc": float(p[idx["obs_homopoly_excess"]]),
                "dmax": float(p[idx["delta_maxhp"]]),
                "dnorm": float(p[idx["delta_norm"]]),
            })

    with open(aa_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "species", "aa", "n_runs",
            "mean_run_len", "median_run_len",
            "mean_obs_max_nt_homopolymer", "mean_best_possible", "mean_random_expected",
            "mean_delta_maxhp", "mean_delta_norm",
            "mean_suppression_index",
            "mean_codon_entropy", "mean_unique_codons",
            "mean_homopoly_excess"
        ])

        for aa in sorted(aa_rows.keys()):
            rows = aa_rows[aa]
            rl = [r["aa_len"] for r in rows]
            w.writerow([
                species, aa, len(rows),
                f"{statistics.mean(rl):.6f}", f"{statistics.median(rl):.6f}",
                f"{statistics.mean([r['obs'] for r in rows]):.6f}",
                f"{statistics.mean([r['best'] for r in rows]):.6f}",
                f"{statistics.mean([r['rnd'] for r in rows]):.6f}",
                f"{statistics.mean([r['dmax'] for r in rows]):.6f}",
                f"{statistics.mean([r['dnorm'] for r in rows]):.6f}",
                f"{statistics.mean([r['sup'] for r in rows]):.6f}",
                f"{statistics.mean([r['ent'] for r in rows]):.6f}",
                f"{statistics.mean([r['nuniq'] for r in rows]):.6f}",
                f"{statistics.mean([r['exc'] for r in rows]):.6f}",
            ])

    # Species summary
    obs_all, best_all, rnd_all, sup_all, len_all, dmax_all, dnorm_all = [], [], [], [], [], [], []
    with open(runs_path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            len_all.append(int(p[idx["aa_len"]]))
            obs_all.append(float(p[idx["obs_max_nt_homopolymer"]]))
            best_all.append(float(p[idx["best_possible_max_nt_homopolymer"]]))
            rnd_all.append(float(p[idx["random_expected_max_nt_homopolymer"]]))
            sup_all.append(float(p[idx["suppression_index"]]))
            dmax_all.append(float(p[idx["delta_maxhp"]]))
            dnorm_all.append(float(p[idx["delta_norm"]]))

    def safe_mean(x): return float(statistics.mean(x)) if x else 0.0
    def safe_median(x): return float(statistics.median(x)) if x else 0.0

    with open(sp_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "species", "n_runs",
            "mean_run_len", "median_run_len",
            "mean_obs_max_nt_homopolymer", "mean_best_possible", "mean_random_expected",
            "mean_delta_maxhp", "mean_delta_norm",
            "mean_suppression_index"
        ])
        w.writerow([
            species, len(len_all),
            f"{safe_mean(len_all):.6f}", f"{safe_median(len_all):.6f}",
            f"{safe_mean(obs_all):.6f}",
            f"{safe_mean(best_all):.6f}",
            f"{safe_mean(rnd_all):.6f}",
            f"{safe_mean(dmax_all):.6f}",
            f"{safe_mean(dnorm_all):.6f}",
            f"{safe_mean(sup_all):.6f}",
        ])

    return runs_path, aa_path, sp_path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compute codon heterogeneity and DNA homopolymer suppression metrics from genome FASTA + GTF."
    )
    ap.add_argument("--genome", required=True, help="Genome FASTA")
    ap.add_argument("--gtf", required=True, help="Annotation GTF (needs CDS with transcript_id)")
    ap.add_argument("--species", required=True, help="Species identifier")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--min_run", type=int, default=6, help="Minimum AA homopolymer length (default 6)")
    ap.add_argument("--homopoly_excess_thresh", type=int, default=3, help="DNA homopolymer excess threshold (default 3)")
    ap.add_argument("--mc", type=int, default=300, help="Monte Carlo samples for random expectation (default 300)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed")
    args = ap.parse_args()

    try:
        runs, aa, sp = run_species(
            genome_fa=args.genome,
            gtf=args.gtf,
            species=args.species,
            outdir=args.outdir,
            min_run=args.min_run,
            homopoly_excess_thresh=args.homopoly_excess_thresh,
            mc=args.mc,
            seed=args.seed
        )
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    print("Wrote:")
    print(" ", runs)
    print(" ", aa)
    print(" ", sp)


if __name__ == "__main__":
    main()