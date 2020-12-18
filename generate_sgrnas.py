import argparse
import logging
import os
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

_PACKAGEDIR = pathlib.Path(__file__).parent
_TESTDIR = _PACKAGEDIR / "testdata"
DNA_PAIRINGS = str.maketrans("atcgATCG", "tagcTAGC")
# TODO(jsh): surely this is in biopython?
DNA_BASES = ["A", "T", "C", "G"]


def revcomp(x):
    x = x.upper()
    return x.translate(DNA_PAIRINGS)[::-1]


def write_params_to_log(args):
    logging.info(
        f"""
    ************
    SOURCE FILES
    ************
    genbank: {args.genbank}
    model_param: {args.model_param}
    file_find: {args.file_find}
    file_filter: {args.file_filter}
    file_mismatch: {args.file_mismatch}

    ******************
    PARAMETER SETTINGS
    ******************
    locus_tag(s): {",".join(args.locus_tag)} ({len(args.locus_tag)})
    steps: {",".join(args.step)}
    downstream (nt): {args.downstream}
    off-target seed (nt): {args.off_target_seed}
    off-target ignore (#): {args.off_target_ignore}

    ***************
    FILTER CRITERIA
    ***************
    PAM remove: {",".join(args.pam_remove)} ({len(args.pam_remove)})
    sgRNA remove: {",".join(args.sgrna_remove)} ({len(args.sgrna_remove)})
    downstream remove: {",".join(args.downstream_remove)} ({len(args.downstream_remove)})
    GC lower(%): {args.gc_lower}
    GC upper(%): {args.gc_upper}
    offset upper (nt): {args.offset_upper}
    off-target upper (#): {args.off_target_upper}
    """
    )


def construct_perfect_sgrnas(args):
    # find all sgRNAs
    genes = list()
    for rec in SeqIO.parse(args.genbank, "gb"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    gene_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
                    gene_name = feature.qualifiers.get("gene", ["NA"])[0]
                    loc = feature.location
                    genes.append((gene_tag, gene_name, loc.start, loc.end))
    genes = pd.DataFrame(genes, columns=["locus_tag", "gene_name", "start", "end"])
    logging.info("...")
    logging.info(f"{genes.head(10)}")
    logging.info(f"{genes.shape}")
    logging.info("...")

    hits = list()
    gene_count = 1
    for rec in SeqIO.parse(args.genbank, "gb"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    gene_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
                    if gene_tag in args.locus_tag or "all" in args.locus_tag:
                        if "all" in args.locus_tag:
                            logging.info(
                                "Progress finding sgRNAs (n=%i): %i (%s)"
                                % (
                                    len(genes["start"].to_list()),
                                    gene_count,
                                    gene_tag,
                                ),
                            )
                        else:
                            logging.info(
                                "Progress finding sgRNAs (n=%i): %i (%s)"
                                % (len(args.locus_tag), gene_count, gene_tag),
                            )
                        nc_seq = str(feature.extract(rec.seq))
                        gene_name = feature.qualifiers.get("gene", ["NA"])[0]
                        loc = feature.location
                        m_i = 0
                        gene_count = gene_count + 1
                        for m in re.finditer(
                            "CC", str(nc_seq.upper())[0 : (len(nc_seq) - 20 - 1)]
                        ):
                            # Collect sgRNA information, PAM sequence and downstream sequence
                            m_i = m_i + 1
                            sgrna_pam = str(nc_seq[m.start() : (m.start() + 23)])
                            sgrna = revcomp(sgrna_pam)[0:20]
                            pam = revcomp(sgrna_pam)[20:23]
                            gc = len(
                                [
                                    i
                                    for (i, x) in enumerate(list(sgrna))
                                    if x in ["C", "G"]
                                ]
                            ) / len(sgrna)
                            if loc.strand == 1:
                                loc_genome = (
                                    int(loc.start) + m.start() + 3 + 1
                                )  # loc.start is nucleotide before gene, such that loc.end-loc.start = gene length
                                rev_genome = int(loc.start) + m.start()
                                down_seq = revcomp(
                                    str(rec.seq)[
                                        (rev_genome - args.downstream) : rev_genome
                                    ]
                                )
                            if loc.strand == -1:
                                loc_genome = int(loc.end) - m.start() - 3
                                rev_genome = int(loc.end) - m.start()
                                down_seq = str(rec.seq)[
                                    rev_genome : (rev_genome + args.downstream)
                                ]

                            off_target_n = 0
                            off_target_cds = 0
                            off_target_gene = 0
                            off_target_locus_tag = 0
                            if not args.off_target_ignore:
                                # Determine template or non-template off-target hits of sgRNA-seed (length = off_target_seed)
                                off1 = sgrna[(20 - args.off_target_seed) : 20]
                                off2 = revcomp(off1)
                                off_target_n = 0
                                off_target_gene = []
                                off_target_locus_tag = []
                                for t in re.finditer(r"(?=" + off1 + ")", str(rec.seq)):
                                    off1_pam = str(rec.seq)[
                                        (t.start() + args.off_target_seed + 1) : (
                                            t.start() + args.off_target_seed + 3
                                        )
                                    ]
                                    if (
                                        off1_pam == "GG"
                                    ):  # only consider hit, when there is a potential PAM sequence
                                        off_target_n = off_target_n + 1
                                        # Determine in what genes the matching seed sequences appear
                                        start = [
                                            i
                                            for i, x in enumerate(
                                                genes["start"].to_list()
                                            )
                                            if x < (t.start() + 1)
                                        ]  # + 1 is to correct for index difference genbank and python
                                        end = [
                                            i
                                            for i, x in enumerate(
                                                genes["end"].to_list()
                                            )
                                            if x
                                            >= (t.start() + args.off_target_seed + 1)
                                        ]  # + 1 is to correct for index difference genbank and python
                                        index = list(set(start).intersection(set(end)))
                                        if len(index) > 0:
                                            for i in range(len(index)):
                                                off_target_gene.append(
                                                    genes["gene_name"].to_list()[
                                                        index[i]
                                                    ]
                                                )
                                                off_target_locus_tag.append(
                                                    genes["locus_tag"].to_list()[
                                                        index[i]
                                                    ]
                                                )
                                # TODO(jsh): This is (maybe) where the off-target scan is happening.
                                for t in re.finditer(r"(?=" + off2 + ")", str(rec.seq)):
                                    off2_pam = str(rec.seq)[
                                        (t.start() - 3) : (t.start() - 1)
                                    ]
                                    if (
                                        off2_pam == "CC"
                                    ):  # only consider hit, when there is a potential PAM sequence (here reverse taken, since reverse complement is searched)
                                        off_target_n = off_target_n + 1
                                        # Determine in what genes the matching seed sequences appear
                                        start = [
                                            i
                                            for i, x in enumerate(
                                                genes["start"].to_list()
                                            )
                                            if x < (t.start() + 1)
                                        ]  # + 1 is to correct for index difference genbank and python
                                        end = [
                                            i
                                            for i, x in enumerate(
                                                genes["end"].to_list()
                                            )
                                            if x
                                            >= (t.start() + args.off_target_seed + 1)
                                        ]  # + 1 is to correct for index difference genbank and python
                                        index = list(set(start).intersection(set(end)))
                                        if len(index) > 0:
                                            for i in range(len(index)):
                                                off_target_gene.append(
                                                    genes["gene_name"].to_list()[
                                                        index[i]
                                                    ]
                                                )
                                                off_target_locus_tag.append(
                                                    genes["locus_tag"].to_list()[
                                                        index[i]
                                                    ]
                                                )

                                off_target_n = off_target_n - 1  # substract self-hit
                                off_target_gene.remove(gene_name)
                                off_target_locus_tag.remove(gene_tag)
                                off_target_cds = len(off_target_gene)
                                off_target_gene = "; ".join(off_target_gene)
                                off_target_locus_tag = "; ".join(off_target_locus_tag)

                            hits.append(
                                (
                                    args.genbank,
                                    gene_tag,
                                    gene_name,
                                    loc.strand,
                                    m_i,
                                    m.start() + 3 + 1,
                                    loc_genome,
                                    gc,
                                    sgrna,
                                    pam,
                                    down_seq,
                                    sgrna + pam + down_seq,
                                    off_target_n,
                                    off_target_cds,
                                    off_target_gene,
                                    off_target_locus_tag,
                                )
                            )

    hitframe = pd.DataFrame(
        hits,
        columns=[
            "genbank",
            "locus_tag",
            "gene_name",
            "strand",
            "n",
            "offset",
            "position_offset",
            "gc",
            "sgrna",
            "pam",
            "downstream",
            "sequence",
            "off_target_n",
            "off_target_cds",
            "off_target_gene",
            "off_target_locus_tag",
        ],
    )
    hitframe.to_csv(args.file_find, index=False)

    logging.info("...")
    logging.info(f"{hitframe.head(10)}")
    logging.info(f"{hitframe.shape}")


def str_compare(s1, s2):
    match = [False] * len(s2)
    for s_1 in range(len(s1)):
        i_nt = [m.start(0) for m in re.finditer("N", s1[s_1])]
        if len(i_nt) == len(s1[s_1]):
            # TODO(jsh): fix this case.
            continue  # don't compare against a string of 'N's or we'll drop everything
        short_remove = "".join(
            [char for idx, char in enumerate(s1[s_1]) if idx not in i_nt]
        )
        for s_2 in range(len(s2)):
            short_sequence = "".join(
                [char for idx, char in enumerate(s2[s_2]) if idx not in i_nt]
            )
            if short_remove == short_sequence:
                match[s_2] = True
    return [i for (i, x) in enumerate(match) if x == True]


def filter_guides(args):
    # Filter
    hitframe = pd.read_csv(args.file_find)
    hitframe = hitframe.drop(
        str_compare(args.pam_remove, hitframe["pam"].to_list()), axis=0
    ).reset_index(drop=True)
    hitframe = hitframe.drop(
        str_compare(args.sgrna_remove, hitframe["sgrna"].to_list()), axis=0
    ).reset_index(drop=True)
    hitframe = hitframe.drop(
        str_compare(args.downstream_remove, hitframe["downstream"].to_list()),
        axis=0,
    ).reset_index(drop=True)
    hitframe = hitframe.drop(
        [i for i, x in enumerate(hitframe["gc"].to_list()) if x > args.gc_upper],
        axis=0,
    ).reset_index(drop=True)
    hitframe = hitframe.drop(
        [i for i, x in enumerate(hitframe["gc"].to_list()) if x < args.gc_lower],
        axis=0,
    ).reset_index(drop=True)
    hitframe = hitframe.drop(
        [
            i
            for i, x in enumerate(hitframe["offset"].to_list())
            if x > args.offset_upper
        ],
        axis=0,
    ).reset_index(drop=True)
    hitframe = hitframe.drop(
        [
            i
            for i, x in enumerate(hitframe["off_target_n"].to_list())
            if x > args.off_target_upper
        ],
        axis=0,
    ).reset_index(drop=True)
    logging.info("...")
    logging.info(f"{hitframe.head(10)}")
    logging.info(f"{hitframe.shape}")
    hitframe.to_csv(args.file_filter, index=False)


def construct_mismatch_guides(args):
    # Create mismatch and predictions
    hitframe = pd.read_csv(args.file_filter)
    mm = []
    lookup = pd.read_csv(args.model_param)
    # TODO(jsh): we're computing a linear model ... by hand?!
    for i in range(len(hitframe["n"].to_list())):
        logging.info(
            "Progress generating mismatches (n=%i): %i"
            % (len(hitframe["n"].to_list()), i + 1),
        )
        sgrna = hitframe["sgrna"].to_list()[i]
        gb = hitframe["genbank"].to_list()[i]
        gene_tag = hitframe["locus_tag"].to_list()[i]
        gc = hitframe["gc"].to_list()[i]
        for j in range(len(sgrna)):
            wt = sgrna[j]
            sub = [x for (k, x) in enumerate(DNA_BASES) if x != wt]
            for m in range(len(sub)):
                i_loc = j
                mismatch = sub[m]
                mm_sgrna = list(sgrna)
                mm_sgrna[j] = sub[m]
                mm_sgrna = "".join(mm_sgrna)
                p0 = float(lookup.loc[lookup["feature"] == "intercept"]["weight"])
                p1 = gc * float(lookup.loc[lookup["feature"] == "gc_content"]["weight"])
                p2 = float(lookup.loc[lookup["feature"] == str(i_loc)]["weight"])
                p3 = float(lookup.loc[lookup["feature"] == wt + mismatch]["weight"])
                prediction = p0 + p1 + p2 + p3
                index = (
                    str(hitframe["n"].to_list()[i])
                    + "("
                    + str(j * len(sub) + m + 1)
                    + ")"
                )
                mm.append(
                    (
                        gb,
                        gene_tag,
                        index,
                        str(i_loc + 1),
                        wt,
                        sub[m],
                        sgrna,
                        mm_sgrna,
                        str(prediction),
                    )
                )

    mm = pd.DataFrame(
        mm,
        columns=[
            "genbank_file",
            "locus_tag",
            "n",
            "position",
            "original",
            "substitution",
            "full_sgrna",
            "mismatch_sgrna",
            "predicted efficacy",
        ],
    )
    logging.info("...")
    logging.info(f"{mm.head(10)}")
    logging.info(f"{mm.shape}")
    mm.to_csv(args.file_mismatch, index=False)


def parse_args():
    # Load arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--genbank", type=str, help="path to genbank file", required=True
    )
    parser.add_argument(
        "--locus_tag",
        type=str,
        help="locus tags (i.e. genes) to be analyzed",
        required=True,
    )
    parser.add_argument(
        "--model_param",
        type=str,
        help="path to .csv file containing model parameters",
        required=False,
        default=str(_TESTDIR / "model_param.csv"),
    )
    parser.add_argument(
        "--step",
        type=str,
        help="pipeline components to execute (comma separated)",
        required=False,
        default="find,filter,mismatch",
    )
    parser.add_argument(
        "--downstream",
        type=int,
        # TODO(jsh): what is "analysis"?
        help="number of downstream nucleotides to include in analysis",
        required=False,
        default=15,
    )
    parser.add_argument(
        "--off_target_seed",
        type=int,
        help="how many PAM-adjacent bases to use as seed-sequence for off-target hits",
        required=False,
        default=9,
    )
    parser.add_argument(
        "--off_target_ignore",
        action="store_true",
        help="bypass slow off-target detection step",
        default=False,
    )
    parser.add_argument(
        "--file_find",
        type=str,
        help="location to write/read .csv file of table with all sgRNAs",
        required=False,
        default=str(_TESTDIR / "sgRNA.csv"),
    )
    parser.add_argument(
        "--file_filter",
        type=str,
        help="location to write/read .csv file containing filtered sgRNAs",
        required=False,
        default=str(_TESTDIR / "sgRNA_filtered.csv"),
    )
    parser.add_argument(
        "--file_mismatch",
        type=str,
        help="location to write/read .csv file containing mismatch sgRNAs",
        required=False,
        default=str(_TESTDIR / "sgRNA_mismatched.csv"),
    )
    parser.add_argument(
        "--pam_remove",
        type=str,
        help="list of PAM sequences to filter",
        required=False,
        default="NNN",
    )
    parser.add_argument(
        "--downstream_remove",
        type=str,
        help="list of PAM sequences to filter",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--sgrna_remove",
        type=str,
        help="list of specific sgRNA targets to filter",
        required=False,
        default=("N" * 20),
    )
    parser.add_argument(
        "--gc_upper",
        type=float,
        help="gc cutoff filter (maximum-value)",
        required=False,
        default=1,
    )
    parser.add_argument(
        "--gc_lower",
        type=float,
        help="gc cutoff filter (minimum-value)",
        required=False,
        default=0,
    )
    parser.add_argument(
        "--offset_upper",
        type=int,
        help="maximum distance from start codon past which to filter",
        required=False,
        default=200,
    )
    parser.add_argument(
        "--off_target_upper",
        type=int,
        help="maximum number of seed nucleotides allowed for off-target",
        required=False,
        default=10,
    )
    args = parser.parse_args()
    args.sgrna_remove = args.sgrna_remove.split(",")
    args.pam_remove = args.pam_remove.split(",")
    if args.downstream_remove is None:
        args.downstream_remove = "N" * args.downstream
    args.downstream_remove = args.downstream_remove.split(",")
    args.locus_tag = args.locus_tag.split(",")
    args.step = args.step.split(",")
    for i in range(len(args.pam_remove)):
        problem = None
        if len(args.pam_remove[i]) != 3:
            problem = "was not 3 nucleotides long"
        elif args.pam_remove[i][1] == "G" or args.pam_remove[i][2] == "G":
            problem = "contained a 'G' in 2nd or 3rd position"
        if problem is not None:
            message = (
                f"args.pam_remove value {i} was {args.pam_remove[i]} "
                f"which {problem}."
            )
            raise ValueError(message)
    for i in range(len(args.sgrna_remove)):
        if len(args.sgrna_remove[i]) != 20:
            message = (
                f"args.sgrna_remove included {args.sgrna_remove[i]} "
                f"which is not 20 nt long."
            )
            raise ValueError(message)
    for i in range(len(args.downstream_remove)):
        if len(args.downstream_remove[i]) != args.downstream:
            message = (
                f"args.downstream_remove included {args.downstream_remove[i]} "
                f"which must be args.downstream ({args.downstream}) bases long."
            )
            raise ValueError(message)
    return args


def main():
    args = parse_args()
    write_params_to_log(args)
    if "find" in args.step:
        construct_perfect_sgrnas(args)
    # TODO(jsh): filter_guides precedes (and therefore ignores) mismatched guides.  Is this correct?
    if "filter" in args.step:
        filter_guides(args)
    if "mismatch" in args.step:
        construct_mismatch_guides(args)
    logging.info("****ALL STEPS COMPLETED****")


if __name__ == "__main__":
    main()
