#!/usr/bin/env python
# python Generate_sgRNAs.py -h
# python generate_sgrnas_v3.py --genbank GCF_000009045.1_ASM904v1_genomic.gbff --locus_tag BSU_00010
# python generate_sgrnas_v3.py --genbank GCF_000009045.1_ASM904v1_genomic.gbff --locus_tag all

import argparse
import logging
import pathlib
import re
import math

import pandas as pd
from Bio import SeqIO
from collections import Counter

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

_PACKAGEDIR = pathlib.Path(__file__).parent
DNA_PAIRINGS = str.maketrans("atcgATCG", "tagcTAGC")
DNA_BASES = ["A", "T", "C", "G"]

def revcomp(x):
    x = x.upper()
    return x.translate(DNA_PAIRINGS)[::-1]

def write_params_to_log(args):
    def join_or_drop(maybe_list):
        if maybe_list is None:
            return maybe_list
        elif type(maybe_list) == list:
            return ",".join(maybe_list)
        else:
            raise TypeError('join_or_drop() called on something other than None or list')
    def len_or_drop(maybe_list):
        if maybe_list is None:
            return 0
        elif type(maybe_list) == list:
            return len(maybe_list)
        else:
            raise TypeError('len_or_drop() called on something other than None or list')
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
    locus_tag(s): {join_or_drop(args.locus_tag)} ({len_or_drop(args.locus_tag)})
    steps: {join_or_drop(args.step)}
    downstream (nt): {args.downstream}
    off-target seed (nt): {args.off_target_seed}
    
    ***************
    FILTER CRITERIA
    ***************
    PAM remove: {join_or_drop(args.pam_remove)} ({len_or_drop(args.pam_remove)})
    sgRNA remove: {join_or_drop(args.sgrna_remove)} ({len_or_drop(args.sgrna_remove)})
    downstream remove: {join_or_drop(args.downstream_remove)} ({len_or_drop(args.downstream_remove)})
    GC lower(%): {args.gc_lower}
    GC upper(%): {args.gc_upper}
    offset upper (nt): {args.offset_upper}
    off-target upper (#): {args.off_target_upper}
    """
    )

def construct_perfect_sgrnas(args):
    # find all sgRNAs
	
    genes = list();
    for rec in SeqIO.parse(args.genbank,"gb"):
        upper_seq = str(rec.seq);
        lower_seq = revcomp(str(rec.seq));
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    gene_tag  = feature.qualifiers.get("locus_tag",["NA"])[0];
                    gene_name = feature.qualifiers.get("gene",["NA"])[0];
                    loc       = feature.location
                    genes.append((gene_tag,gene_name,loc.strand,loc.start,loc.end))
    genes = pd.DataFrame(genes,columns=['locus_tag','gene_name','strand','start','end'])
    logging.info(f"\n{genes.head(10)}")
    logging.info(f"\n{genes.shape}")
	
	# Table all offtarget hits
    offtarget = list();
    for m in re.finditer(r'(?=GG)',upper_seq):
        seed_sequence = upper_seq[(m.start()-1-args.off_target_seed):(m.start()-1)];
        offtarget.append(seed_sequence)
    for m in re.finditer(r'(?=GG)',lower_seq):
        seed_sequence = lower_seq[(m.start()-1-args.off_target_seed):(m.start()-1)];
        offtarget.append(seed_sequence)
    offtarget = Counter(offtarget)
			
    hits = list();
    gene_count = 1;
    for rec in SeqIO.parse(args.genbank,"gb"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    gene_tag    = feature.qualifiers.get("locus_tag",["NA"])[0];
                    if gene_tag in args.locus_tag or 'all' in args.locus_tag:
                        if "all" in args.locus_tag:
                            print('Progress finding sgRNAs (n=%i): %i (%s)' % (len(genes['start'].to_list()),gene_count,gene_tag), end='\r')
                        else:
                            print('Progress finding sgRNAs (n=%i): %i (%s)' % (len(args.locus_tag),gene_count,gene_tag), end='\r')
                        nc_seq      = str(feature.extract(rec.seq));
                        gene_name   = feature.qualifiers.get("gene",["NA"])[0];
                        loc         = feature.location;
                        m_i         = 0;
                        gene_count  = gene_count+1;
                        for m in re.finditer(r'(?=CC)',str(nc_seq.upper())[0:(len(nc_seq)-20-1)]):
                            # Collect sgRNA information, PAM sequence and downstream sequence
                            m_i         = m_i + 1;
                            sgRNA_pam   = str(nc_seq[m.start():(m.start()+23)])
                            sgRNA       = revcomp(sgRNA_pam)[0:20]
                            pam         = revcomp(sgRNA_pam)[20:23]
                            GC          = len([i for (i,x) in enumerate(list(sgRNA)) if x in ['C','G']])/len(sgRNA)
                            if loc.strand == 1:
                                loc_genome = int(loc.start)+m.start()+3+1 # loc.start is nucleotide before gene, such that loc.end-loc.start = gene length
                                rev_genome = int(loc.start)+m.start()
                                down_seq   = revcomp(str(rec.seq)[(rev_genome-args.downstream):rev_genome])
                            if loc.strand == -1:
                                loc_genome = int(loc.end)-m.start()-3
                                rev_genome = int(loc.end)-m.start()
                                down_seq = str(rec.seq)[rev_genome:(rev_genome+args.downstream)]
                            off_target_n = offtarget[sgRNA[(20-args.off_target_seed):20]]-1;
                            hits.append((args.genbank,gene_tag,gene_name,loc.strand,m_i,m.start()+3+1,loc_genome,GC,sgRNA,pam,down_seq,sgRNA+pam+down_seq,off_target_n))
        
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
            "off_target_n"
		],
	)
    hitframe.to_csv(args.file_find, index=False)

    logging.info(f"\n{hitframe.head(10)}")
    logging.info(f"\n{hitframe.shape}")

def filter_guides(args):
    _FLATTEN = str.maketrans("Nn", "..")

    def filter_patterns(frame,column,patterns):
        """Remove all rows from frame where frame[column] matches any provided patterns"""
        # bail if there are no patterns
        if patterns is None or len(patterns) == 0:
            return frame
        # First, replace all 'Nn's with '.'s, then make a big OR pattern
        regex ="|".join([p.translate(_FLATTEN) for p in patterns])
        return frame.loc[~frame[column].str.match(regex)]

    hitframe = pd.read_csv(args.file_find)
    hitframe = filter_patterns(hitframe, "pam", args.pam_remove)
    hitframe = filter_patterns(hitframe, "sgrna", args.sgrna_remove)
    hitframe = filter_patterns(hitframe, "downstream", args.downstream_remove)
    hitframe = hitframe.loc[hitframe.gc >= args.gc_lower]
    hitframe = hitframe.loc[hitframe.gc <= args.gc_upper]
    hitframe = hitframe.loc[hitframe.offset <= args.offset_upper]
    hitframe = hitframe.loc[hitframe.off_target_n <= args.off_target_upper]

    logging.info(f"\n{hitframe.head(10)}")
    logging.info(f"\n{hitframe.shape}")
    hitframe.to_csv(args.file_filter, index=False)

def construct_mismatch_guides(args):
    # Create mismatch and predictions
    hitframe = pd.read_csv(args.file_filter)
	# PREPROCESS LOOKUP TABLE
    lookup   = pd.read_csv(args.model_param)
    features = lookup['feature'].to_list();
    weights  = lookup['weight'].to_list();
    lookuptable = []
    for loc in range(20):
        for wt in range(len(DNA_BASES)):
            n_wt = DNA_BASES[wt]
            for mu in range(len(DNA_BASES)):
                n_mu = DNA_BASES[mu]
                if n_wt != n_mu:
                    P1 = float(weights[features.index(str(loc))])+float(weights[features.index(n_wt+n_mu)]);
                    lookuptable.append((loc,n_wt,n_mu,P1));
    lookuptable = pd.DataFrame(lookuptable,columns=['loc','wt','mu','P1'])
    loc = lookuptable['loc'].to_list()
    wt  = lookuptable['wt'].to_list()   
    mu  = lookuptable['mu'].to_list()
    P1  = lookuptable['P1'].to_list()
    intercept = float(weights[features.index('intercept')])
    GC_weight = float(weights[features.index('GC_content')])   
		  
    mm = [];
    for i in range(len(hitframe['n'].to_list())):
        if (i+1)%10==0:
            print('Progress generating mismatches (n=%i): %i' % (len(hitframe['n'].to_list()),i+1),end='\r')
        sgRNA    = hitframe['sgrna'].to_list()[i]
        hitlist  = hitframe['n'].to_list()[i]
        gb       = hitframe['genbank'].to_list()[i]
        gene_tag = hitframe['locus_tag'].to_list()[i]
        GC       = hitframe['gc'].to_list()[i]
        P0       = intercept + GC*GC_weight;
        count    = 0;
        for j in range(60):
            index = 12*math.floor(j/3) + 3*DNA_BASES.index(sgRNA[math.floor(j/3)]) + j%3
            count = count + 1;
            mm_sgRNA = sgRNA[:loc[index]]+mu[index]+sgRNA[(loc[index]+1):]
            mm.append((gb,gene_tag,str(hitlist)+'('+str(count)+')',str(loc[index]+1),wt[index],mu[index],sgRNA,mm_sgRNA,str(P0+P1[index])))
	
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
    
    logging.info(f"\n{mm.head(10)}")
    logging.info(f"\n{mm.shape}")
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
        default=str(_PACKAGEDIR / "model_param.csv"),
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
        help="number of downstream nucleotides to consider when filtering sgRNAs",
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
        "--file_find",
        type=str,
        help="location to write/read .csv file of table with all sgRNAs",
        required=False,
        default=str(_PACKAGEDIR / "sgRNA.csv"),
    )
    parser.add_argument(
        "--file_filter",
        type=str,
        help="location to write/read .csv file containing filtered sgRNAs",
        required=False,
        default=str(_PACKAGEDIR / "sgRNA_filtered.csv"),
    )
    parser.add_argument(
        "--file_mismatch",
        type=str,
        help="location to write/read .csv file containing mismatch sgRNAs",
        required=False,
        default=str(_PACKAGEDIR / "sgRNA_mismatched.csv"),
    )
    parser.add_argument(
        "--pam_remove",
        type=str,
        help="list of PAM sequences to filter",
        required=False,
        default=None,
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
        default=None,
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
        default=1000,
    )
    parser.add_argument(
        "--off_target_upper",
        type=int,
        help="maximum number of seed nucleotides allowed for off-target",
        required=False,
        default=10,
    )
    args = parser.parse_args()
    if args.sgrna_remove is not None:
        args.sgrna_remove = args.sgrna_remove.split(",")
        for i in range(len(args.sgrna_remove)):
            if len(args.sgrna_remove[i]) != 20:
                message = (
                    f"args.sgrna_remove included {args.sgrna_remove[i]} "
                    f"which is not 20 nt long."
                )
                raise ValueError(message)
            if args.sgrna_remove[i].upper().count('N')==len(args.sgrna_remove[i]):
                message = (
                    f"args.sgrna_remove included {args.sgrna_remove[i]} "
                    f"which filters all sgRNAs."
                )
                raise ValueError(message)
    if args.pam_remove is not None:
        args.pam_remove = args.pam_remove.split(",")
        for i in range(len(args.pam_remove)):
            problem = None
            if len(args.pam_remove[i]) != 3:
                problem = "was not 3 nucleotides long"
            elif args.pam_remove[i][1] == "G" or args.pam_remove[i][2] == "G":
                problem = "contained a 'G' in 2nd or 3rd position"
            elif args.pam_remove[i].upper().count('N')==len(args.pam_remove[i]):
                problem = "filtered out all sgRNAs"
            if problem is not None:
                message = (
                    f"args.pam_remove value {i} was {args.pam_remove[i]} "
                    f"which {problem}."
                )
                raise ValueError(message)
    if args.downstream_remove is not None:
        args.downstream_remove = args.downstream_remove.split(",")
        for i in range(len(args.downstream_remove)):
            if len(args.downstream_remove[i]) != args.downstream:
                message = ( f"args.downstream_remove included {args.downstream_remove[i]} "
                    f"which is not {args.downstream} nt long."
                )
                raise ValueError(message)
            if args.downstream_remove[i].upper().count('N')==len(args.downstream_remove[i]):
                message = (
                    f"args.downstream_remove included {args.downstream_remove[i]} "
                    f"which filters all sgRNAs."
                )
                raise ValueError(message)
    args.locus_tag = args.locus_tag.split(",")
    args.step = args.step.split(",")
    return args

def main():
    args = parse_args()
    write_params_to_log(args)
    if "find" in args.step:
        construct_perfect_sgrnas(args)
    if "filter" in args.step:
        filter_guides(args)
    if "mismatch" in args.step:
        construct_mismatch_guides(args)
    logging.info("****ALL STEPS COMPLETED****")

if __name__ == "__main__":
    main()
