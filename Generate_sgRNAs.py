### JORDI VAN GESTEL
### 2020.11.16 Generate sgRNAs, filter and mismatch sgRNA

# python Generate_sgRNAs.py -h
# python Generate_sgRNAs.py --genbank GCF_000009045.1_ASM904v1_genomic.gbff --locus_tag BSU_00010

import os
import sys
import re
import argparse
import pandas as pd 
from Bio import SeqIO

# global variables
downstream        = 15;
PAM_remove        =['NNN']
sgRNA_remove      =['N'*20]
downstream_remove =['N'*downstream]
GC_lower          = 0
GC_upper          = 1
offset_upper      = 1000
off_target_seed   = 9
off_target_upper  = 10
off_target_ignore = "FALSE"

s1_file   = 'sgRNA.csv';
s2_file   = 'sgRNA_filter.csv';
s3_file   = 'sgRNA_mismatch.csv';

DNA_bp = ['A','T','C','G']
DNA_PAIRINGS = str.maketrans('atcgATCG', 'tagcTAGC')

def revcomp(x):
    x = x.upper();
    return x.translate(DNA_PAIRINGS)[::-1]

def str_compare(s1,s2):
    match = [False]*len(s2)
    for s_1 in range(len(s1)):
        i_nt = [m.start(0) for m in re.finditer('N',s1[s_1])]
        if len(i_nt) != len(s1[s_1]):
            short_remove = "".join([char for idx, char in enumerate(s1[s_1]) if idx not in i_nt])
            for s_2 in range(len(s2)):
                short_sequence = "".join([char for idx, char in enumerate(s2[s_2]) if idx not in i_nt])  
                if short_remove == short_sequence:
                    match[s_2] = True
    return [i for (i,x) in enumerate(match) if x == True]

# Get gRNA from locus_tag
def get_gRNAs_tag(genbank,locus_tag,model_param,step):

    # print parameter settings to screen
    print('\nsource files:')
    print('genbank: '       + genbank)
    print('locus_tag(s): '     + ",".join(locus_tag) + " ("+str(len(locus_tag))+")")
    print('model_param: '   + model_param)
    print('file_find: '   + s1_file)
    print('file_filter: '   + s2_file)
    print('file_mismatch: '   + s3_file)
    print('\nparameter settings:')
    print('steps: '         + ",".join(step))
    print('downstream (nt): '+ str(downstream))
    print('off-target seed (nt): '+ str(off_target_seed))
    print('off-target ignore (#): ' + str(off_target_ignore))
    print('\nfilter criteria:')
    print('PAM remove: '   + ",".join(PAM_remove) + " ("+str(len(PAM_remove))+")")
    print('sgRNA remove: ' + ",".join(sgRNA_remove) + " ("+str(len(sgRNA_remove))+")")
    print('downstream remove: '+ ",".join(downstream_remove) + " ("+str(len(downstream_remove))+")")
    print('GC lower(%): '   + str(GC_lower))
    print('GC upper(%): '   + str(GC_upper))
    print('offset upper (nt): '   + str(offset_upper))
    print('off-target upper (#): ' + str(off_target_upper))

    for i in range(len(PAM_remove)):
        if len(PAM_remove[i]) != 3 or PAM_remove[i][1] == 'G' or PAM_remove[i][2] == 'G':
            step = []
            print("...\nERROR: PAM filter is wrong\n...")
    
    for i in range(len(sgRNA_remove)):
        if len(sgRNA_remove[i]) != 20:
            step = []
            print("...\nERROR: sgRNA filter is wrong\n...")
    
    for i in range(len(downstream_remove)):
        if len(downstream_remove[i]) != downstream:
            step = []
            print("...\nERROR: downstream filter is wrong\n...")
    
    # find all sgRNAs
    if 'find' in step:
        genes = list();
        for rec in SeqIO.parse(os.getcwd()+ '\\' + genbank,"gb"):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        gene_tag  = feature.qualifiers.get("locus_tag",["NA"])[0];
                        gene_name = feature.qualifiers.get("gene",["NA"])[0];
                        loc       = feature.location
                        genes.append((gene_tag,gene_name,loc.start,loc.end))
        genes = pd.DataFrame(genes,columns=['locus_tag','gene_name','start','end'])
        print("...")
        print(genes.head(10))
        print(genes.shape)
        print("...")
        
        hits = list();
        gene_count = 1;
        for rec in SeqIO.parse(os.getcwd()+ '\\' + genbank,"gb"):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        gene_tag    = feature.qualifiers.get("locus_tag",["NA"])[0];
                        if gene_tag in locus_tag or 'all' in locus_tag:
                            if 'all' in locus_tag:
                                print('Progress finding sgRNAs (n=%i): %i (%s)' % (len(genes['start'].to_list()),gene_count,gene_tag), end='\r')
                            else:
                                print('Progress finding sgRNAs (n=%i): %i (%s)' % (len(locus_tag),gene_count,gene_tag), end='\r')
                            nc_seq      = str(feature.extract(rec.seq));
                            gene_name   = feature.qualifiers.get("gene",["NA"])[0];
                            loc         = feature.location;
                            m_i         = 0;
                            gene_count  = gene_count+1;
                            for m in re.finditer('CC',str(nc_seq.upper())[0:(len(nc_seq)-20-1)]):
                                # Collect sgRNA information, PAM sequence and downstream sequence
                                m_i         = m_i + 1;
                                sgRNA_pam   = str(nc_seq[m.start():(m.start()+23)])
                                sgRNA       = revcomp(sgRNA_pam)[0:20]
                                pam         = revcomp(sgRNA_pam)[20:23]
                                GC          = len([i for (i,x) in enumerate(list(sgRNA)) if x in ['C','G']])/len(sgRNA)
                                if loc.strand == 1:
                                    loc_genome = int(loc.start)+m.start()+3+1 # loc.start is nucleotide before gene, such that loc.end-loc.start = gene length
                                    rev_genome = int(loc.start)+m.start()
                                    down_seq   = revcomp(str(rec.seq)[(rev_genome-downstream):rev_genome])
                                if loc.strand == -1:
                                    loc_genome = int(loc.end)-m.start()-3
                                    rev_genome = int(loc.end)-m.start()
                                    down_seq = str(rec.seq)[rev_genome:(rev_genome+downstream)]
                                
                                if off_target_ignore == "FALSE":
                                    # Determine template or non-template off-target hits of sgRNA-seed (length = off_target_seed)
                                    off1 = sgRNA[(20-off_target_seed):20]
                                    off2 = revcomp(off1)
                                    off_target_n = 0
                                    off_target_gene = []
                                    off_target_locus_tag = []
                                    for t in re.finditer(r'(?='+off1+')',str(rec.seq)):
                                        off1_PAM = str(rec.seq)[(t.start()+off_target_seed+1):(t.start()+off_target_seed+3)]
                                        if off1_PAM == 'GG': # only consider hit, when there is a potential PAM sequence
                                            off_target_n = off_target_n + 1
                                            # Determine in what genes the matching seed sequences appear
                                            start  = [i for i, x in enumerate(genes['start'].to_list()) if x < (t.start()+1)]   # + 1 is to correct for index difference genbank and python
                                            end    = [i for i, x in enumerate(genes['end'].to_list()) if x >= (t.start()+off_target_seed+1)]       # + 1 is to correct for index difference genbank and python
                                            index  = list(set(start).intersection(set(end)))
                                            if len(index) > 0:
                                                for i in range(len(index)):
                                                    off_target_gene.append(genes['gene_name'].to_list()[index[i]])
                                                    off_target_locus_tag.append(genes['locus_tag'].to_list()[index[i]])
                                    for t in re.finditer(r'(?='+off2+')',str(rec.seq)):
                                        off2_PAM = str(rec.seq)[(t.start()-3):(t.start()-1)]
                                        if off2_PAM == 'CC': # only consider hit, when there is a potential PAM sequence (here reverse taken, since reverse complement is searched)
                                            off_target_n = off_target_n + 1
                                            # Determine in what genes the matching seed sequences appear
                                            start  = [i for i, x in enumerate(genes['start'].to_list()) if x < (t.start()+1)]   # + 1 is to correct for index difference genbank and python
                                            end    = [i for i, x in enumerate(genes['end'].to_list()) if x >= (t.start()+off_target_seed +1)]       # + 1 is to correct for index difference genbank and python
                                            index  = list(set(start).intersection(set(end)))
                                            if len(index) > 0:
                                                for i in range(len(index)):
                                                    off_target_gene.append(genes['gene_name'].to_list()[index[i]])
                                                    off_target_locus_tag.append(genes['locus_tag'].to_list()[index[i]])
    
                                    off_target_n = off_target_n - 1 # substract self-hit
                                    off_target_gene.remove(gene_name)
                                    off_target_locus_tag.remove(gene_tag)
                                    off_target_CDS = len(off_target_gene)
                                    off_target_gene = "; ".join(off_target_gene)
                                    off_target_locus_tag = "; ".join(off_target_locus_tag)
                                    
                                if off_target_ignore == "TRUE":
                                    off_target_n = 0
                                    off_target_CDS = 0
                                    off_target_gene = 0
                                    off_target_locus_tag = 0
                                    
                                hits.append((genbank,gene_tag,gene_name,loc.strand,m_i,m.start()+3+1,loc_genome,GC,sgRNA,pam,down_seq,sgRNA+pam+down_seq,off_target_n,off_target_CDS,off_target_gene,off_target_locus_tag))
        
        hits = pd.DataFrame(hits,columns=['genbank','locus_tag','gene_name','strand','n','offset','position_offset','GC','sgRNA','PAM','downstream','sequence','off_target_n','off_target_CDS','off_target_gene','off_target_locus_tag'])
        hits.to_csv(os.getcwd() + '\\' + s1_file,index=False)

        print("...")
        print(hits.head(10))
        print(hits.shape)
    
    # Filter
    if 'filter' in step:
        hits = pd.read_csv(os.getcwd() + '\\' + s1_file) 
        hits = hits.drop(str_compare(PAM_remove,hits['PAM'].to_list()),axis=0).reset_index(drop=True)
        hits = hits.drop(str_compare(sgRNA_remove,hits['sgRNA'].to_list()),axis=0).reset_index(drop=True)
        hits = hits.drop(str_compare(downstream_remove,hits['downstream'].to_list()),axis=0).reset_index(drop=True)
        hits = hits.drop([i for i, x in enumerate(hits['GC'].to_list()) if x > GC_upper],axis=0).reset_index(drop=True)
        hits = hits.drop([i for i, x in enumerate(hits['GC'].to_list()) if x < GC_lower],axis=0).reset_index(drop=True)
        hits = hits.drop([i for i, x in enumerate(hits['offset'].to_list()) if x > offset_upper],axis=0).reset_index(drop=True)
        hits = hits.drop([i for i, x in enumerate(hits['off_target_n'].to_list()) if x > off_target_upper],axis=0).reset_index(drop=True)
        hits.to_csv(os.getcwd() + '\\' + s2_file,index=False)

        print("...")
        print(hits.head(10))
        print(hits.shape)
        
    # Create mismatch and predictions
    if 'mismatch' in step:
        hits = pd.read_csv(os.getcwd() + '\\' + s2_file) 
        mm   = [];
        lookup = pd.read_csv(os.getcwd() + '\\' + model_param)  
        for i in range(len(hits['n'].to_list())):
            print('Progress generating mismatches (n=%i): %i' % (len(hits['n'].to_list()),i+1),end='\r')
            sgRNA    = hits['sgRNA'].to_list()[i]
            gb       = hits['genbank'].to_list()[i]
            gene_tag = hits['locus_tag'].to_list()[i]
            GC       = hits['GC'].to_list()[i]
            for j in range(len(sgRNA)):
                wt  = sgRNA[j]
                sub = [x for (k,x) in enumerate(DNA_bp) if x != wt]
                for m in range(len(sub)):
                    i_loc = j
                    mismatch    = sub[m]
                    mm_sgRNA    = list(sgRNA)
                    mm_sgRNA[j] = sub[m]
                    mm_sgRNA    = ''.join(mm_sgRNA)
                    P0 = float(lookup.loc[lookup['feature'] == 'intercept']['weight'])
                    P1 = GC*float(lookup.loc[lookup['feature'] == 'GC_content']['weight'])
                    P2 = float(lookup.loc[lookup['feature'] == str(i_loc)]['weight'])
                    P3 = float(lookup.loc[lookup['feature'] == wt+mismatch]['weight'])
                    prediction = P0 + P1 + P2 + P3
                    index = str(hits['n'].to_list()[i])+'('+str(j*len(sub)+m+1)+')'
                    mm.append((gb,gene_tag,index,str(i_loc+1),wt,sub[m],sgRNA,mm_sgRNA,str(prediction)))
                    
        mm = pd.DataFrame(mm,columns=['genbank_file','locus_tag','n','position','original','substitution','full_sgRNA','mismatch_sgRNA','predicted efficacy'])
        mm.to_csv(os.getcwd() + '\\' + s3_file,index=False)

        print("...")
        print(mm.head(10))
        print(mm.shape)

def main():

    # Load arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
      '--genbank', type=str,
      help='Name GenBank file',
      required=True)
    parser.add_argument(
      '--locus_tag', type=str,
      action='append',
      help='Locus tags, i.e. genes to be analyzed',
      required=True)

    # optional parameters
    parser.add_argument(
      '--model_param', type=str,
      help='[optional] csv_file_name of model parameters',
      default=None)
    parser.add_argument(
      '--step', type=str,
      help='[optional] [find,filter,mismatch] Perform part of the analysis',
      default=None)
    parser.add_argument(
      '--downstream', type=str,
      help='[optional] number of downstream nt to include in analysis, default = 15',
      default=None)
    parser.add_argument(
      '--off_target_seed', type=str,
      help='[optional] seed-sequence for off-target hits, default = 9',
      default=None)
    parser.add_argument(
      '--off_target_ignore', type=str,
      help='[optional] if TRUE do not determine off-target hits to speed up analysis, default = FALSE',
      default=None)
    
    parser.add_argument(
      '--file_find', type=str,
      help='[optional] csv file name of table with all sgRNAs, default = sgRNA.csv',
      default=None)
    parser.add_argument(
      '--file_filter', type=str,
      help='[optional] csv file name of table with filtered sgRNAs, default = sgRNA_filter.csv',
      default=None)
    parser.add_argument(
      '--file_mismatch', type=str,
      help='[optional] csv file name of table with mismatch sgRNAs, default = sgRNA_mismatch.csv',
      default=None)
    
    # filter criteria
    parser.add_argument(
      '--PAM_remove', type=str,
      help='[filter] remove sgRNA with PAM sequence, default NNN',
      default=None)
    parser.add_argument(
      '--sgRNA_remove', type=str,
      help='[filter] remove sgRNA with sgRNA sequence, default [N]*20',
      default=None)
    parser.add_argument(
      '--downstream_remove', type=str,
      help='[filter] remove sgRNA with upstream sequence, default [N]*downstream',
      default=None)
    parser.add_argument(
      '--GC_upper', type=str,
      help='[filter] remove sgRNA with GC above x, default = 1',
      default=None)
    parser.add_argument(
      '--GC_lower', type=str,
      help='[filter] remove sgRNA with GC below x, default = 0',
      default=None)
    parser.add_argument(
      '--offset_upper', type=str,
      help='[filter] remove sgRNAs more than x nt from start codon, default = 200',
      default=None)
    parser.add_argument(
      '--off_target_upper', type=str,
      help='[filter] remove sgRNAs with off-target hits of seed-sequence, default = 10',
      default=None)
    
    args = parser.parse_args()

    # Determine genbank file and locus tag(s)
    genbank = args.genbank
    locus_tag = "".join(args.locus_tag).split(',');

    global downstream
    global PAM_remove
    global sgRNA_remove
    global downstream_remove
    global GC_upper
    global GC_lower
    global offset_upper
    global off_target_seed
    global off_target_upper
    global off_target_ignore
    global s1_file
    global s2_file 
    global s3_file

    # Check which arguments are given and load as global variables
    if args.downstream is not None:
        downstream = int(args.downstream);
        downstream_remove ='N'*downstream;
    if args.sgRNA_remove is not None:
        sgRNA_remove = "".join(args.sgRNA_remove).split(',');
    if args.PAM_remove is not None:
        PAM_remove = "".join(args.PAM_remove).split(',');
    if args.downstream_remove is not None:
        downstream_remove = "".join(args.downstream_remove).split(',');
    if args.GC_upper is not None:
        GC_upper = float(args.GC_upper);
    if args.GC_lower is not None:
        GC_lower = float(args.GC_lower);
    if args.offset_upper is not None:
        offset_upper = int(args.offset_upper);
    if args.off_target_seed is not None:
        off_target_seed = int(args.off_target_seed);
    if args.off_target_upper is not None:
        off_target_upper = int(args.off_target_upper);
    if args.off_target_ignore is not None:
        off_target_ignore = args.off_target_ignore;

    if args.file_find is not None:
        s1_file = args.file_find;
    if args.file_filter is not None:
        s2_file = args.file_filter;
    if args.file_mismatch is not None:
        s3_file = args.file_mismatch;
    
    if args.model_param is not None:
        model_param = args.model_param;
    else:
        model_param = 'model_param.csv'
    if args.step is not None:
        step = "".join(args.step).split(',');
    else:
        step = ['find','filter','mismatch']

    get_gRNAs_tag(genbank,locus_tag,model_param,step)

if __name__ == "__main__":
  sys.exit(main())




