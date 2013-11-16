#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python ../../Maker.py --config test.config

    '''
    print message

def parse_config(configfile):
    config = defaultdict(str)
    config['cpus'] = 12
    config['model']= 'Oryza'

    with open (configfile, 'r') as filefh:
        for line in filefh:
            line = line.rstrip()
            unit = line.split('\t')
            if len(unit) > 1:
                config[unit[0]]=unit[1]
    return config


def maker_exe(exe):
    message='''
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/opt/ncbi-blast/2.2.28+/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/opt/ncbi-blast/2.2.28+/bin/blastn #location of NCBI+ blastn executable
blastx=/opt/ncbi-blast/2.2.28+/bin/blastx #location of NCBI+ blastx executable
tblastx=/opt/ncbi-blast/2.2.28+/bin/tblastx #location of NCBI+ tblastx executable
formatdb=/usr/bin/formatdb #location of NCBI formatdb executable
blastall=/usr/bin/blastall #location of NCBI blastall executable
xdformat=/usr/local/bin/xdformat #location of WUBLAST xdformat executable
blasta=/usr/local/bin/blasta #location of WUBLAST blasta executable
RepeatMasker=/opt/repeat-masker/4-0-2/RepeatMasker
#/usr/local/bin/RepeatMasker #location of RepeatMasker executable
exonerate=/usr/local/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/opt/snap/2010-07-28/bin/snap # /usr/local/bin/snap #location of snap executable
gmhmme3=#/usr/local/bin/gmhmme3 #location of eukaryotic genemark executable
gmhmmp=#/usr/local/bin/gmhmmp #location of prokaryotic genemark executable
augustus=/opt/stajichlab/augustus/2.7/bin/augustus
fgenesh= #location of fgenesh executable

#-----Other Algorithms
probuild=/usr/local/bin/probuild #location of probuild executable (required for genemark)
'''
    with open (exe, 'w') as filefh:
        print >> filefh, message    

def maker_bopts(bopts):
    message='''
#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'

pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn=1e-10 #Blastn eval cutoff
bit_blastn=40 #Blastn bit cutoff
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)

pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx=1e-06 #Blastx eval cutoff
bit_blastx=30 #Blastx bit cutoff
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)

pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx=1e-10 #tBlastx eval cutoff
bit_tblastx=40 #tBlastx bit cutoff
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)

pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking

ep_score_limit=20 #Exonerate protein percent of maximal score threshold
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold
'''
    with open (bopts, 'w') as filefh:
        print >> filefh, message

def maker_opts(opts, config):
    message='''
#-----Genome (these are always required)
genome=''' + config['genome'] + ''' #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=''' + config['est'] + '''  #set of ESTs or assembled mRNA-seq in fasta format
altest=''' + config['altest'] + ''' #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=''' + config['protein'] + '''  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=''' + config['model'] + ''' #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/opt/maker/2.28/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=/opt/snap/2010-07-28/HMM/O.sativa.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=arabidopsis #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=''' + config['cpus'] + ''' #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=1000000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
'''
    with open (opts, 'w') as filefh:
        print >> filefh, message

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch
   
'''
number is the number of sequence in each subfiles, need to determine first by count the total of sequence in fasta
like we have 2400 sequences here, set number to 60 will create 40 subfiles 
'''
def splitfasta(fastafile, number, outdir):
    filelist = []
    record_iter = SeqIO.parse(open(fastafile),"fasta")
    for i, batch in enumerate(batch_iterator(record_iter, number)) :
        filename = "group_%i.fasta" % (i+1)
        filelist.append(outdir + '/' + filename)
        handle = open(outdir + '/' + filename, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
    return filelist


def jobarray(outdir):
    cmd='''
#!/bin/bash
#PBS -l mem=10gb
#PBS -l walltime=200:00:00
#PBS -t 1-10

NUM=10

module load stajichlab
module load stajichlab-perl
module load maker/2.28
module load augustus/2.7
module load repeat-masker/4-0-2
module load ncbi-blast


cd $PBS_O_WORKDIR''' + '/' + outdir + '''

#OPTS="-TMP /rhome/cjinfeng/BigData/00.RD/Annotation/HEG4/Maker/bin/HEG4_allpathlgv1/temp_dir"
OPTS=""

RUNID=$PBS_ARRAYID
if [ ! "$RUNID" ]; then
 RUNID=$$
fi

maker $OPTS group_$PBS_ARRAYID.maker_opts.ctl group_$PBS_ARRAYID.maker_bopts.ctl group_$PBS_ARRAYID.maker_exe.ctl >& group_$PBS_ARRAYID.maker.out

'''
    with open ('maker.jobarray.sh', 'w') as filefh:
        print >> filefh, cmd


'''create ctl file for each small fasta file and generate shell to run by qsub'''
def setMaker(subfile, config):
    gff3 = ''
    filefh = open ('maker.sh', 'w')
    for sub in subfile:
        dirname = os.path.dirname(os.path.abspath(sub))
        config['genome'] = os.path.abspath(sub)
        prefix  = os.path.splitext(os.path.abspath(sub))[0]
        base    = os.path.basename(prefix) 
        maker_exe(prefix + '.maker_exe.ctl') 
        maker_bopts(prefix + '.maker_bopts.ctl')
        maker_opts(prefix + '.maker_opts.ctl', config)
        cmd = '/opt/maker/2.28/bin/maker ' + prefix + '.maker_opts.ctl ' + prefix + '.maker_bopts.ctl ' + prefix + '.maker_exe.ctl > ' + prefix + '.maker.out\n'
        cmd += '/opt/maker/2.28/bin/gff3_merge -d ' + dirname + '/' + base + '.maker.output/' +  base + '_master_datastore_index.log' 
        #gff3.append(dirname + '/' + base + '.all.gff')
        gff3 = dirname + '/' + '*.all.gff'
        print >> filefh, cmd
    filefh.close()
    return gff3

def mergegff3(gff3):
    #gff3line = ' '.join(gff3)
    '''Or just use *.gff in /opt/maker/2.28/bin/gff3_merge -o genome.all.gff3 *.gff'''
    cmd = '/opt/maker/2.28/bin/gff3_merge -o genome.all.gff3 ' + gff3
    with open ('maker.gff3.sh', 'w') as filefh:
        print >> filefh, cmd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.config) > 0
    except:
        usage()
        sys.exit(2)

    '''read config file'''
    config = parse_config(args.config)
    '''create directory for maker run'''
    outdir = 'Maker.run.cut'
    os.system('mkdir ' + outdir)
    '''split fasta file into small sub files'''
    subfile = splitfasta(config['genome'], 60, outdir)
    '''create ctl files for each small sub files'''
    gff3    = setMaker(subfile,config)
    '''run maker.sh using Maker/bin/qsub-pbs.pl'''
    '''merge genome gff3'''
    mergegff3(gff3)
    '''run maker.gff3.sh using Maker/bin/qsub-pbs.pl'''
    
if __name__ == '__main__':
    main()

