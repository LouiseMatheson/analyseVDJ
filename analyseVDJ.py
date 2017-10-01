import pysam
import HTSeq as ht
import numpy as np
import glob, os
import re
import pandas as pd
pd.set_option('display.max_rows', 1000) #if more than 1000 annotations then increase to print complete read count tables to output
import statistics as st
import scipy as sc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import argparse
import sys
import io

parser = argparse.ArgumentParser()
parser.add_argument('-s','--strand',help = 'specify whether library is strand specific, one of: \'same\' (default), \'opposing\' or \'unstranded\'')
parser.add_argument('-a','--annotation',help = 'annotation file to use; should be a .txt file with columns for Chromosome, Start, End, Name and Strand (if applicable). Please provide complete path unless in current directory. By default will search for a .txt file containing the word \'annotation\' or \'Annotation\'. Chromosome names in annotation file must match those in BAM file, and if format is not 1, 2, 3...X, Y or chr1, chr2 etc please provide a list using -c flag')
parser.add_argument('-r' '--reduceNames',dest = 'reduceNames', help = 'reduce length of sample names: will remove lane number, extension, barcode, L00n.R, unique_V from sample names before plotting', action = 'store_true',default = False )
parser.add_argument('-k','--kappa', dest = 'kappa',help = 'specify F if this is not a VkJk-seq library, or if reports are not present')
parser.add_argument('-b','--bam',dest = 'bam',nargs='*',help = 'optionally list all bam files that should be used (if this option is specified, please list every file required and provide complete paths unless in current directory). If none are listed, all .bam files in the current directory and those listed with -d will be used')
parser.add_argument('-d','--directories',dest = 'directories',nargs = '*', help = 'Provide paths to any additional directories to search for BAM/annotation/report files. Does not apply for files listed with -a or -b: these must be in the current directory or the complete path must be specified. The current directory will always be searched, in addition to any specified here')
parser.add_argument('-p','--prefix',dest = 'prefix',help='Please specify a prefix to append to files generated, if required')
parser.add_argument('-l','-log',dest = 'log',help = 'Specify F if read counts used in heatmap should not be log transformed')
parser.add_argument('-c','--chromosomeNames', dest = 'chromosomeNames',nargs = '*', help = 'If format of chromosome names is not 1,2,3...X,Y,MT or chr1, chr2 etc (max number 22), please provide the names of each chromosome name present in the annotation or BAM file, eg seq1 seq2 seq3 etc')

args = parser.parse_args()

#get current directory and add other directories to search
current_dir = os.getcwd()
if args.directories is None:
    args.directories = [current_dir]
else:
    args.directories.append(current_dir)

#define prefix for files (or '' if none provided)
if args.prefix:
    prefix = args.prefix+'_'
else:
    prefix=''

#file to print to for .txt output
f = open(prefix+'out.txt','w')
f.write('analyseVDJ output')
f.close()

outfile = open(prefix+'out.txt','a')
print(sys.argv,file = outfile)
print('The following directories will be searched:\n',args.directories, file = outfile)
#identify all bam files in the folder
if args.bam:
    bam_files = args.bam
else:
    bam_files = list()
    for directory in args.directories:
        for file in glob.glob(directory+'/*.bam'):
            bam_files.append(file)

#find bam files with the same prefix - if everything before '_L00' is the same they are assumed to be from the same sample eg different J genes
samples = list()
for bam in bam_files:
    try:
        s = re.search('^(.*)_L00',re.sub('.*/','',bam)).group(1)#This assumes that files from the same sample will be in the same directory
        if s not in samples:
            samples.append(s)
    except:
        print('Samples do not contain _L00 so samples identified may not be truly unique',file = outfile)
        samples = bam_files
print('The following unique samples were identified:',file = outfile)
i = 1
for sample in samples:
    print('Sample',i,':',sample,file = outfile)
    i = i+1
sample_groups = {}
bam_to_sample = {}
for sample in samples:
    bam_in_sample = list()
    for bam in bam_files:
        if re.search(sample,bam) is not None:
            bamfile = re.sub('.*/','',bam)
            bam_in_sample.append(bamfile)
            bam_to_sample[bamfile] = sample
    sample_groups[sample] = bam_in_sample

for k,v in sample_groups.items():
    print(k,':',v,file = outfile)


#then read in annotation file: if not provided with -a, first search in directories provided + current directory for .txt file containing 'a/Annotation', otherwise ask user
annotation_file = list()
if args.annotation:
    annotation_file.append(args.annotation)
else:
    for directory in args.directories:
        for file in glob.glob(directory+'/*.txt'):
            a = re.search('annotation',file)
            if a is not None:
                annotation_file.append(file)
            A = re.search('Annotation',file)
            if A is not None:
                annotation_file.append(file)

if len(annotation_file) == 0:
    print('No annotation file found!\nPlease provide a .txt file containing the string \'annotation\' or \'Annotation\', or specify name of file with -a --annotation flag')
    quit()
#in case more than one file containing 'a/Annotation' identified:
elif len(annotation_file) > 1:
    print('More than one annotation file found:')
    i = 1
    for file in annotation_file:
        print(i,':',file)
        i = i+1
    annotation_file = annotation_file[int(input('Please enter the number of the file to use:'))-1]
print('Annotation file identified:\n',annotation_file[0],file = outfile)
annotation_file = annotation_file[0]

annotations = pd.read_table(annotation_file)

#check if there are identifiable columns for chromosome, start, end, strand (if applicable)
annotations.rename(columns = {'chr':'Chromosome','name':'Name','end':'End','strand':'Strand','start':'Start','Chr':'Chromomsome'},inplace = True)
for c in ['Chromosome','Start','End']:
    if c in annotations.columns:
        continue
    else:
        print('Unable to identify column for',c,'\nPlease use annotation file with columns labelled \'Chromosome\', \'Start\', \'End\'')
        quit()

if args.strand != 'unstranded':
    if 'Strand' not in annotations.columns:
        print('Unable to identify column for Strand\nPlease use annotation file with column labelled \'Strand\', or if unstranded specify...')
        quit()
#change strand so it's +/- instead of common alternatives
    annotations['Strand'] = annotations['Strand'].replace('F','+')
    annotations['Strand'] = annotations['Strand'].replace('R','-')
    annotations['Strand'] = annotations['Strand'].replace('-1','+')
    annotations['Strand'] = annotations['Strand'].replace('1','-')
else:
    annotations['Strand'] = '.'

#reverse strand designation for annotations if opposing strand-specific
if args.strand == 'opposing':
    annotations['Strand'] = annotations['Strand'].replace('+', 'R')
    annotations['Strand'] = annotations['Strand'].replace('-', 'F')
    annotations['Strand'] = annotations['Strand'].replace('F', '+')
    annotations['Strand'] = annotations['Strand'].replace('R', '-')

#work out which column is the name
if 'Name' not in annotations.columns:
    namecol = list(annotations.columns[-annotations.columns.isin(['Strand','Chromosome','Start','End'])])
    if len(namecol) == 1:
        annotations.rename(columns = {namecol[0]:'Name'},inplace = True)
        print('Using',namecol[0],'as name of annotations', file = outfile)
    elif len(namecol) == 0:
        annotations['Name'] = range(len(annotations))
        print('No annotation names found: using index as name',file = outfile)
    else:
        print('The following columns could be used as annotation names:')
        for item in namecol:
            print(item)
        correct_annotation = input('Please enter the name of the column you would like to use:')
        print('Using column specified by user as name:', correct_annotation,file = outfile)
        annotations.rename(columns = {correct_annotation:'Name'},inplace =True)

# provide all chromosome names for genomic array: default is to provide 1-22,X,Y, MT to cover mouse and human genome
if args.chromosomeNames:
    chrs = args.chromosomeNames
else:
    chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrMT","1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]

#Now convert to genomic array
if args.strand == 'unstranded':
    geneArr = ht.GenomicArray(chrs, stranded = False, typecode='O')
else:
    geneArr = ht.GenomicArray(chrs,typecode='O')
for r in range(len(annotations)):
        iv = ht.GenomicInterval(str(annotations['Chromosome'][r]),int(annotations['Start'][r]),int(annotations['End'][r]),annotations['Strand'][
r])
        geneArr[iv] = str(annotations['Name'][r])

#now can read in bam files and quantify reads for each
read_counts = pd.DataFrame(index = annotations['Name'])

for bam in bam_files:
    counts = {}
    for item in annotations['Name']:
        counts[item] = 0
    bfile = ht.BAM_Reader(bam)
    for a in bfile:
        if a.aligned:
            if a.aQual > 20:
                already_counted = list()
                for iv2, ann in geneArr[a.iv].steps():
                    if ann is not None:
                        if ann not in already_counted:
                            counts[ann] +=1
                            already_counted.append(ann) # ensures that a read is over more than one annotation, each is counted only once
    read_counts[re.sub('.*/','',bam)] = counts.values()

#calculate total counts for normalisation
total_counts = {}
for col in read_counts.columns:
    total_counts[col] = read_counts[col].sum()
print('Total read counts for each bam file:\n',total_counts, file = outfile)
#also across all bam files in a sample
total_sample_counts = {}
for k,v in sample_groups.items():
    total_sample_counts[k] = read_counts[v].sum().sum()
if len(samples) < len(bam_files): #ie if there is more than one bam file for at least one sample
    print('Total read counts for each sample:\n',total_sample_counts, file = outfile)

#extracting data from reports for QC - only applicable if this is VkJk
if args.kappa != 'F':
    print('Assuming data is VkJk-seq - searching for reports',file = outfile)
    sorting_reports = list()
    V_reports = list()
    misprime_reports = list()
    for directory in args.directories:
        for file in glob.glob(directory+'/*sorting_report.txt'):
            sorting_reports.append(file)
        for file in glob.glob(directory+'/*V_report.txt'):
            V_reports.append(file)
        for file in glob.glob(directory+'/*mispriming_correction.txt'):
            misprime_reports.append(file)
    #now extract relevant numbers and generate plots
    if len(sorting_reports) > 0:
        print('Identified sorting report(s):\n',sorting_reports,file = outfile)
        Jk_before_pipeline = pd.DataFrame(index = ['Jk1','Jk2','Jk4','Jk5'])
        for file in sorting_reports:
            sJ = re.findall('lane[0-9]_(.+)_[ACTG]+_L00',file)[0]
            Jk_before_pipeline[sJ] = re.findall('was detected in ([0-9]+) cases',open(file).read())
        #pie chart showing number of reads assigned to each Jk gene at start of pipeline
        pdf_pages = PdfPages(prefix+'Jk_pie_before_pipeline.pdf')
        for col in Jk_before_pipeline.columns:
            fig = plt.figure(figsize=(8,8),dpi=100)
            plt.pie(Jk_before_pipeline[col], labels=Jk_before_pipeline.index, colors = ['lightskyblue','palegoldenrod','darkturquoise','mediumorchid'])
            plt.title(col)
            pdf_pages.savefig(fig)
        pdf_pages.close()
    else: print('No sorting report identified',file = outfile)

    Pipeline = pd.DataFrame()
    if len(misprime_reports) > 0:
        print('Identified misprime correction report(s):\n', misprime_reports,file = outfile)
        for file in misprime_reports:
            sname = re.findall('lane[0-9]_(.*)_[CAGT]+_L', file)[0]
            Pipeline.loc[sname, 'Misprime corrected'] = \
                re.findall('mispriming in total: ([0-9]+) ', open(file).read())[0]
    else:
        print('No misprime correction report identified',file = outfile)
    if len(sorting_reports) > 0:
        for file in sorting_reports:
            sname = re.findall('lane[0-9]_(.*)_[CAGT]+_L', file)[0]
            Pipeline.loc[sname,'Total reads'] = \
                re.findall('Total number of sequences processed: ([0-9]+)',open(file).read())[0]
            Pipeline.loc[sname, 'Total with bait'] = \
                sum([int(x) for x in re.findall('detected in ([0-9]+) cases', open(file).read())])
            Pipeline.loc[sname, 'Downstream bait >20bp'] = \
                sum([int(x) for x in re.findall('.20bp or longer.: ([0-9]+) ', open(file).read())])
            Pipeline.loc[sname, 'Molecular identifier'] = \
                sum([int(x) for x in re.findall('known sequence [12] .[ACGT]+.: ([0-9]+)', open(file).read())])
    if len(V_reports) > 0:
        print('Identified V report(s):\n',V_reports,file = outfile)
        for file in V_reports:
            sname = re.findall('lane[0-9]_(.*)_[CAGT]+_L', file)[0]
            Pipeline.loc[sname, 'V read alignment'] = \
                re.findall('V alignments stored: ([0-9]+)', open(file).read())[0]
            Pipeline.loc[sname, 'Deduplicated'] = \
                re.findall('Unique V-alignments written out: ([0-9]+) ', open(file).read())[0]
    else:
        print('No V report identified',file = outfile)
    for k,v in total_sample_counts.items():
        sname = re.findall('lane[0-9]_(.*)_[CAGT]+',k)[0]
        Pipeline.loc[sname, 'Total unique over V genes'] = v
    Pipeline['Sample'] = Pipeline.index
    print('Read counts at each stage of pipeline:\n',Pipeline, file = outfile)
    Pipeline = pd.melt(Pipeline,id_vars = ['Sample'],value_name = 'Read counts',var_name = 'Pipeline stage')
    Pipeline['Read counts'] = pd.to_numeric(Pipeline['Read counts'])
    sns.factorplot(x = 'Sample',y = 'Read counts',hue = 'Pipeline stage',data = Pipeline,kind = 'bar')
    sns.set_style('white')
    sns.set_style('ticks',)
    sns.despine()
    sns.factorplot(x='Sample', y='Read counts', hue='Pipeline stage', data=Pipeline, kind='bar',alpha = 1, palette = sns.color_palette(['tomato','mediumspringgreen','mediumpurple','gold','cornflowerblue','lightpink', 'darkturquoise','orange']))
    plt.xticks(rotation = 40, fontsize = 8,ha='right')
    plt.savefig(prefix+'pipeline.pdf')


#pie charts for Jk representation after pipeline
    Jk_after_pipeline = pd.DataFrame()
    for k,v in sample_groups.items():
        for l in v:
            Jgene = re.findall('_R.(Jk[0-9]).unique',l)[0]
            Jk_after_pipeline.loc[Jgene,re.findall('lane[0-9]_(.*)_[ACTG]+$',k)[0]] = total_counts[l]
    print('Reads per Jk gene identified from baits:\n',Jk_before_pipeline,file = outfile)
    print('Unique Vk-recombined reads per Jk gene:\n',Jk_after_pipeline,file = outfile)
    pdf_pages = PdfPages(prefix+'Jk_pie_after_pipeline.pdf')
    for col in Jk_after_pipeline.columns:
        fig = plt.figure(figsize=(8, 8), dpi=100)
        plt.pie(Jk_after_pipeline[col], labels=Jk_after_pipeline.index,
                colors=['lightskyblue', 'palegoldenrod', 'darkturquoise', 'mediumorchid'])
        plt.title(col)
        pdf_pages.savefig(fig)
    pdf_pages.close()

    #IMGT reports
    IMGTfiles = list()
    for directory in args.directories:
        for file in glob.glob(directory+'/*Summary.txt'):
            IMGTfiles.append(file)
    if len(IMGTfiles) > 0:
        print('IMGT summary files identified:\n',IMGTfiles,file = outfile)
        IMGTsamples = {}
        for item in IMGTfiles:
            s = re.findall('(.*).[0-9].1_Summ',item)[0]
            if s in IMGTsamples:
                IMGTsamples[s].append(item)
            else:
                IMGTsamples[s] = [item]
        print('The following sets of samples were identified for IMGT reports:',file = outfile)
        for k,v in IMGTsamples.items():
            print(k,v,file = outfile)
        productivity = pd.DataFrame(index = ['Productive','Non-productive'])
        pdf_pages = PdfPages(prefix + 'IMGT_functionality_and_CDR3.pdf')
        for s in IMGTsamples:
            IMGTdf = pd.DataFrame()
            sName = re.sub('.1_Summary','',s)
            sName = re.sub('.*/','',sName)
            for subset in IMGTsamples[s]:
                IMGTdf = IMGTdf.append(pd.read_table(subset))
            IMGTdf = IMGTdf.loc[:,['Functionality','V-REGION identity %','V-GENE and allele','J-GENE and allele', 'CDR3-IMGT length']]
            IMGTdf = IMGTdf[IMGTdf['Functionality'] != 'No result']
            IMGTdf['Functionality'].replace(to_replace=" \(.*\)", value='', regex=True, inplace=True)
            IMGTdf = IMGTdf[IMGTdf['Functionality'] != 'unknown']
            IMGTdf = IMGTdf[IMGTdf['V-REGION identity %'] > 98]
            IMGTdf['V-GENE and allele'].replace(to_replace='\\*.*$', value='', regex=True, inplace=True)
            IMGTdf['J-GENE and allele'].replace(to_replace="\\*.*$", value='', regex=True, inplace=True)
            IMGTdf['V-GENE and allele'].replace(to_replace='Musmus IGKV', value='', regex=True, inplace=True)
            IMGTdf['J-GENE and allele'].replace(to_replace='Musmus IGKJ', value='Jk', regex=True, inplace=True)
            IMGTdf = IMGTdf[IMGTdf['J-GENE and allele'] != 'Jk3']
            IMGTdf = IMGTdf[IMGTdf['V-GENE and allele'] != '4-52']
            productivity.loc['Non-productive',sName] = IMGTdf.where(IMGTdf['Functionality']=='unproductive')['Functionality'].count()
            productivity.loc['Productive', sName] = \
                IMGTdf.where(IMGTdf['Functionality'] == 'productive')['Functionality'].count()
            fig = plt.figure(figsize=(8, 8), dpi=100)
            plt.pie(productivity[sName],labels=productivity.index,
                    colors=['darkturquoise', 'mediumorchid'])
            plt.title(sName)
            pdf_pages.savefig(fig)
            IMGTdf = IMGTdf[IMGTdf['CDR3-IMGT length'] != 'X']
            IMGTdf['CDR3-IMGT length'] = IMGTdf['CDR3-IMGT length'].apply(int)
            fig = plt.figure(figsize=(8, 8), dpi=100)
            sns.distplot(IMGTdf['CDR3-IMGT length'],kde = False, color = 'slateblue')
            sns.set_style('white')
            sns.set_style('ticks', )
            sns.despine()
            plt.xlim(0,18)
            plt.title(sName)
            pdf_pages.savefig(fig)
        pdf_pages.close()

    else:
        print('No IMGT summary files identified', file = outfile)

#now normalise the data
median_sample_count = st.median(list(total_sample_counts.values()))
norm_counts = pd.DataFrame(index = annotations['Name'])
log_counts = pd.DataFrame(index=annotations['Name'])
for c in read_counts.columns:
    cSample = bam_to_sample[c]
    norm_factor = total_sample_counts[cSample]
    norm_counts[c] = median_sample_count * read_counts[c] / norm_factor
    log_counts[c] = np.log2(norm_counts[c]+1)

# if -r flag specified, common strings that are predicted to not be part of the sample name are removed
if args.reduceNames:
    log_counts = log_counts.rename(columns = {col:re.sub('lane[0-9]_','',re.sub('_L00[0-9]_','_',re.sub('_[ACGT]+_','_',re.sub('_L00[0-9]_R.','_',re.sub('.unique_V','',re.sub('.bam$','',col)))))) for col in log_counts.columns})
    norm_counts = norm_counts.rename(columns={col: re.sub('lane[0-9]_', '', re.sub('_L00[0-9]_', '_',
                                                                                 re.sub('_[ACGT]+_', '_',
                                                                                        re.sub('_L00[0-9]_R.', '_',
                                                                                               re.sub('.unique_V', '',
                                                                                                      re.sub('.bam$',
                                                                                                             '',
                                                                                                             col))))))
                                            for col in norm_counts.columns})
    read_counts = read_counts.rename(columns={col: re.sub('lane[0-9]_', '', re.sub('_L00[0-9]_', '_',
                                                                                   re.sub('_[ACGT]+_', '_',
                                                                                          re.sub('_L00[0-9]_R.', '_',
                                                                                                 re.sub('.unique_V', '',
                                                                                                        re.sub('.bam$',
                                                                                                               '',
                                                                                                               col))))))
                                              for col in read_counts.columns})

#print read counts to output file
print('Raw read counts over annotations:\n',read_counts,file = outfile)
print('Normalised read counts over annotations:\n',norm_counts,file = outfile)
print('Log2 transformed, normalised read counts over annotations:\n',log_counts,file = outfile)

#now generate clustered heatmap
if args.log == 'F':
    counts_use = norm_counts
    colorbar = 'Read counts'
    print('Using normalised, untransformed read counts for heatmap',file = outfile)
else:
    counts_use = log_counts
    colorbar = 'Log2(read counts)'
    print('Using log2-transformed, normalised read counts for heatmap', file=outfile)
sns.set(color_codes=True,font_scale=0.7)
sns_plot = sns.clustermap(counts_use,method='complete',cmap='Spectral_r',cbar_kws = {'label':colorbar})
if len(annotations) > 100:
    plt.setp(sns_plot.ax_heatmap.get_yticklabels(),rotation = 0, fontsize = 3)
elif len(annotations) > 50:
    plt.setp(sns_plot.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
elif len(annotations) > 25:
    plt.setp(sns_plot.ax_heatmap.get_yticklabels(), rotation=0, fontsize=5)
else:
    plt.setp(sns_plot.ax_heatmap.get_yticklabels(), rotation=0, fontsize=7)
plt.setp(sns_plot.ax_heatmap.get_xticklabels(),rotation = 45,ha='right',fontsize=10)
sns_plot.savefig(prefix+"heatmap.pdf")



