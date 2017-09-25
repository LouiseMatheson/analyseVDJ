import pysam
import HTSeq as ht
import numpy as np
import glob, os
import re
import pandas as pd
import statistics as st
import scipy as sc
import matplotlib as mpl
import math

stranded = True
#identify all bam files in the folder
bam_files = list()
for file in glob.glob('*.bam'):
    bam_files.append(file)

#find bam files with the same prefix - if everything before '_L00' is the same
samples = list()
for bam in bam_files:
    try:
        s = re.search('^(.*)_L00',bam).group(1)#would it be better to use _L00 for this? (in case .J in actual name!)
        if s not in samples:
            samples.append(s)
    #but this will fail if there is no .J
    except:
        print('Samples do not contain _L00 so samples identified may not be truly unique')
        samples = bam_files
print('The following unique samples were identified:')
i = 1
for sample in samples:
    print('Sample',i,':',sample)
    i = i+1

sample_groups = {}
bam_to_sample = {}
for sample in samples:
    bam_in_sample = list()
    for bam in bam_files:
        if re.search(sample,bam) is not None:
            bam_in_sample.append(bam)
            bam_to_sample[bam] = sample
    sample_groups[sample] = bam_in_sample



for k in sample_groups.keys():
    print(k)
for v in sample_groups.values():
    print(v)

#then read in annotation file:
annotation_file = list()
for file in glob.glob('*.txt'):
    a = re.search('annotation',file)
    if a is not None:
        annotation_file.append(file)
    A = re.search('Annotation',file)
    if A is not None:
        annotation_file.append(file)

if len(annotation_file) == 0:
    print('No annotation file found!\nPlease provide a .txt file containing the string \'annotation\' or \'Annotation\'')
    quit()
elif len(annotation_file) > 1:
    print('More than one annotation file found:')
    i = 1
    for file in annotation_file:
        print(i,':',file)
        i = i+1
    annotation_file = annotation_file[int(input('Please enter the number of the file to use:'))-1]
else:
    print('Annotation file identified:\n',annotation_file[0])
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

if stranded:
    if 'Strand' not in annotations.columns:
        print('Unable to identify column for Strand\nPlease use annotation file with column labelled \'Strand\', or if unstranded specify...')
        quit()
#change strand so it's +/- instead of common alternatives
    annotations['Strand'] == str(annotations['Strand'])
    annotations['Strand'] = annotations['Strand'].replace('F','+')
    annotations['Strand'] = annotations['Strand'].replace('R','-')
    annotations['Strand'] = annotations['Strand'].replace('-1','+')
    annotations['Strand'] = annotations['Strand'].replace('1','-')

#work out which column is the name
if 'Name' not in annotations.columns:
    namecol = list(annotations.columns[-annotations.columns.isin(['Strand','Chromosome','Start','End'])])
    if len(namecol) == 1:
        annotations.rename(columns = {namecol[0]:'Name'},inplace = True)
        print('Using',namecol[0],'as name of annotations')
    elif len(namecol) == 0:
        annotations['Name'] = range(len(annotations))
        print('No annotation names found: using index as name')
    else:
        print('The following columns could be used as annotation names:')
        for item in namecol:
            print(item)
        annotations.rename(columns = {input('Please enter the name of the column you would like to use:'):'Name'},inplace =True)
#Now convert to genomic array
#can define human too and give this as an option?
mmChr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X"
,"Y","MT"]
geneArr = ht.GenomicArray(mmChr,typecode='O')
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
                            already_counted.append(ann)
    read_counts[bam] = counts.values()

print(read_counts)

#calculate total counts for normalisation
total_counts = {}
for col in read_counts.columns:
    total_counts[col] = read_counts[col].sum()
print(total_counts)
total_sample_counts = {}
for k,v in sample_groups.items():
    total_sample_counts[k] = read_counts[v].sum().sum()
#print(total_sample_counts)

#now normalise the data
median_sample_count = st.median(list(total_sample_counts.values()))
print(median_sample_count)
norm_counts = pd.DataFrame(index = annotations['Name'])
log_counts = pd.DataFrame(index=annotations['Name'])
for c in read_counts.columns:
    cSample = bam_to_sample[c]
    print(cSample)
    norm_factor = total_sample_counts[cSample]
    norm_counts[c] = median_sample_count * read_counts[c] / norm_factor
    log_counts[c] = np.log2(norm_counts[c]+1)
#print(log_counts)

#now cluster the data on samples and Vs
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

HCann = linkage(log_counts)
#print(HCann)
HCsam = linkage(log_counts.T)
#print(HCsam)

clustann = fcluster(HCann,t=0)
#print(clustann)
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(color_codes=True,font_scale=0.7)
sns_plot = sns.clustermap(log_counts,method='complete',cmap='Spectral_r')
plt.setp(sns_plot.ax_heatmap.get_yticklabels(),rotation = 0)
plt.setp(sns_plot.ax_heatmap.get_xticklabels(),rotation = 45,ha='right')
sns_plot.savefig("heatmap.pdf")
