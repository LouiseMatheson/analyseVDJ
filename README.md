# analyseVDJ
Tool for initial QC and analysis of VDJ-seq repertoire data

This tool takes in a set of bam files and an annotation file, and quantifies reads over the annotations. It then normalises, optionally log transforms, and clusters those datasets based on both samples (eg for checking whether replicates cluster) and annotations, and generates a heatmap to give an overview of the read counts for each sample. 

Whilst it was designed for immunoglobulin repertoire data, it could be used for other types of (single end) data, and options are given for same (default), opposing or non-strand-specific data.

If the data is VkJk-seq data, and has been analysed with the pipeline that we have developed (https://github.com/FelixKrueger/BabrahamLinkON/blob/master/run_VkSk-Seq_pipeline.md), then several reports will be generated that can additionally be analysed here to provide QC plots. Summary files produced by running the data through the IMGT High/V-QUEST tool  will also be identified, and plots generated for the CDR3 length distribution and proportion of productive vs non-productive rearrangements.

Since VkJk-seq bam files are typically split based on the Jk gene, sets of BAM files will also be identified that are predicted to be part of the same sample: the total read counts across all BAM files will then be used for normalisation. This identification is based on the string '\_L00': since by default the Jk gene is given after this string, files for which the name preceding it are identical are considered part of the same sample. If not included in the file name, sample groups cannot be identified.

The tool can be run from the command line with the command python analyseVDJ.py

By default, all .bam files in the current directory will be identified and used as input. The current directory will also be searched for .txt files that include the string 'annotation' or 'Annotation' to use as the annotation file. If more than one file is found, the user will be prompted to select the correct one. If none are found an error will be reported. Report/IMGT files will be identified from unique strings present if generated using our VkJk-seq analysis pipeline, or for IMGT Summary files \(\*1_Summary.txt) which should have an identifier appended to the start of the name to distinguish samples.

Additional directories to search can be listed with the -d flag. All files fulfilling the criteria above that are in the current directory or a directory listed here will be used as input. To specify a particular annotation file (need not contain 'A/annotation') or specific BAM files, please use the -a and -b flags respectively. These must either be in the current directory, or the complete path must be specified for each: only files that are listed will then be used as input and the -d flag will be ignored for that file type. 

Annotation files should contain columns labelled 'Chromosome', 'Start', 'End', 'Name' and (if applicable) 'Strand'. Some common variants will be corrected, but if the columns cannot be identified, except in the case of 'Name', an error will be generated. 

If the 'Name' column is not specified: if no other columns (in addition to those expected) are present, the row index will be used as the name. If only one additional column is present then this will be used. If more than one additional column is available, the user will be prompted to select one of these. 

Total read counts for each sample will be calculated, and samples normalised the the median total count. Log transformed read counts will also be generated as log2(read counts + 1).

The tool will provide a text file as output, listing eg the samples identified, total read counts. The raw read counts, normalised read counts and log transformed read counts will also be included. A prefix for this file, and the PDF files below, can be specified with the -p flag; the text output will be named prefix_out.txt and saved to the current directory.

In addition, several PDF files will be generated (depending on whether the sample is a VkJk-seq sample, and whether specific reports can be found).

For all:
prefix_heatmap.pdf - a clustered heatmap of read counts for each sample and annotation. These can be log2
transformed or untransformed (specify -l F)

For VkJk-seq:

prefix_pipeline.pdf - read counts over each stage of the VkJk-seq analysis pipeline (depending on reports found may not include all possible stages)

prefix_Jk_before_pipeline.pdf - proportion of reads associated with each Jk 'bait' sequence before running through pipeline (if sorting report provided)

prefix_Jk_after_pipeline.pdf - proportion of unique Vk-recombined reads associated with each Jk gene

prefix_IMGT_functionality_and_CDR3.pdf - proportion of productive vs non-productive rearrangements and distribution of CDR3 lengths (if IMGT summary files provided)

Four test bam files are provided, with corresponding annotation and report files, and the output files generated if run with the command: python analyseVDJ.py -r -p test

All possible input options are provided here:

optional arguments:

  -h, --help            show this help message and exit
  
  -s STRAND, --strand STRAND
                        specify whether library is strand specific, one of: 'same'
                        (default), 'opposing' or 'unstranded'
  
  -a ANNOTATION, --annotation ANNOTATION
                        annotation file to use; should be a .txt file with
                        columns for Chromosome, Start, End, Name and Strand
                        (if applicable). Please provide complete path unless
                        in current directory. By default will search for a
                        .txt file containing the word 'annotation' or
                        'Annotation'. Chromosome names in annotation file must
                        match those in BAM file, and if format is not 1, 2,
                        3...X, Y or chr1, chr2 etc please provide a list using -c flag
  
  -r--reduceNames       reduce length of sample names: will remove lane
                        number, extension, barcode, L00n.R, unique_V from
                        sample names before plotting
  
  -k KAPPA, --kappa KAPPA
                        specify F if this is not a VkJk-seq library, or if
                        reports are not present
  
  -b [BAM [BAM ...]], --bam [BAM [BAM ...]]
                        optionally list all bam files that should be used (if
                        this option is specified, please list every file
                        required and provide complete paths unless in current
                        directory). If none are listed, all .bam files in the
                        current directory and those listed with -d will be
                        used
  
  -d [DIRECTORIES [DIRECTORIES ...]], --directories [DIRECTORIES [DIRECTORIES ...]]
                        Provide paths to any additional directories to search
                        for BAM/annotation/report files. Does not apply for
                        files listed with -a or -b: these must be in the
                        current directory or the complete path must be
                        specified. The current directory will always be
                        searched, in addition to any specified here  

  -p PREFIX, --prefix PREFIX
                        Please specify a prefix to append to files generated,
                        if required
  
  -l LOG, -log LOG      Specify F if read counts used in heatmap should not be
                        log transformed
                        
   -c [CHROMOSOMENAMES [CHROMOSOMENAMES ...]], --chromosomeNames [CHROMOSOMENAMES [CHROMOSOMENAMES ...]]
                        If format of chromosome names is not 1,2,3...X,Y,MT or
                        chr1, chr2 etc (max number 22), please provide the
                        names of each chromosome name present in the
                        annotation or BAM file, eg seq1 seq2 seq3 etc
