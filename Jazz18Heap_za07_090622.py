#!/usr/bin/env python -i

##      Jazz18Heap.py version za07 090622
## Jazz18Heap a fast Python Script for finding matches to a probe sequence in large numbers of
##    High throughput sequencing output files
## This program implrements a simple index search that will look for instances of sequence reads matching
##    a reference in at least one k-mer.
## Jazz18Heap is intended to look for evidence of matches to a relatively short reference sequence
##    e.g., less than 100KB, but probably workable up to several MB in a large number of high throughput sequencing experiments
## Jazz18Heap is intended for finding relatively rare sequence matches (not common matches)
## Jazz18Heap doesn't substitute for the many tools that align and track coverage.  
## Jazz18Heap main goal is rapid identification of homologous sequences (homology defined as a perfect match to a long sequence-e.g. 32b for klen=32)
##
## Inputs are as follows (command line, Key=Value syntax)
##     ReferenceFile = <FastA file with list of sequences to match k-mers from>
##     ExcludeFile = <FastA file with sequences that will be excluded from matches>
##     DataFiles = <List of .fastq files, .fasta files, or NCBA-SRA accessions [e.g, SRR####]
##        Lists are comma delimited with no spaces in list
##           .fasta or .fastq files can be gzip compressed, although this somewhat slows the program down
##           .fasta and .fasta.gz files must have exactly one sequence per line (no multiline sequences)
##           For .fasta files to be downloaded from NCBI-SRA archive this entails the command line parameter --fasta 0
##        * Wildcards are allowed here as well, or list of files in a file with the extension .files
##        Providing a directory here will search all files in this directory or subdirectory for fasta and fastq data files
##
## Optional Parameters (will default to reasonable values if not set)
##     OutFileBase = <Character String to Label Output Files with>
##     ReportGranularity = <How many reads to process before reporting hit numbers (default is 1Million)>
##     SearchGranularity = <How much distance between k-mers to be examined in each data read
##         setting SearchGranularity=1 makes Jass18Heap look at every k-mer in every read
##         setting SearchGranularity=8 makes Jass18Heap look at every 8th k-mer in every read
##         higher numbers may miss a few hits but can greatly improve speed.
##         setting to a large number (999999) ensures only one k-mer will be looked up per read
##     SearchOffset = <Where to start jumping through each read for potential k-mers (zero means start at first base)
##     fastqdump = <Where to look for fastq-dump -- this is only important if you are downloading files from NCBI SRA during execution of Jazz18Heap>
##     klen = <How long are the k-mers used>.  Default klen = 32 
##     snpAllow = <Set to True to allow a single mismatch in each k-mer (default is False)
##     SeparateOutput = <Set to True to output CaughtRead fastA separately for each input file
##     Circular = <Set to True to force every Reference sequence to be treated as a circle>
##        Default : Uses the FastA name line-- if this line contains "Circular", the sequence is treated as a circle
##     FileByFileCounts = <Set True [default] to have a number of matches per reference file recorded for each sequence data file>
##     RefByRefCounts = <Set True to obtain a 'Ref-By-Ref' output that lists the numbers of matches for each reference sequence in each scanned file>
##        There are various setting for RefByRefCounts and FileByFileCounts
##          RefReads=True counts the number of matching reads for each reference sequence
##          RefKmers=True counts the number of matching kmers for each reference sequence
##          RefKSpecies=True counts the number of matching kmer species for each reference sequence
##           This is generally useful if the reference list is a set of fixed length k-mers and the intended output is for a count of each of these k-mers in each of the datafiles
##           Generally this will also involve setting SearchGranularity = 1 (so searching for every occurence of every k-mer)
##           This is generally useful if the reference list is a set of fixed length k-mers and the intended output is for a count of each of these k-mers in each of the datafiles
##     For RefByRef Output you can set two modes
##        RefColumns=True yields an output where each reference sequence is represented by one column (or two if both reads and kmers are counted)
##        RefRows=True yields an output where each reference sequence is represented by one row
##     RunsToSamples = Allows multiple runs to be associated with sample names.
##        Input is a file for which each line consists of a sample name (e.g. N2) followed by a tab and a comma-delimited list of run names (e.g. SRR####)
##     CaughtReads = <Set True [default] to obtain a list of reads matching a k-mer during the scan as a fasta file>
##     FileByFileMatches = <Set True [default] to have a list of individual matches in the FileByFile output (useful but somewhat typographically cumbersome>
##     MaxFileReads = <Default is no maximum>  This sets the maximum number of reads to process for each data file (e.g. MaxFileReads=100 looks at 100 reads)
##     SubKLen = <Default is (3,6,9)> Specifies the scales on which the individual k-mers are filtered to remove noisy degenerate sequences.  
##         SubKLen=3 instructs the program to check all 3-mers in each searched word to see how many match each other.  If they all match (or a high fraction do
##         Then the k-mer is not particularly complex (ie the same 3-mers over and over again) and will likely get some background.  Thresholds by default are 
##         No more than 18% of possible coincidences between 3-mers are indeed the same, no more than 5% of 6-mers are coincidental and no more than 2% of 9-mers.  
##         These numbers were chosen so random sequences with reduced base composition (e.g. 100%AT) are still likely (>95%) to be marked as non-degenerate, while
##         Many different kinds of degenerate sequence (sequences with long homopolymer runs or repeats of a simple sequence) are likely to be identified as degenerate
##         and filetered out of the relevant dataset.  The Tolerances (default (0.18, 0.05, 0.02) for 3-mer, 6-mer, and 9-mer coincidence frequency thresholds can be adjusted 
##         with the command line SimpleFilter=(new values -- same number as the number of different SubKLengths)
##
## Multithreading: Jazz18Heap has the very primative multitasking ability to spawn a number
##     of derivative processes for a large number of data files to be scanned.  To use 16 Threads
##     Set Threads=16 in the command line.
##
## Output:
##     Output consists of a file-by-file list of hits, a log file with information on the run, and a FastA file with "Caught" Reeads
##     File-By-File Output Format
##       Column 1: Dataset_ReadNumber (e.g. SRRXXX_2 is the read 2 file from SRXXX)
##       Column 2: Number of Reads extracted and analyzed
##       Column 3: Number of bases analzyed
##       Column 4: Number of hits from the viral sampler query
##       Column 5: Total number of hit k-mers
##       Column 6: Summary of hits.  Each element has 
##            -> a query name (e.g. X51522.1_phageP4)
##            -> a position (e.g. p5012 starts at position 5012)
##            -> a strand ('s' indicates a sense match starting at base , 'a' indicates an antisense match, and 
##            -> a k-mer match count  (e.g. m4 indicates that four different k-mers were matched every nth k-mer being checked [n is the search granularity])
##
##     A second optional tabular output can be obtained either with each query sequence getting one line or with each biological sample getting one line
##       User request of this optional second tabular output is through the command line setting RefByRef=True
##       This output also allows the user to input a table which maps read file names to samples which may have been looked at in several different read files
##           Such an output file is specified by the option RunsToSamples=<FileName>.
##           The file format for RunsToSamples is to have a series of lineswhich have samplenames at a start followed by tab and a comma delimted series of read file names
##                example N2<tab>SRR133242,SRR42567,SRR51242
##                Spaces should be avoided in this list and spaces and commas avoided in the file names
##                Read File names are automatically grouped accoring to the "basename" of the read file (everything before a period or underscore)
##                Thus files SRR133242_1.fasta.gz and SRR133242_2.fasta.gz will be grouped and should be referred to as simply
##                SRR133242 for the purposes of the RunsToSamples input file
##       RefByRef output format lists can be in a format where each reference name is a column (RefColumns=True)
##           or in a format where each reference name is a row (RefRows=True)
##
##     FastA "Caught Read" files are a second output, with information about each read in the ID line
##         ID Structure:
##              DataFileName;LineNumberInDataFile;ReferenceFileName;ReferenceSequenceName;NumberOfMatchedKmers
##    
## LineEndings: User can set line endings in the command line for Unix or Windows
##     For Unix/Linux line endings set LineEndings=Unix ('\n')
##     For PC/Windows line endings set LineEndings=Windows ('\r\n')
##     For '\r' line endings set LineEndings=r
##
## Running the program:
##   Jazz18Heap only runs at full speed with a variant of Python (pypy) that included a just-in-time compiler
##   Syntax Jazz18Heap RefFile=<MyRefFile> DataFiles=MyFastA1.fasta,MyFastA2.fasta.gz,MyFastQ*.fastq <Other_Parameters> 
##
##  Copyright 2020-2022, Andrew Fire and Stanford University
##  With thanks to Dae-Eun Jeong, Loren Hansen, Matt McCoy, Nimit Jain, Massa Shoura, Karen Artiles, Lamia Wahba, Nelson Hall
##  Revision List (partial)
##
##  9/6/20 Version YA has some major shifts-- see below
##    Added count of K-mer species that were hit to delineate "one hit wonder" cases where a single query kmer (or a few) were hit many times 
##    Added a simple sequence filter using an algorithm inspired by "Dust" (David Lipman and Roman Tatusov)
##      Each k-mer is broken into consecutive sub k-mer of fixed length SubKLen1.  Default is SubKLen1 is 3
##      Then calculate for each k-mer S[k]=
##        Probability that two (different) and randomly chosen SubKLen substrings of K will be identical
##           (Number of coincidences between difference SubKLen substrings of K) / (Total Possible Coincidences)
##           If this is above a Threshhold (SimpleFilter1), then reject K-mer K as being simple-sequence
##    Changed output for multiple sequence queries to credity k-mer matches to each query instead of just first match per read
##  7/11/21 Fixed a bug in reporting column by column k-mers (no issues if the code ran previously... this prevented from running with certain settings)
##  7/12/21 Fixed a major bug in dust filtering.  The resultant is a tiny bit slower by default but much better at filtering out spurious matches
##  7/12/21 Reset the reporting of positions to report the first position in the reference sequence where the indicated read matches (previous versions
##            extrapolated back from the first match by the number of unmatched bases at the beginning of the relevant read.
##            Currently a value of 12245s means that the first perfectly matching k-mer starts at position 12245 on the sense strand
##            A value of 12245a means that the first perfectly matching k-mer is the reverse complement of the k-mer starting at position 12245
##            And reported values are 1-based (not zero based), so matching at the first base will have a value at 1a.
##  7/12/21 The program will now do a k-mer by k-mer analysis of any read that has a positive (so the granularity is relevant to identify reads that might be positive
##            Then the program goes back through on a k-mer by k-mer basis to know exactally how many k-mers match and precisely where they start

import gzip, os
from sys import argv, version, executable
from time import time, strftime, localtime, asctime
from itertools import chain, product
from glob import glob
import subprocess
from random import choice


t0 = time()
now1 = strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())
Fn0 = []  ## List of reference files to search for sequences from
Gn0 = []  ## List of files containing sequences to exclude from k-mer match 
DataFile0 = []
OutFileBase0 = ''
klen1 = 0
ReportGranularity1 = 1000000  ## How often to report granularity, zero to provide less info in reports
SearchGranularity1 = 8 ## Granularity of search-- setting this to 1 searches ever k-mer in data files, setting to 10 will search every 10th k-mer
SearchOffset1 = 0 ## Granularity of search-- setting this to 1 searches ever k-mer in data files, setting to 10 will search every 10th k-mer
OneSnp1 = False ## Allow a single snp in each k-mer?
Circular1 = False ## Assume every sequence is circular?

MyCode1 = argv[0]
ai1 = 1

FastQDumpProgram1 = 'fastq-dump' ##'fasterq-dump' ## Location of a program that will download files from SRA if needed, can be reset by setting FastQDump=<path to fast(er)q-dump executable>.

MultiTask1 = False
BlockMultiCall1 = False
MultiCall1 = False
Multi1 = 0
Multi2 = 1
RefByRef1 = False
RefRows1 = False
RefColumns1 = True
RefReads1 = True
RefKmers1 = False
RefKSpecies1 = True
RunToSample1 = False
RunToSampleNumberD1 = {}
MaxFileReads1 = 0
SampleNumberToColumnNumberD1 = {}
AllSampleNameArray1 = []
SimpleFilter1 = (0.18,0.05,0.02) ## A simple filter that removes k-mers with limited complexity
SubKLen1 = (3,6,9)
## The simplicity of any given k-mer is set to the probablility of any two arbitrary subk-mers from the k-mer being equal
## This is 1.0 for a homopolymer and 0 for a k-mer where no sub-k-mer is repeated
SeparateOutput1 = False ## Set to true to provide separate "CaughtReads" output for each Read file
CaughtReads1 = True
FileByFileMatches1 = True
FileByFileCounts1 = False

LineEnding1 = '\n'

while ai1<len(argv):
    a1=argv[ai1]
    ai1+=1
    a1=a1.replace('"','').replace("'","")
    if a1[0]=='-':
        a11=a1.strip()[1:].lower()
        a22=argv[ai1].strip().replace('"','').replace("'","")
        ai1+=1
    else:
        a11=a1.split('=')[0].strip().lower()
        a22=a1.split('=')[-1].strip()
    if (a1.startswith('help') or ('-help' in a1)) and Multi1==0:
        mycode1 = open(argv[0]).read().splitlines()
        for h1 in mycode1:
            if h1.startswith('##'):
                print(h1[2:])
            elif h1.startswith('import'):
                break
        exit()
    if a11.startswith('reference') or a11.startswith('include'):
        Fn0.extend(a22.split(','))
    if a11.startswith('exclude'):
        Gn0.extend(a22.split(','))
    elif a11.startswith('data'):
        for a222 in a22.split(','):
            if os.path.isdir(a222):
                for D1,S1,FList1 in os.walk(a222):
                    for fc1 in FList1:
                        if fc1.endswith('.fastq') or fc1.endswith('.fastq.gz') or fc1.endswith('.fastq') or fc1.endswith('.fasta.gz'):
                            DataFile0.append(os.path.join(D1,fc1))
            elif '*' in a222:
                DataFile0.extend(glob(a222))
            elif a222.endswith('.files') or a222.endswith('.files.txt'):
                a2222 = open(a222,mode='r').read().split()
                for fc1 in a2222:
                    fc1 = fc1.strip()
                    if not(fc1):
                        continue
                    elif '*' in fc1:
                        DataFile0.extend(glob(fc1))
                    elif os.path.isdir(a222):
                        for D1,S1,FList1 in os.walk(a222):
                            for fc11 in FList1:
                                if fc11.endswith('.fastq') or fc11.endswith('.fastq.gz') or fc11.endswith('.fastq') or fc11.endswith('.fasta.gz'):
                                    DataFile0.append(os.path.join(D1,fc11))                    
                    else:
                        DataFile0.append(fc1)
            else:
                DataFile0.append(a222)
    elif a11.startswith('out'):
        OutFileBase0 = a22
    elif a11.startswith('lineend') or a11.startswith('linend'):
        if a22:
            c = a22[0].lower()
            if c in 'lumnafo': ## \n, linux, unix, mac os x, apple, android, freebsd
                LineEnding1 = '\n'
            elif c in 'pwdsb': ## \both, pc, windows, symbian, palm, dos
                LineEnding1 = '\r\n'
            elif c in 'zcrt92': ## \r, Commodore, Clasic Mac (OS9), Apple 2, zx-spectrum, TRS-80
                LineEnding1 = '\r'
    elif a11.startswith('k'):
        klen1 = int(a22)
    elif a11.startswith('reportgran'):
        ReportGranularity1 = int(a22)
    elif a11.startswith('maxfile') or a11.startswith('maxread'):
        MaxFileReads1 = int(a22)
    elif a11.startswith('searchgran'):
        SearchGranularity1 = int(a22)
    elif a11.startswith('searchoff'):
        SearchOffset1 = int(a22)
    elif a11.startswith('cir') :
        if a22.lower().startswith('f'):
            Circular1 = False
        else:
            Circular1 = True
    elif a11.startswith('caught') :
        if a22.lower().startswith('f'):
            CaughtReads1 = False
        else:
            CaughtReads1 = True
    elif a11.startswith('filebyfilematches') :
        if a22.lower().startswith('f'):
            FileByFileMatches1 = False
        else:
            FileByFileMatches1 = True
    elif a11.startswith('filebyfilecounts') :
        if a22.lower().startswith('f'):
            FileByFileCounts1 = False
        else:
            FileByFileCounts1 = True
    elif a11.startswith('runtosample') or a11.startswith('runstosample') or a11.startswith('samplestorun') or a11.startswith('sampletorun'):
        for L0 in open(a22,mode='rt'):
            L1 = L0.strip().split('\t')
            SampleName1 = L1[0]
            AllSampleNameArray1.append(SampleName1)
            for L11 in L1[1:]:
                for run1 in L11.split(','):
                    if run1:
                        RunToSampleNumberD1[run1] = len(AllSampleNameArray1)-1
        RunToSample1 = True
    elif a11.startswith('snp') :
        if a22.lower().startswith('f'):
            OneSnp1 = False
        else:
            OneSnp1 = True
    elif a11.startswith('simple') :
        if a22.lower().startswith('f'):
            SimpleFilter1 = []
        elif not(a22.lower().startswith('t')):
            SimpleFilter1 = [float(a222) for a222 in a22.split(',')]
    elif a11.startswith('subk') :
        SubKLen1 = [int(a222) for a222 in a22.split(',')]
    elif a11.startswith('seperate') or a11.startswith('separate'):
        if a22.lower().startswith('f'):
            SeparateOutput1 = False
        else:
            SeparateOutput1 = True
    elif a11.startswith('matchtable') or a11.startswith('refbyref'):
        if a22.lower().startswith('f'):
            RefByRef1 = False
        else:
            RefByRef1 = True
    elif a11.startswith('refreads'):
        if a22.lower().startswith('f'):
            RefReads1 = False
        else:
            RefReads1 = True
    elif a11.startswith('refkmers'):
        if a22.lower().startswith('f'):
            RefKmers1 = False
        else:
            RefKmers1 = True
    elif a11.startswith('refkspecies'):
        if a22.lower().startswith('f'):
            RefKSpecies1 = False
        else:
            RefKSpecies1 = True
    elif a11.startswith('refcolumns'):
        if a22.lower().startswith('f'):
            RefColumns1 = False
            RefRows1 = True
        else:
            RefColumns1 = True
            RefRows1 = False
            RefByRef1 = True
    elif a11.startswith('refrows'):
        if a22.lower().startswith('f'):
            RefRows1 = False
            RefColumns1 = True
        else:
            RefRows1 = True
            RefColumns1 = False
            RefByRef1 = True
    elif a11.startswith('multi') or  a11.startswith('thread') :
        Multi1 = 0
        Multi2 = int(a22)
        MultiCall1 = True
    elif a11.startswith('instance'): 
        if '/' in a22:
            Multi1 = int(a22.split('/')[0])-1
            Multi2 = int(a22.split('/')[1])
            MultiTask1 = True
            MultiCall1 = False
            BlockMultiCall1 = True
    elif a11.startswith('fastqdump') or a11.startswith('fastq-dump'):
        if os.path.isfile(a22) and not(os.path.isdir(a22)):
            FastQDumpProgram1 = a22
        elif os.path.isdir(a22):
            FastQDumpProgram1 = os.path.join(a22,'fastq-dump')
            if not(os.path.isfile(FastQDumpProgram1)):
                FastQDumpProgram1 = os.path.join(a22,'bin','fastq-dump')
                if not(os.path.isfile(FastQDumpProgram1)):
                    os.environ["PATH"]=os.getenv("PATH")+':'+a22
    elif a11.startswith('fasterqdump') or a11.startswith('fasterq-dump'):
        if os.path.isfile(a22) and not(os.path.isdir(a22)):
            FastQDumpProgram1 = a22
        elif os.path.isdir(a22):
            FastQDumpProgram1 = os.path.join(a22,'fasterq-dump')
            if not(os.path.isfile(FastQDumpProgram1)):
                FastQDumpProgram1 = os.path.join(a22,'bin','fasterq-dump')
                if not(os.path.isfile(FastQDumpProgram1)):
                    os.environ["PATH"]=os.getenv("PATH")+':'+a22
MultiMnemonic1 = 'J18H '
if MultiTask1 or MultiCall1:
    MultiMnemonic1 = 'J18H '+str(Multi1+1)+'/'+str(Multi2)+' '
if MultiCall1 and not(BlockMultiCall1):
    for j in range(1,Multi2):
        subprocess.Popen([executable,]+argv+['instance='+str(j+1)+'/'+str(Multi2)])
        
                                 
##Default Values for Debugging
if not(OutFileBase0):
    OutFileBase0 = 'Jazz18Heap_I'

OutFileBase0 += '{:02d}'.format(Multi1)+'_'+now1

LogFile1="LogFile_"+OutFileBase0+'.tdv'
def LogNote1(note):
    LogFile=open(LogFile1,mode='a')
    LogFile.write(note+'\t'+'; Time='+"{0:.2f}".format(time()-t0)+' sec'+'\t'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())+' '+LineEnding1)
    LogFile.close()
    print(note.split('#')[0].replace('\t',' ').strip(',') + '; Time='+"{0:.2f}".format(time()-t0)+' sec')
##Default Values for Debugging
def FileInfo1(FileID):
    if type(FileID)==str:
        Fn1=FileID
    else:
        Fn1=FileID.name
    s11=os.stat(Fn1)
    return ','.join([Fn1,
                     '#',
                      'Path='+os.path.abspath(Fn1),
                        'Size='+str(s11[6]),
                        'Accessed='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[7])),
                        'Modified='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[8])),
                        'Created='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[9])),
                        'FullInfo='+str(s11)])
LogNote1("Running Jazz18Heap with parameters:"+' '.join(argv)+' #Python Version'+version)
LogNote1("Python Flavor/Version: "+version)

if not('pypy' in version.lower()):
    LogNote1("Note that you should run Jazz18Heap with PyPy ('www.pypy.org') for full speed (otherwise speeds may be 10x sloer)") 
def LogOpeningFile1(FileID):
    LogNote1("Opening "+FileInfo1(FileID))
def LogClosingFile1(FileID):
    LogNote1("Closed "+FileInfo1(FileID))
def LogRunningFile1(FileID):
    LogNote1("Running "+FileInfo1(FileID))
LogRunningFile1(MyCode1)
if not Fn0:
    Fn0 = ['nematode_polintons_reduced060921.fasta',]
    LogNote1('Reference File List set as debugging default to: '+str(Fn0)+ '. <RefFile= option on command line>')
if not(DataFile0):
    DataFile0 = ['/Users/firelab08/Desktop/SRR3083808.fasta.gz',]
    LogNote1('Data File List set as debugging default to: '+str(DataFile0)+ '. <DataFile= option on command line>')
if not(klen1):
    klen1 = 32
    LogNote1('k-mer length set as default to: '+str(klen1)+ '. <Klen= option on command line>')


D0 = {}  ## keys are k-mer sequences, values are tuples
         ## 0: Ordinal number of sequence that is hit
         ## 1: Ordinal position in sequence (one based +1 for sense, - for antisense)
         ## 2: Number of mismatches (presently 0 or 1)
NameList1 = []
FileList00 = []
SeqList1 = []
LenList1 = []
AntiList1 = []
snpD1 = {'G':'ATC','A':'TCG','T':'CGA','C':'GAT'}

## All possible SubKLen1 kmers are keys in the dictionary with the values being the ordinal
def versioner1(v88):
    s88 = ['0']
    for c1 in v88:
        if c1.isdigit():
            s88[-1]+=c1
        else:
            if s88[-1]!='0':
                s88.append('0')
    return(map(int,s88))
fqdVersion1 = ['0']

def SimpleDust1(k,SubKLen,SimpleFilter):
	## k is any k-mer (e.g., 32-mer), SubKLen is a list of lengths to do dust filtering at (e.g., 3,6,9), SimpleFilter is a set of thresholds (coincidence frequencies) above which a value of True  will be returned
	## The function calculates how many sub-k=mers of each length match other sub-k-mers (total coincidences as a function of the maximum possible number
	## a homopolymer k-mer k will have coincidences==possible-coincidences and c/mc will be 1.0.  A sequence with no sub-k-mer appearing more than once will have a c/mc of 0.0
	## sequences of intermediate complexity will have an intermediate c/mc
	## Heuristic exploration leads to a provisional recommendation of subkmer lengths of SubKLen=(3,6,9) for 32-mer words, with thresholds >0.18, >0.05, >0.02 for calling a given k-mer degenerate
	if not(SimpleFilter):
		return False
	for z in range(len(SubKLen)):
		d = {}
		n = len(k)-SubKLen[z]+1 ## number of SubKmers in each Kmer
		for j in range(n):
			cd1 = k[j:j+SubKLen[z]]
			if not cd1 in d:
				d[cd1] = 0
			d[cd1]+=1
		c = sum([x**2 for x in d.values()])-n ## coincidences
		mc = n**2-n ## possible coincidences
		if c>SimpleFilter[z]*mc:
			return True
	return False

def antisense1(seq):
    return seq.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()[::-1]
for nF0 in Fn0:
    Fx0 = open(nF0,mode='rU')
    LogOpeningFile1(Fx0)
    F0 = chain(Fx0,['>',])
    s1 = []
    n1 = os.path.basename(Fx0.name).split('.')[0]
    iN0 = 0
    toss1 = 0 ## k-mer sequences that flunked the degeneracy filter 
    for L0 in F0:
        if L0[0]=='>':           
            S1 = ''.join(s1).upper().replace('N','')  ##for now delete n's.  this willalso change some base numbers
            if S1:
                NameList1.append(n1)
                FileList00.append(nF0)
            if len(S1)>=klen1:
                if ("circular" in n1.lower()) or Circular1:
                    S1 = S1+S1[:klen1-1]  #assume everything is a circle
                lD1 = len(S1)
                A1 = antisense1(S1)
                N1 = len(NameList1)-1
                LenList1.append(lD1)
                SeqList1.append(S1)
                AntiList1.append(A1)
                for i in range(len(S1)-klen1+1):
                    K1 = S1[i:i+klen1]
                    if not(SimpleDust1(K1,SubKLen1,SimpleFilter1)) and (not(K1 in D0) or (D0[K1][2]!=0)):
                        D0[K1] = (N1,i+1,0)
                        D0[antisense1(K1)] = (N1,-i-1,0) ## antisense of a k-mer at position p will have the position -p.  Note all positions are one-based 
                    if OneSnp1:
                        for j in range(klen1):
                            for k in range(3):
                                mutK1 = K1[:j]+snpD1[K1[j]][k]+K1[j+1:]
                                if not(mutK1 in D0) and not(SimpleDust1(mutK1,SubKLen1,SimpleFilter1)):
                                    D0[mutK1] = (N1,i+1,1)
                                    D0[antisense1(mutK1)] = (N1,-i-1,1)
            n1a = L0[1:].strip().split()
            n2a = []
            for n0a in n1a:
                if not("=" in n0a) or ("range=" in n0a):
                    n2a.append(n0a.replace("range=",""))
            n1 = '_'.join(n2a)
            s1 = []
        else:
            s1.append(L0.strip())                         
    Fx0.close()
    LogNote1('Added k-mer dictionary from '+nF0)
LogNote1('Finished making initial index')

for nF0 in Gn0:  ## currently set up to only remove perfect matches
    Fx0 = open(nF0,mode='rU')
    LogOpeningFile1(Fx0)
    F0 = chain(Fx0,['>',])
    s1 = []
    n1 = 'Default_'+nF0
    iN0 = 0
    for L0 in F0:
        if L0[0]=='>':           
            S1 = ''.join(s1)
            if len(S1)>klen1:
                if ("circular" in n1.lower()) or Circular1:
                    S1 = S1+S1[:klen1-1]  #assume everything is a circle
                A1 = S1.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()[::-1]
                for i in range(len(S1)-klen1+1):
                    if S1[i:i+klen1] in D0:
                        del(D0[S1[i:i+klen1]])
                    if A1[i:i+klen1] in D0:
                        del(D0[A1[i:i+klen1]])
        else:
            s1.append(L0.strip())                         
    Fx0.close()
    LogNote1('Masked from k-mers from '+nF0)
LogNote1('Finished filtering index')

TempFileUID1 = ''
for i in range(12):
    TempFileUID1 += choice('QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm')
ScratchPath1 = os.path.join(os.getcwd(), 'Jazz18HeapTempFiles'+TempFileUID1)

DataFile1 = []
if RunToSample1:
    for df0 in DataFile0:
        Run1 = os.path.basename(df0).split('.')[0].split('_')[0]
        if Run1 in RunToSampleNumberD1:
            SampleNumber1 = RunToSampleNumberD1[Run1]
            if SampleNumber1 % Multi2==Multi1:
                if os.path.isfile(df0) and not(df0.endswith('.sra')):
                    DataFile1.append(df0)
                else:
                    LogNote1(df0+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
                    LogNote1("Preparing to download sequence read set "+df0+" from NCBI")
                    try:
                        TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-V'])
                    except:
                        LogNote1("Searching for a version of fast(er)q-dump that will run; if this fails, you may need to redownload the program and unzip the archive, also add FastQDumpProgram=<path to program> to command line")
                        os.environ["PATH"]=os.getenv("PATH")+':./:/opt/local/bin:/opt/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin:/Applications:~/Downloads'
                        for fqd1 in os.getenv("PATH").split(':'):
                            if os.path.isdir(fqd1):
                                for fqd2 in os.listdir(fqd1):
                                    if fqd2.startswith('sratoolkit'):
                                        fqd3 = os.path.join(fqd1,fqd2)
                                        if os.path.isdir(fqd3):
                                            fqd4 = os.path.join(fqd3,'bin','fastq-dump')
                                            if os.path.isfile(fqd4):
                                                if versioner1(fqd2) > versioner1(fqdVersion1):
                                                    fqdVersion1 = fqd2
                                                    FastQDumpProgram1 = fqd4
                                                    TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-V'])
                    LogNote1("Trying presumed fast(er)q-dump program file located at "+FastQDumpProgram1)
                    if not(os.path.isdir(ScratchPath1)):
                        os.mkdir(ScratchPath1)
                    TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,
                                                             '--fasta',
                                                             '0',
                                                             '--origfmt',
                                                             '--outdir',
                                                             ScratchPath1,
                                                             df0])
                    LogNote1("Result of "+df0+" NCBI Download " + TryFastQDump1)
                    PresumptiveFilePath1 = os.path.join(ScratchPath1,df0+'.fasta')
                    if os.path.isfile(PresumptiveFilePath1):
                        DataFile1.append(PresumptiveFilePath1)
                    else:
                        LogNote1('Looks Like fast(er)q-dump failed for '+df0)
else:                
    for df0 in DataFile0[Multi1::Multi2]:
        if os.path.isfile(df0) and not(df0.endswith('.sra')):
            DataFile1.append(df0)
        else:
            LogNote1(df0+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
            LogNote1("Preparing to download sequence read set "+df0+" from NCBI")
            try:
                TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-V'])
            except:
                LogNote1("Searching for a version of fast(er)q-dump that will run; if this fails, you may need to redownload the program and unzip the archive, also add FastQDumpProgram=<path to program> to command line")
                os.environ["PATH"]=os.getenv("PATH")+':./:/opt/local/bin:/opt/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin:/Applications:~/Downloads'
                for fqd1 in os.getenv("PATH").split(':'):
                    if os.path.isdir(fqd1):
                        for fqd2 in os.listdir(fqd1):
                            if fqd2.startswith('sratoolkit'):
                                fqd3 = os.path.join(fqd1,fqd2)
                                if os.path.isdir(fqd3):
                                    fqd4 = os.path.join(fqd3,'bin','fastq-dump')
                                    if os.path.isfile(fqd4):
                                        if versioner1(fqd2) > versioner1(fqdVersion1):
                                            fqdVersion1 = fqd2
                                            FastQDumpProgram1 = fqd4
                                            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-V'])
            LogNote1("Trying presumed fast(er)q-dump program file located at "+FastQDumpProgram1)
            if not(os.path.isdir(ScratchPath1)):
                os.mkdir(ScratchPath1)
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,
                                                     '--fasta',
                                                     '0',
                                                     '--origfmt',
                                                     '--outdir',
                                                     ScratchPath1,
                                                     df0])
            LogNote1("Result of "+df0+" NCBI Download " + TryFastQDump1)
            PresumptiveFilePath1 = os.path.join(ScratchPath1,df0+'.fasta')
            if os.path.isfile(PresumptiveFilePath1):
                DataFile1.append(PresumptiveFilePath1)
            else:
                LogNote1('Looks Like fast(er)q-dump failed for '+df0)

if not(SeparateOutput1) and CaughtReads1:
    OutFile1 = open('CaughtReads_'+OutFileBase0+'.fa',mode='w') ## FastA File
    LogOpeningFile1(OutFile1)
F2 = open('FileByFile_'+OutFileBase0+'.tdv',mode='w') ## FastA File
F2.write('File')
if RunToSample1:
    F2.write('\tSample')
F2.write('\tTotalReads\tTotalBases\tHitReads\tHitKmers')
if FileByFileCounts1:
    for n1 in NameList1:
        if RefReads1:
            F2.write('\t'+n1+'_reads')
        if RefKmers1:
            F2.write('\t'+n1+'_kmers')
        if RefKSpecies1:
            F2.write('\t'+n1+'_kmerSpecies')
if FileByFileMatches1:
    F2.write('\tMatchList')
F2.write(LineEnding1)
OrientA1 = ['b','s','a']
TotalReads1 = 0
TotalHitReads1 = 0
TotalBases1 = 0
TotalHitKmers1 = 0
MatchList0 = []
MatchList1 = []
MatchList2 = []
SampleNumberList1 = []
for iD0,nD0 in enumerate(DataFile1):
    Roton1 = 1
    HotLine1 = 2
    if nD0.lower().endswith('fastq') or nD0.lower().endswith('fastq.gz'):
        Roton1 = 4
        HotLine1 = 2
    if nD0.lower().endswith('fasta') or nD0.lower().endswith('fasta.gz'):
        Roton1 = 2
        HotLine1 = 2
    if nD0.endswith('gz'):
        if version.startswith('2.'):
            F1 = gzip.open(nD0,mode='r')            
        else:
            F1 = gzip.open(nD0,mode='rt')
    else:
        if version.startswith('2.'):
            F1 = open(nD0,mode='rU')
        else:
            F1 = open(nD0,mode='r')
    if ReportGranularity1:
        LogOpeningFile1(nD0)
    F1r = F1 ## may be faster with .read().splitlines()
    ReportInterval1 = ReportGranularity1 * Roton1
    FileReads1 = 0
    FileHitReads1 = 0
    FileBases1 = 0
    FileHitKmers1 = 0
    RotonCounter1 = 0
    LineCounter1 = 0
    HitD1 = {}
    ReadFileMnemonic1 = os.path.basename(F1.name).split('.')[0]
    Run1 = ReadFileMnemonic1.split('_')[0]
    if Run1 in RunToSampleNumberD1:
        SampleNumber1 = RunToSampleNumberD1[Run1]
    else:
        RunToSampleNumberD1[Run1]=len(RunToSampleNumberD1)
        SampleNumber1 = RunToSampleNumberD1[Run1]
        AllSampleNameArray1.append(Run1)
    if SampleNumber1 in SampleNumberToColumnNumberD1:
        Column1 = SampleNumberToColumnNumberD1[SampleNumber1]
    else:
        Column1 = len(SampleNumberToColumnNumberD1)
        SampleNumberList1.append(SampleNumber1)
        SampleNumberToColumnNumberD1[SampleNumber1] = Column1
    if RefReads1:
        MatchList0.append([0] * len(NameList1))
        FileList0 = [0] * len(NameList1)
    if RefKmers1:
        MatchList1.append([0] * len(NameList1))
        FileList1 = [0] * len(NameList1)
    if RefKSpecies1:
        MatchList2.append([0] * len(NameList1))
        FileList2 = [0] * len(NameList1)
    if SeparateOutput1 and CaughtReads1:
        OutFile1 = open("JazzCaughtReads_"+ReadFileMnemonic1+"_"+now1+".fa", mode='w')
    nL1=-1
    FileReads1 = 0
    KMatchDict1 ={}
    for L1 in F1r:
        nL1+=1
        RotonCounter1 += 1
        if RotonCounter1 == HotLine1:
            RotonCounter1 = 0
            R1 = L1.strip()
            len1 = len(R1)
            m1 = 0
            FileReads1+=1
            KMatchDict2 = {}
            if MaxFileReads1 and (FileReads1>MaxFileReads1):
                break
            Caught1 = False
            for i in range(SearchOffset1,len1-klen1+1,SearchGranularity1):
                q1 = R1[i:i+klen1]
                if q1 in D0:
                    Caught1 = True
                    break
            if Caught1:
                R11 = [False]*len(R1)
                for i in range(SearchOffset1,len1-klen1+1):
                    q1 = R1[i:i+klen1]
                    if q1 in D0:
                        h,p1,snp1 = D0[q1]
                        if not h in KMatchDict1:
                            KMatchDict1[h]={}
                        if not(p1 in KMatchDict1[h]):
                            KMatchDict1[h][p1] = 0
                        KMatchDict1[h][p1] += 1
                        if not h in KMatchDict2:
                            KMatchDict2[h]=0
                        KMatchDict2[h] += 1
                        if m1 == 0:
                            if p1>0:
                                pa1 = p1
                            for m in range(i,i+klen1):
                                R11[m] = True
                        else:
                            if R11[i+klen1-1]:
                                R11[i+klen1] = True
                            else:
                                for m in range(i,i+klen1):
                                    R11[m] = True                               
                        if p1<0:
                            pa1 = p1
                        m1 += 1
                if RefReads1:
                    for h in KMatchDict2:
                        MatchList0[Column1][h]+=1
                        FileList0[h]+=1
                if RefKmers1:
                    for h in KMatchDict2:
                        MatchList1[Column1][h] += KMatchDict2[h]
                        FileList1[h] += KMatchDict2[h]
                FileHitReads1 += 1
                FileHitKmers1 += m1
                if CaughtReads1:
                    R111 = ''
                    for i in range(len(R1)):
                        if R11[i]:
                            R111+=R1[i].upper()
                        else:
                            R111+=R1[i].lower()
                    if p1>0:
                        OutFile1.write('>'+
                                   ReadFileMnemonic1+
                                   '_Read='+str(nL1//Roton1)+
                                   '_Source='+FileList00[h]+
                                   '_Seq='+NameList1[h]+
                                   '_Pos='+str(pa1)+'s'+
                                   '_KMatches='+str(m1)+
                                   LineEnding1)
                    else:
                        OutFile1.write('>'+
                                   ReadFileMnemonic1+
                                   '_Read='+str(nL1//Roton1)+
                                   '_Source='+FileList00[h]+
                                   '_Seq='+NameList1[h]+
                                   '_Pos='+str(-pa1)+'a'+
                                   '_KMatches='+str(m1)+
                                   LineEnding1)
                    OutFile1.write(R111+LineEnding1)
                if not h in HitD1:
                    HitD1[h]={}
                if p1>0:
                    if not((pa1,m1,1) in HitD1[h]):
                        HitD1[h][(pa1,m1,1)] = 0
                    HitD1[h][(pa1,m1,1)] += 1
                else:
                    if not((-pa1,m1,-1) in HitD1[h]):
                        HitD1[h][(-pa1,m1,-1)] = 0
                    HitD1[h][(-pa1,m1,-1)] += 1
            FileReads1 += 1
            FileBases1 += len1
        LineCounter1 += 1
        if ReportGranularity1 and LineCounter1==ReportInterval1:
            LineCounter1 = 0
            LogNote1(MultiMnemonic1+'Completed Read '+
                     str(1+nL1//Roton1)+
                     ' from '+
                     ReadFileMnemonic1 +
                     '  Bases='+
                     str(FileBases1)+
                     '  HitReads='+
                     str(FileHitReads1)+
                     '  HitKmers='+
                     str(FileHitKmers1))                     
    if ReportGranularity1:
        LogClosingFile1(F1)
    F1.close()
    if RefKSpecies1:
        for h in KMatchDict1:
            MatchList2[Column1][h] += len(KMatchDict1[h])
            FileList2[h] += len(KMatchDict1[h])
    if SeparateOutput1 and CaughtReads1:
        OutFile1.close()
    rdT1 = ''
    if FileByFileMatches1:
        for i00 in sorted(list(HitD1.keys())):
            rdT1+=NameList1[i00]+':'
            for (p1,m1,o1) in sorted(list(HitD1[i00].keys())):
                rdT1 += 'p'+str(p1)+OrientA1[o1]+'_m'+str(m1)
                if HitD1[i00][(p1,m1,o1)]>1:
                    rdT1 += '('+str(HitD1[i00][(p1,m1,o1)])+')'
                rdT1 += ','
            if rdT1[-1]==',':
                rdT1 = rdT1[:-1]
            rdT1 += '/'
    rdT1 = rdT1.strip().strip('/')
    rdT2 =''
    if rdT1:
        rdT2 = ' Matches='+rdT1
    LogNote1(MultiMnemonic1+'Finishing File '+str(iD0)+': '
             +ReadFileMnemonic1+
             ' Reads='+
             str(FileReads1)+
             ' Bases='+
             str(FileBases1)+
             ' HitReads='+
             str(FileHitReads1)+
             ' HitKmers='+
             str(FileHitKmers1)+
             rdT2)
    if ReportGranularity1:
        LogNote1(MultiMnemonic1+'Project Progress: '
         +'AllFilesSoFar'+
         ' Reads='+
         str(TotalReads1)+
         ' Bases='+
         str(TotalBases1)+
         ' HitReads='+
         str(TotalHitReads1)+
         ' HitKmers='+
         str(TotalHitKmers1))
    F2.write(nD0)
    if RunToSample1:
        F2.write('\t'+AllSampleNameArray1[SampleNumber1])
    F2.write('\t'+
             str(FileReads1)+
             '\t'+
             str(FileBases1)+
             '\t'+
             str(FileHitReads1)+
             '\t'+
             str(FileHitKmers1))
    if FileByFileCounts1:
        for i in range(len(NameList1)):
            if RefReads1:
                F2.write('\t'+str(FileList0[i]))
            if RefKmers1:
                F2.write('\t'+str(FileList1[i]))
            if RefKSpecies1:
                F2.write('\t'+str(FileList2[i]))
    if FileByFileMatches1:
        F2.write('\t'+rdT1)
    F2.write(LineEnding1)
    TotalReads1 += FileReads1
    TotalBases1 += FileBases1
    TotalHitReads1 += FileHitReads1
    TotalHitKmers1 += FileHitKmers1
    if nD0.startswith(ScratchPath1) and nD0.endswith('.fasta'):
        os.remove(nD0)
    
LogNote1(MultiMnemonic1+'Finishing Run: '
         +'AllFiles'+
         ' Reads='+
         str(TotalReads1)+
         ' Bases='+
         str(TotalBases1)+
         ' HitReads='+
         str(TotalHitReads1)+
         ' HitKmers='+
         str(TotalHitKmers1))
F2.write('AllFiles'+
             '\t'+
             str(TotalReads1)+
             '\t'+
             str(TotalBases1)+
             '\t'+
             str(TotalHitReads1)+
             '\t'+
             str(TotalHitKmers1)+LineEnding1)
LogClosingFile1(F2)
F2.close()
if RefByRef1:
    F3 = open('RefByRef_'+OutFileBase0+'.tdv',mode='w') ## List of Hit reference counts
    if RefRows1:
        F3.write('SeqID\tSequence')
        for sam1 in SampleNumberList1:
            if RefReads1:
                F3.write('\t'+AllSampleNameArray1[sam1]+'_reads')
            if RefKmers1:
                F3.write('\t'+AllSampleNameArray1[sam1]+'_kmers')
            if RefKSpecies1:
                F3.write('\t'+AllSampleNameArray1[sam1]+'_kmers')
        F3.write(LineEnding1)
        for i,(n,s) in enumerate(zip(NameList1,SeqList1)):
            F3.write(n+'\t'+s)
            for j in range(len(SampleNumberList1)):
                if RefReads1:
                    F3.write('\t'+str(MatchList0[j][i]))
                if RefKmers1:
                    F3.write('\t'+str(MatchList1[j][i]))
                if RefKSpecies1:
                    F3.write('\t'+str(MatchList2[j][i]))
            F3.write(LineEnding1)       
    if RefColumns1:
        F3.write('Sample\tReads\t')
        for n1 in NameList1:
            if RefReads1:
                F3.write('\t'+n1+'_reads')
            if RefKmers1:
                F3.write('\t'+n1+'_kmers')
            if RefKSpecies1:
                F3.write('\t'+n1+'_kmerSpecies')
        F3.write(LineEnding1)
        for j,sam1 in enumerate(SampleNumberList1):
            F3.write(AllSampleNameArray1[sam1]+'\t')
            for i in range(len(NameList1)):
                if RefReads1:
                    F3.write('\t'+str(MatchList0[j][i]))
                if RefKmers1:
                    F3.write('\t'+str(MatchList1[j][i]))
                if RefKSpecies1:
                    F3.write('\t'+str(MatchList2[j][i]))
            F3.write(LineEnding1)       
    F3.close()
		
if not(SeparateOutput1) and CaughtReads1:
    LogClosingFile1(OutFile1)
    OutFile1.close()        
LogNote1('Finished running '+' '.join(argv))
try:
    LogNote1(open(MyCode1,mode='r').read())
except:
    pass

     
    
        


        
        
        
