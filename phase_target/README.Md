# Phase_target

You'll have a fasta file with the target loci that you are interested in pulling down from your sample. This might be a file that doesn't have multiple taxa (i.e. each target locus is present just once), or it may have multiple taxa per locus. This first little script is going to split out the fasta file into locus specific files.

### 1. Getting your locus-specific fasta files together
First, define your file name
```
fastafilename=target_Zhang_et_al.fa
```
These steps are the same regardless of whether your target fasta has the same locus more than once or not
```
nolines=`wc -l $allelefilename | awk '{print $1}'`
mkdir fasta_files
```
For files where the target locus is present just once
```
headercontents=`tail -n+1 $allelefilename | head -n1`
locus_name=`echo $headercontents | sed 's/>//g'`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g'`
echo $linecontents >> fasta_files/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta
```
For files where the target locus is present more than once you'll need to define how the sequence header combines the 
the sample name and the gene/locus name e.g. what delimits them, and what is the first field
```
delimiter="-"
first_field="sample"
# or alternatively first_field="locus"

if [ $first_field = "sample" ]; then

headercontents=`tail -n+1 $allelefilename | head -n1`
locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta

else

headercontents=`tail -n+1 $allelefilename | head -n1`
locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

fi
```

### 2. Using your locus specific fasta files as references to pull down the specific sequence for each sample
Stuff that needs to be installed: Java, GATK v4.0, GATK v3.8, picard (paths to these will be defined below). Also on your path you'll need bwa, samtools, and R. You will then need to set up a 'phasing_settings' file. On each separate line, in this order, you will need:

Line 1: path to your up-to-data jdk/bin/java. If your default 'java' command is pointing to the right version, you can just put in java on this line (I have had trouble with $JAVA_HOME on systems where the native shell is not bash)

Line 2: path to folder containing the gatk executable

Line 3: path to picard

Line 4: paired or single, depending on your sequencing method

Line 5: the pathway to your 'cleaned' F reads (or just your cleaned reads if single end). "${name}" should be used as a placeholder for the actual sample name wherever this occurs in the pathway to the reads (as the actual ${name} will be inserted by the script based on samples.txt. This program expects the path to the reads to be standard across your samples, so you will need to rename them if this is not the case

Line 6: the same thing for your 'cleaned' R reads if you have paired sequence data (blank if you don't)

Line 7: the path to GATK3.8/GenomeAnalysisTK.jar (unfortunately not all tools needed were ported to GATK v4)

Line 8: the number of threads you'd like bwa to use

e.g.
```
/nfs1/FW_HMSC/Baker_Lab/bin/jdk1.8.0_72/bin/java
/nfs1/FW_HMSC/Baker_Lab/bin/gatk_folder
/nfs1/FW_HMSC/Baker_Lab/bin/picard/dist/picard.jar
paired
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.1.fq.gz
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.2.fq.gz
/nfs1/FW_HMSC/Baker_Lab/bin/GenomeAnalysisTK.jar
4
```
Another example of the pathway to the reads where these are in a folder with the name of the sample:
```
/home/a499a400/Kaloula/cleaned-reads/${name}/split-adapter-quality-trimmed/${name}-READ1.fastq.gz
/home/a499a400/Kaloula/cleaned-reads/${name}/split-adapter-quality-trimmed/${name}-READ2.fastq.gz
```
An example with single end sequencing:
```
/usr/bin/java
/usr/local/gatk-4.0.2.1
/usr/local/picard-tools-2.18.0/picard.jar
single
/mnt/hcs-gemmell/beetles/reads/${name}-pe100-reads.fq.gz

/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
8
```

As well as your phasing_settings file, you'll need a file ("samples.txt") which lists the samples you are interested in assembling, one on each line. The names given here should match the ${name} given in phasing_settings e.g.
```
amphizoa
bemHap1
chlSer1
lioTuu1
omoHam1
pterMel1
traGib1
```
If you have paired end data, before you run the script for this step, make sure that your F and R files have the same read names in both (i.e. you cannot have a suffix of _1 in the read1 file, and _2 in the read2 file). If they have different names in each file, use sed to rename them e.g.
```
for i in `ls *.2.fq.gz`;
do gunzip $i;
newname=`echo $i | sed 's/.gz//g'`;
sed -i 's/_2$/_1/g' $newname;
gzip $newname;
done
```

-- allow for multiple mapping and then use the consensus to pull down the final locus




