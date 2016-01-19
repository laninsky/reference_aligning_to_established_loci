#Reference_aligning_to_established_loci

Starting with the *.alleles file produced by pyRAD, (1) produces one file per locus containing all the alleles present for each sample at that locus; (2) pulls one allele from each locus file for the species/samples you are interested in to act as a reference and stores these in a single fasta file; (3) carries out a reference-based alignment using bwa, gatk, samtools and picard using these reference loci; (4) combines these alleles for your new samples with the previous samples from step #2; (5) generates a summary of the data and SNP files for downstream use using code similar to that found in: https://github.com/laninsky/pyRAD_alleles_into_structure

#Step 1
Make sure to modify your allelefilename to what your allele file is actually called.
```
bash
allelefilename=c88d6m4p3.alleles
nolines=`wc -l $allelefilename | awk '{print $1}'`
i=1

for j in `seq 1 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $linecontents | awk '{print $1}' >> $i.fasta;
echo $linecontents | awk '{print $2}' >> $i.fasta;
else
i=$(expr $i + 1);
fi
done;
```

For step 2, you need a text file of the samples/species you want to use as references called "ref_samples.txt". For each locus, the "top-ranked" sample (the sample in line 1) will have one of its first allele pulled out and placed into a reference fasta file ("reference.fa"). If the top-ranked sample is not present, the script will check for the next sample and so on.

An example ref_samples.txt file (these samples had numerical codes)
```
28311
13193
28427
```

#Step 2
The script expects the first allele for each sample to be named samplename.assembled_0. If you have a different coding system for your alleles, make sure to tweak "samplenamesuffix" in the code below.
```
samplenamesuffix=.assembled_0

norefs=`wc -l ref_samples.txt | awk '{print $1}'`

for i in `ls *.fasta`;
do for j in `seq 1 $norefs`;
do refname=`tail -n+$j ref_samples.txt | head -n1`;
if grep --quiet $refname $i;
then echo ">"$i >> reference.fa; 
temp=`grep -A1 $refname$samplenamesuffix $i`;
echo $temp | awk '{print $2}' >> reference.fa;
break;
fi
done
done
```

#Step 3
This step is going to carry out the reference-guided assembly for your samples. You will need to have bwa, samtools, gatk and picard installed, with bwa and samtools in your path. To install gatk and picard, I did the following:

-- installed up to date java (jdk), declared JAVA_HOME in .bashrc, .cshrc etc. Added full path to java to $path/$PATH.

-- installed up to date apache ant, declared ANT_HOME. Added full path to ant to $path/$PATH

-- installed htsjdk using apache ant, and then copied directory to picard's folder after installation

-- installed picard using apache ant

-- unzipped the gatk binary

You will also need to set up a 'phasing_settings' file as per https://github.com/laninsky/phase_everyone#phasing_settings-file-example. Briefly, on each separate line, in this order, you will need:

Line 1: path to your up-to-data jdk/bin/java. If your default 'java' command is pointing to the right version, you can just put in java on this line (I have had trouble with $JAVA_HOME on systems where the native shell is not bash)

Line 2: path to gatk 

Line 3: path to picard

Line 4: paired or single, depending on your sequencing method

Line 5: the pathway to your 'cleaned' F reads (or just your cleaned reads if single end). Substitute "${name}" for the actual sample name wherever this occurs in the pathway to the reads. This program expects the to the reads to be standard across your samples, so you will need to rename them if this is not the case

Line 6: the same thing for your 'cleaned' R reads if you have paired sequence data
```
/nfs1/FW_HMSC/Baker_Lab/bin/jdk1.8.0_72/bin/java
/nfs1/FW_HMSC/Baker_Lab/bin/GenomeAnalysisTK.jar
/nfs1/FW_HMSC/Baker_Lab/bin/picard/dist/picard.jar
paired
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.1.gq.gz
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.2.gq.gz
```
Another example of the pathway to the reads where these are in a folder with the name of the sample:
```
/home/a499a400/Kaloula/cleaned-reads/${name}/split-adapter-quality-trimmed/${name}-READ1.fastq.gz
/home/a499a400/Kaloula/cleaned-reads/${name}/split-adapter-quality-trimmed/${name}-READ2.fastq.gz
```
