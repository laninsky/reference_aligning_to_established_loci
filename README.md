#Reference_aligning_to_established_loci

Starting with the *.alleles file produced by pyRAD, (1) produces one file per locus containing all the alleles present for each sample at that locus; (2) pulls one allele from each locus file for the species/samples you are interested in to act as a reference and stores these in a single fasta file; (4) carries out a reference-based alignment using bwa, gatk, samtools and picard using these reference loci; (5) combines these alleles for your new samples with the previous samples from step #2; (6) generates a summary of the data and SNP files for downstream use using code similar to that found in: https://github.com/laninsky/pyRAD_alleles_into_structure

#Step 1
Make sure to modify your allelefilename to what your allele file is actually called.
```
bash
allelefilename=c88d6m4p3.alleles
nolines=`wc -l $allelefilename | awk '{print $1}'`
i=1

for j in `seq 1 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep ">";
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
