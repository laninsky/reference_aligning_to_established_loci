# Reference_aligning_to_established_loci

Starting with the *.alleles file produced by pyRAD, (1) produces one fasta file per locus containing all the alleles present for each sample at that locus; (2) pulls one allele from each locus file for the species/samples you are interested in to act as a reference and stores these in a single fasta file; (3) carries out a reference-based alignment using bwa, gatk, samtools and picard using these reference loci; (4) generates a summary of the data and SNP files for downstream use using code similar to that found in: https://github.com/laninsky/pyRAD_alleles_into_structure; (5) combines these alleles for your new samples with the previous samples from step #2; (6) aligns the alleles for your new samples with your previous samples using MAFFT; (7) pulls out the SNPs from each locus for each of your fasta files; (8) filters these SNPs by lineages you require to be present/missing data to just one SNP per locus; and finally (9) (optional) further filters the SNPs for your reference-aligned samples to meet a minimum depth requirement.

The phase_everyone folder is a little different in that it assumes you are reference-aligning samples to their own reference (e.g. aligning reads from Sample_A to loci derived from Sample_A in order to phase them), and doing this for all samples in the dataset (rather than the instructions in this README, where you are aliging Sample_B to loci derived from Sample_A and then adding Sample_A's phased loci to your dataset).

# Step 1
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

# Step 2
For step 2, you need a text file of the samples/species you want to use as references called "ref_samples.txt". For each locus, the "top-ranked" sample (the sample in line 1) will have one of its first allele pulled out and placed into a reference fasta file ("reference.fa"). If the top-ranked sample is not present, the script will check for the next sample and so on (therefore you want to rank the samples from most preferred reference through to least that you put in the ref_samples.txt file). If the samples you want to reference align require different references (e.g. they are different species), it would be better to run steps 2-4 separately for each sample.

An example ref_samples.txt file (these samples had numerical codes)
```
28311
13193
28427
```

Before you start this script, backup your locus specific fasta files to another folder (call it backup/*.fasta because we'll reference this in a few steps time), just in case. The script expects the first allele for each sample to be named samplename.assembled_0. If you have a different coding system for your alleles, make sure to tweak "samplenamesuffix" in the code below (e.g. changing it to 'unassembled' etc!)
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

# Step 3
If everything went OK with the previous script, and you have a copy of your locus specific fasta files in another location (call it backup/*.fasta because we'll use this further down), you can delete the locus specific fasta files in this folder (just to make ls etc a little quicker). This step is going to carry out the reference-guided assembly for your samples. You will need to have bwa, samtools, gatk and picard installed, with bwa and samtools in your path. To install picard, I did the following:

-- installed up to date apache ant, declared ANT_HOME. Added full path to ant to $path/$PATH

-- installed htsjdk using apache ant, and then copied directory to picard's folder after installation

-- installed picard using apache ant

-- unzipped the gatk folder

You will also need to set up a 'phasing_settings' file. On each separate line, in this order, you will need:

Line 1: path to your up-to-data jdk/bin/java. If your default 'java' command is pointing to the right version, you can just put in java on this line (I have had trouble with $JAVA_HOME on systems where the native shell is not bash)

Line 2: path to folder containing the gatk executable

Line 3: path to picard

Line 4: paired or single, depending on your sequencing method

Line 5: the pathway to your 'cleaned' F reads (or just your cleaned reads if single end). "${name}" should be used as a placeholder for the actual sample name wherever this occurs in the pathway to the reads (as the actual ${name} will be inserted by the script based on samples.txt (phasing.sh) or samplenames.txt (phase_everyone/phase_everyone.sh). This program expects the path to the reads to be standard across your samples, so you will need to rename them if this is not the case

Line 6: the same thing for your 'cleaned' R reads if you have paired sequence data

e.g.
```
/nfs1/FW_HMSC/Baker_Lab/bin/jdk1.8.0_72/bin/java
/nfs1/FW_HMSC/Baker_Lab/bin/gatk_folder
/nfs1/FW_HMSC/Baker_Lab/bin/picard/dist/picard.jar
paired
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.1.fq.gz
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.2.fq.gz
```
Another example of the pathway to the reads where these are in a folder with the name of the sample:
```
/home/a499a400/Kaloula/cleaned-reads/${name}/split-adapter-quality-trimmed/${name}-READ1.fastq.gz
/home/a499a400/Kaloula/cleaned-reads/${name}/split-adapter-quality-trimmed/${name}-READ2.fastq.gz
```

Example if using pre-version 4 of gatk (make sure to use the phase_everyone_pre_v4_gatk.sh versions)
```
/nfs1/FW_HMSC/Baker_Lab/bin/jdk1.8.0_72/bin/java
/nfs1/FW_HMSC/Baker_Lab/bin/GenomeAnalysisTK.jar
/nfs1/FW_HMSC/Baker_Lab/bin/picard/dist/picard.jar
paired
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.1.fq.gz
/nfs1/FW_HMSC/Baker_Lab/emma_temp/ddRAD_for_Alana/QC_Phred20_discardN/${name}.2.fq.gz
```

As well as your phasing_settings file, you'll need a file ("samples.txt") which lists the samples you are interested in assembling, one on each line ("samplenames.txt if running the phase_everyone/phase_everyone.sh script). The names given here should match the ${name} given in phasing_settings e.g. 
```
MQ23qc
Eau09AI233qc
Eau09AI164qc
PV56sqc
PV70mqc
Eau09AI136qc
PV50mqc
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

To execute the script, make sure you have your phasing_settings file, samples.txt file, reference.fa file, and the Rscripts (modref.R and coverage.R) in the directory where you execute the phasing.sh script. You can do this by:
```
bash phasing.sh
```

Output will include a {samplename}*.1.fa and {samplename}*.2.fa file for each of the samples you listed in samples.txt (containing the first and second allele for each locus, respectively), and a file giving the coverage per sample for each locus ("coverage_summary.txt"). The coverage file will be input into Step 4 to summarize over all the samples.

# Step 4
Summarizing coverage.  These commands will obtain a locus-by-locus summary of coverage and sample size. The output file ("locus_coverage_data.txt") will give the locus name in the first column, how many samples that locus was present in, and then the average coverage across these samples. Make sure summarize_coverage.R is in the folder with your samplename.1.fa and samplename.2.fa files.
```
grep ">" reference.fa | sed 's/>//g' > fasta_names
Rscript summarize_coverage.R
```

# Step 5
Combining reference-aligned samples with original fasta files. If you are combining your data with previously phased data (e.g. fasta files created from a previous *.alleles pyRAD file etc) you'll need to tweak the samplenamesuffix1 and 2 to match what you have in your original fasta files (e.g. changing it between assembled/unassembled etc). Otherwise, if you are rejoining here from the UCE pipepine and have phased every individual in your dataset, then you can change these suffixes to whatever you like (you should also make the combined_fasta directory, but leave it blank as it will be populated by fasta from your samples).

After this step, MAFFT is run to re-align the loci in the files, as our reference-aligned loci will have gaps etc stripped.
```
mkdir combined_fasta
### This step should NOT be performed for UCE pipeline samples where you've phased everyone ###
cp backup/*.fasta combined_fasta

### Tweak these suffixes if desired
samplenamesuffix1=".unassembled_0";
samplenamesuffix2=".unassembled_1";

for i in *.1.fa;
do samplename=`echo $i | sed 's/.1.fa//g'`;
second=`echo $i | sed 's/.1.fa/.2.fa/g'`;
norefs=`wc -l $i | awk '{print $1}'`;
for j in `seq 1 2 $norefs`;
do k=$(($j+1));
locusname=`tail -n+$j $i | head -n1| sed 's/>//g'`;
echo $locusname >> loci_to_align.txt;
echo ">"$samplename$samplenamesuffix1 >> combined_fasta/$locusname; 
firstseq=`tail -n+$k $i | head -n1`;
echo $firstseq >> combined_fasta/$locusname; 
echo ">"$samplename$samplenamesuffix2 >> combined_fasta/$locusname; 
secondseq=`tail -n+$k $second | head -n1`;
echo $secondseq >> combined_fasta/$locusname; 
done;
done;

sort loci_to_align.txt | uniq > uniq_loci_to_align.txt;
rm loci_to_align.txt;

nofasta=`wc -l uniq_loci_to_align.txt | awk '{print $1}'`;
for j in `seq 1 $nofasta`;
do locusname=`tail -n+$j uniq_loci_to_align.txt | head -n1`;
mv combined_fasta/$locusname temp;
mafft temp > combined_fasta/$locusname;
done;

rm temp;
```

# STEP 6
Cleaning up the MAFFT alignments (which wrap text over multiple lines and use lowercase letters). The first sed command is only needed if you have stuffed up the suffixes on Step 5, otherwise you can comment it out and use the cp line instead (if you do need to correct the suffixes, make sure to comment out the cp line instead).
```
nofasta=`wc -l uniq_loci_to_align.txt | awk '{print $1}'`;
for j in `seq 1 $nofasta`;
do locusname=`tail -n+$j uniq_loci_to_align.txt | head -n1`;
#sed 's/\.assembled/\.unassembled/g' combined_fasta/$locusname > temp;
cp combined_fasta/$locusname temp;
rm combined_fasta/$locusname;
Rscript mafft_reformat.R;
mv temp.fa combined_fasta/$locusname;
rm -rf temp;
done;
```

If your downstream analyses (e.g. https://github.com/laninsky/Phase_hybrid_from_next_gen) involve using *.fasta files rather than SNPs, then you are done with this repository! Otherwise proceed to Step 7.

# STEP 7
Pulling out the SNPs from each locus. Create a species_assignments file following the instructions at: https://github.com/laninsky/pyRAD_alleles_into_structure#species_assignments-file
Note, for this species assignments file, do not include suffixes e.g.
```
ChecKK058     hectors
Chem15NZ35    hectors
ChecKK063     hectors
28311 rights
13193 rights
```
Then run the following code
```
Rscript creating_frame.R;

for i in combined_fasta/*.fasta;
do cp $i temp;
echo $i > name.txt;
Rscript fasta_into_SNPs.R;
rm temp;
done;
```
The output will be as follows:

-- locus_summary.txt: a tab-delimited file with locus name, length of locus, number of samples with data for locus, number of samples within each species (defined in your species_assignment) folder, total number of alleles found across your samples, total number of alleles found within each species, total number of bi-allelic SNPs variable across all samples, total number of bi-allelic SNPs observed to be variable within a species, total proportion of heterozygotes in sample, total proportion of heterozygotes within each species, and finally whether sequence is present (1) or absent (0) for each of your samples. For the calculations of numbers of alleles, SNPs and heterozygosity, only samples where data is present are used.

-- full_SNP_record.txt: a tab-delimited file with locus name, position of SNP within locus, and nucleotide base observed at each allele for each sample (therefore sample name given twice). Missing data or gaps in alignment denoted as zero. This file will be the input for further downstream processing of the data to produce a structure file. Samples are grouped by species as defined in your species_assignment file.

-- allele_record.txt: a tab-delimited file with locus name, and haplotype combination across bi-allelic SNPs for each allele for each sample (therefore sample name given twice. Haplotype names arbritrarily numbered starting at 1). Missing data denoted as zero. This file will be the input for further downstream processing of the data to produce a structure file. Samples are grouped by species as defined in your species_assignment file.

-- frame_record.txt: an internal file that the current and next step use to see which individual corresponds to which species.

# STEP 8
So now you have all the SNPs for all the loci, it is time to filter these down and make you a structure file. To decide how to filter your data, this step requires two files: lineage.txt and missing.txt. In lineage.txt, you can define which species (corresponding to your species_assignment file) MUST have data in the structure file (i.e. only SNPs which have at least one individual with data for that species will be selected). If lineage.txt does not exist, then the SNP with the least amount of missing data over all the individuals in the dataset per locus will be selected (defaulting to the first SNP if they all have equal amounts of mising data). The second file, missing.txt, defines which species you want to minimize the amount of missing data in (e.g. if you have ingroup and outgroup species, potentially minimizing the amount of missing data in your ingroup). If no missing.txt file exists, the program defaults to selecting the SNP per locus with the least amount of missing data (taking into account lineage.txt if it exists). The format for each of these files is identical, just each species (in order of preference for the missing.txt file) on a new line e.g.
```
NPrights
NArights
```
If you want to change up the species_asignment file, you'll need to go back and run Step 7 again.

To run step 8:
```
Rscript SNPs_to_structure.R;
```

If you don't care about the coverage of the loci from the samples that you referenced map, you can leave it here (you should have set the cut-off for your loci that came from pyRAD, through pyRAD), otherwise proceed to Step 9. There will be two output files from Step 8:

-- full_SNP_record_step8.txt: same columns etc as 'full_SNP_record.txt' but constrained just to one SNP per locus following filtering for lineage.txt, missing.txt and overall missing data

-- structure_step8.txt: a file containing the same loci listed in full_SNP_record_step8.txt formatted for structure. If you'd like to further filter for missing data, follow https://github.com/laninsky/ambigoos_into_structure#what-if-you-want-to-tweak-the-individuals-in-the-filechange-completeness-of-datasetsnp-selection-criteria

# STEP 9
If you would like to filter for sequencing depth in your reference-aligned samples, you'll need a file with the sequencing depth cut-off you would like to use, named mindepth.txt. In this file, there should just be a single number, the read depth you require for a locus to be included e.g. for 30x sequencing depth
```
30
```
Any of your reference-aligned samples that do not meet this sequencing depth requirement will have their base calls changed to '0'. This step will then run through your lineage.txt file again (if you made one) to make sure all the SNPs still have data for at least one individual for each of the lineages you define here. The other files required for this step are:

-- samples.txt (from Step 3)

-- full_SNP_record_step8.txt (from Step 8)

-- coverage_summary.txt (output from Step 3)

-- frame_record.txt (from Step 8)

To run this step:
```
Rscript filter_for_depth.R
```

If you are happy with the samples in your "ingroup" that you've defined the SNPs over, then you are done!

# STEP 10
If you would like to pare down the samples that your SNPs are defined over (say, you have multiple different species, and you are interested in defining SNPs just over one of those species), you'll want to move all your files from STEP 7 onwards into a backup folder. From that folder, you'll want to copy out full_SNP_record.txt and frame_record.txt. To define the lineages that you want to focus on, you'll need a file called species_filtering.txt. In this file, you'll list the species (corresponding to your species_assignments file) that you would like to restrict the SNP dataset to, e.g.
```
NArights
NPrights
SHrights
```
Then run this step by:
```
mv full_SNP_record.txt tempfull_SNP_record.txt
mv frame_record.txt tempframe_record.txt
Rscript filter_for_species.R
rm tempfull_SNP_record.txt
rm tempframe_record.txt
```

The script will modify your full_SNP_record.txt and frame_record.txt files to only include the taxa of interest. You will then need to run Steps 8-9 again to get your final structure file. Make sure to modify your samples.txt, missing.txt and lineage.txt files for these steps if you have removed taxa that used to be in them.

After defining the SNPs over your ingroups of interest you can then add your outgroups back in by running the following Rscript. You'll need your original full_SNP_record.txt (containing all of the samples) and the full_SNP_record_step9[...] file corresponding to your filtered ingroups SNPs of choice (if you haven't filtered for minimum depth, then just rename your full_SNP_record_step8.txt file to full_SNP_record_step9.txt). Put these in a folder with the Rscript and execute by:
```
Rscript insertingoutgroups.R
```

# STEP 11
If you want to filter down for an amount of missing data in your final file, you'll need a file called 'proportion.txt' with a number - the proportion of missing data you wish to allow - on the first line. On the second line, whether you would like the input structure file that is modified to be one that followed step8 or step9 e.g.
```
0.2
8
```
All you need in the folder is your structurefile_step8.txt or structurefile_step9.txt, proportion.txt and the Rscript, completeness.R.

Then run:
```
Rscript completeness.R
```

# Version history
v0.0.2: GATK v 4.0 is a standalone executable rather than a *.jar file, so this was tweaked in the code on the 26-Mar-2018. The previous phase_everyone.sh files are available as phase_everyone_pre_v4_gatk.sh

v0.0.1: The -stand_emit_conf 30 option is deprecated in GATK v 3.7 and was removed from this code on the 5-June-2017
