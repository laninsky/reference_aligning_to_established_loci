# Phase_target

You'll have a fasta file with the target loci that you are interested in pulling down from your sample. This might be a file that doesn't have multiple taxa (i.e. each target locus is present just once), or it may have multiple taxa per locus. You'll also have a bunch of samples that you want to fish these loci out of.

You'll need to know the paths to Java, GATK v4.0, GATK v3.8, and picard, and on your path you'll need bwa, samtools, mafft, R, and EMBOSS.

### 1. Using sequences in locus specific fasta file as a reference for genomes of interest
For this step you'll need Java, bwa, samtools and R. If you know the path to the other programs then you might as well chuck them into your 'phasing_settings' file. In addition, you'll need a samples.txt file (both described below), and to have the associated Rscripts in this repository in your working directory as well. 

For phasing_settings, on each separate line, in this order, you will need:

Line 1: path to your up-to-data jdk/bin/java. If your default 'java' command is pointing to the right version, you can just put in java on this line (I have had trouble with $JAVA_HOME on systems where the native shell is not bash)

Line 2: path to folder containing the gatk executable

Line 3: path to picard

Line 4: paired or single, depending on your sequencing method

Line 5: the pathway to your 'cleaned' F reads (or just your cleaned reads if single end). "${name}" should be used as a placeholder for the actual sample name wherever this occurs in the pathway to the reads (as the actual ${name} will be inserted by the script based on samples.txt. This program expects the path to the reads to be standard across your samples, so you will need to rename them if this is not the case

Line 6: the same thing for your 'cleaned' R reads if you have paired sequence data (blank if you don't)

Line 7: the path to GATK3.8/GenomeAnalysisTK.jar (unfortunately not all tools needed were ported to GATK v4)

Line 8: the number of threads you'd like to use PER SAMPLE (e.g. if you have 5 samples, and you select 4 threads, 5\*4 = 20 threads total)

Line 9: the name of the fasta file with the reference loci you want to fish out.

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
target_Zhang_et_al.fa
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
target_Zhang_et_al.fa
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

To run this step, make sure you've installed all the programs, have created the phasing_settings and samples.txt files, and have copied the \*.sh and \*.R scripts to your working directory. Then to run:
```
bash generate_reference.sh &> generate_reference.log
```

You should now have a fasta file for each of your samples. In it the sequences will be named with the original target sequence names

### 2. Generating single consensus sequence for each locus for each sample if multiple samples were present for a given locus in original reference file
Your original fasta file might have had multiple sequences per locus (e.g. a gene characterized across multiple individuals). In this case, we want to make a consensus sequence for each locus for each sample. If each locus was only present once in the target fasta you can skip this step (instead just copy youe pileup.fasta files to reference.fasta files e.g. mv amphizoa_pileup.fasta to amphizoa.reference.fasta). 

You'll need mafft and EMBOSS on your path for this step. You are also going to need to tack on another couple of lines to your phasing_settings file so we can figure out how to parse out the sequences by locus. To do this, you are going to need to know how your locus names are delimited from the whole sequence title e.g. "Gene49" in the following dataset is represented four times in the original target file (by Carabus_sp, Omoglymmius_sp, Pheropsophus_sp, and Agabus_sp)
```
[1] "Up to Carabus_spGene49 for sample bemHap1.pileup"
[1] "Up to Omoglymmius_spGene49 for sample bemHap1.pileup"
[1] "Up to Pheropsophus_spGene49 for sample bemHap1.pileup"
[1] "Up to Agabus_spGene49 for sample bemHap1.pileup"
```

The delimiter that would allow us to identify each of these sequence names as belonging to the same locus is "Gene", taking the second half of the string after it is split on this (e.g. 49).

You could also have the locus identifier come first e.g. GeneA_species1, GeneA_species2, GeneA_species3. In this case "\_" would be the delimiter, and then we'd take the _first_ part of the string after it is split.

One thing just to watch is that the previous scripts will have stripped out any "-" from the names of your sequences. It might be a good idea to check the sample-specific fastas before starting this step to make sure your delimiter is going to work. After you've sussed this out, in phasing_seetings:

Line 10: The delimiter that will allow you to pull out the locus name e.g. Gene, _

Line 11: After splitting on this, whether the locus identified is in the 'first' or 'second' part of the string e.g. first, second

phasing_settings now looks like something this:
```
/usr/bin/java
/usr/local/gatk-4.0.2.1
/usr/local/picard-tools-2.18.0/picard.jar
single
/mnt/hcs-gemmell/beetles/reads/${name}-pe100-reads.fq.gz

/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
2
target_Zhang_et_al.fa
Gene
second
```
To run this step, make sure you've installed all the programs, have created the phasing_settings and samples.txt files, and have copied the \*.sh and \*.R scripts to your working directory. Then to run:
```
bash condense_reference.sh &> condense_reference.log
```

After generating this "final reference" for each sample, we'll use this as the reference to reference map back to once more in Step 3. 

### 3. Carrying out reference mapping against sample-specific references, and phasing each sample for each locus

Using the {name}.reference.fa files generated in the previous step, your phasing_settings and samples.txt files, and the R scripts in this repository, phasing.sh will use gatk v3.8, v4.0, bwa, samtools, and picard to generate phased alignments for each locus for each of your samples/

Output will include a {name}.1.fa and {name}.2.fa file for each of the samples you listed in samples.txt (containing the first and second allele for each locus, respectively), and a file giving the coverage per sample for each locus ("coverage_summary.txt").

To run this step, make sure you've installed all the programs, have created the phasing_settings and samples.txt files, and have copied the \*.sh and \*.R scripts to your working directory. Then to run:
```
bash phasing.sh &> phasing.log
```

### 4. Creating locus-specific fasta files for samples with data for those loci, aligned via mafft
Pulling out the fasta for each locus for each sample:
```
mkdir locus_specific_fasta
for i in *.reference.fa;
  do name=`echo $i | sed 's/.reference.fa//g'`;
  linecontents=`tail -n+1 $i | head -n1`;
  locus_name=`echo $linecontents | sed 's/>//g'`
  echo ">"$name >> locus_specific_fasta/$locus_name.fasta
  nolines=`wc -l $i | awk '{print $1}'`

  for j in `seq 2 $nolines`;
    do linecontents=`tail -n+$j $i | head -n1`;
      if echo $linecontents | grep --quiet ">";
        then echo $fasta_sequence | awk '{print $1}' >> locus_specific_fasta/$locus_name.fasta;
        locus_name=`echo $linecontents | sed 's/>//g'`
        echo ">"$name  >> locus_specific_fasta/$locus_name.fasta
        fasta_sequence=""
      else
        fasta_sequence="$fasta_sequence$linecontents"
      fi
    done;
 done;   
```
Aligning samples for loci with 2+ samples present:
```
  for j in locus_specific_fasta/*.fasta;
    do nolines=`wc -l $j | awk '{print $1}'`
    if [ "$nolines" -gt "2" ]
      then
      mv $j $j.temp
      mafft --thread $numberofcores $j.temp > $j.aligned
      alignedlines=`wc -l $j.aligned | awk '{print $1}'`
      fasta_sequence=""
      linecontents=`tail -n+1 $j.aligned | head -n1`;
      echo $linecontents >> $j
      
      for k in `seq 2 $alignedlines`;
      do linecontents=`tail -n+$k $j.aligned | head -n1`;
        if echo $linecontents | grep --quiet ">";
          then echo $fasta_sequence | awk '{print $1}' >> $j;
          echo $linecontents >> $j;
          fasta_sequence="";
        else
          fasta_sequence="$fasta_sequence$linecontents";
        fi
      done;
    fi
  done
```



