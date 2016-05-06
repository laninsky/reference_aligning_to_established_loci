Converting pyRAD *.loci file to a whole ton of fasta files
What you'll need in your folder:

-- your *.loci file from your pyRAD output.

-- a 'species_assignments' file (see https://github.com/laninsky/Phase_hybrid_from_next_gen/blob/master/post-processing/README.md). Note - this file must contain all the species that occur anywhere in your *.loci file, even if they have a lot of mising loci.

-- the converting_pyrad.R script

To convert the file run the following command:
```
for i in `ls *.loci`;
do mv $i temp;
Rscript converting_pyrad.R;
mv temp $i
done;
```
Then you are ready to copy the other scripts you need for the phasing_everyone step to this folder

Version history
14-Oct-15: updated script because it wasn't handling samplenames correctly if they had an 's' in them.
