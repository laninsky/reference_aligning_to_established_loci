numberofcores=`tail -n+8 phasing_settings | head -n1`;
nosamples=`wc -l samples.txt | awk '{print $1}'`;
delimiter=`tail -n+10 phasing_settings | head -n1`;
locusfield=`tail -n+11 phasing_settings | head -n1`;

# going to run samples in parallel by defining function bysample below
bysample () {
  name=`tail -n+$i samples.txt | head -n1`;
  mkdir $name
  nolines=`wc -l ${name}_pileup.fasta | awk '{print $1}'`
  headercontents=`tail -n+1 ${name}_pileup.fasta | head -n1`

  if [ $locusfield == "second" ]; then
    locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
    headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
    echo $headercontents >> $name/$locus_name.fasta
    fasta_sequence=""

    for j in `seq 2 $nolines`;
      do linecontents=`tail -n+$j ${name}_pileup.fasta | head -n1`;
      if echo $linecontents | grep --quiet ">";
        then echo $fasta_sequence | awk '{print $1}' >> $name/$locus_name.fasta;
        locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
        headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
        echo $headercontents >> $name/$locus_name.fasta
        fasta_sequence=""
      else
        fasta_sequence="$fasta_sequence$linecontents"
      fi
    done;

    echo $fasta_sequence | awk '{print $1}' >> $name/$locus_name.fasta

  else
    locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
    headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
    echo ">"$headercontents >> $name/$locus_name.fasta
    fasta_sequence=""

    for j in `seq 2 $nolines`;
      do linecontents=`tail -n+$j ${name}_pileup.fasta  | head -n1`;
      if echo $linecontents | grep --quiet ">";
        then echo $fasta_sequence | awk '{print $1}' >> $name/$locus_name.fasta;
        locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
        headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
        echo ">"$headercontents >> $name/$locus_name.fasta
        fasta_sequence=""
      else
        fasta_sequence="$fasta_sequence$linecontents"
      fi
    done;
  fi

# Aligning the sequences in the output file
  for j in ${name}/*.fasta;
    do nolines=`wc -l $j | awk '{print $1}'`
    locusname=`echo $j | sed "s'$name/''g" | sed 's/\.fasta//g'`
    echo ">"$locusname >> $name.reference.fa
    if [ "$nolines" -gt "2" ]
      then
      mafft --thread $numberofcores $j > $j.aligned
      cons -sequence $j.aligned -outseq $j.temp -name $locusname -setcase 0
      alignedlines=`wc -l $j.temp | awk '{print $1}'`
      fasta_sequence=""

      for k in `seq 2 $alignedlines`;
        do linecontents=`tail -n+$k $j.temp  | head -n1`;
        fasta_sequence="$fasta_sequence$linecontents";
      done;
      echo $fasta_sequence >> $name.reference.fa

    else
      tail -n+2 $j >> $name.reference.fa
    fi
  done
  
  rm ${name}_pileup.fasta
  rm -rf ${name}
  
}

# for-loop allowing bysample to be parallelized
for i in `seq 1 $nosamples`;
  do bysample "$i" & done
wait

