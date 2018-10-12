#!/bin/sh


#Filtering with CD-HIT & Exonerate
#For Otophysi markers

for directory in *;
do
	
	if  [  -d $directory  ];
	then
		if [  ! -e $directory.exonfiltering.txt  ];
		then
		echo Filtering exons started $directory > $directory.exonfiltering.txt;
		cd $directory;
			
			for f in *filtered_contigs.fasta;
			do
			# change the similarity threshold here at -c if you want
			cd-hit-est -i $f -o $f.cdhit -c 0.98;
			done;
			
			# Filter nuclear exons with exonerate
			while read -r exon;
			do
			filterfile = $( echo *$exon*.cdhit );
			exonerate --model coding2genome -t $filterfile -q ../../ReadingFramesOtophysi/$exon.fasta --ryo ">%ti\t%qab-%qae\n%tas" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 2 > $filterfile.exonerate.fasta;
			sed -i 's/-- completed exonerate analysis//g' $filterfile.exonerate.fasta;
			done < ../../FishLifeExonCapture/OtophysiExons.txt
			
			for f in *exonerate.fasta;
			do
			python ../../FishLifeExonCapture/filterExons.py -f $f -o $f.exonerate_filtered.fa -t $directory -l 100 -c 1.5;
			done;

			# Filter mitochondrial exons with exonerate
			filterfile = $( echo *G0001*.cdhit );
			exonerate --model coding2genome -t $filterfile -q ../../ReadingFramesOtophysi/G0001.fasta --ryo ">%ti\t%qab-%qae\n%tas" --geneticcode 2 --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F --bestn 2 > $filterfile.exonerateMito.fasta;
			sed -i 's/-- completed exonerate analysis//g' $filterfile.exonerateMito.fasta;

			for f in *exonerateMito.fasta;
			do
			python ../../FishLifeExonCapture/filterExons.py -f $f -o $f.exonerate_filtered.fa -t $directory -l 100 -c 1.5 -m True;
			done;

			
		cd ../;
		echo Filtering exons completed $directory > $directory.exonfiltering.txt;	
		fi;
	fi;	
done	