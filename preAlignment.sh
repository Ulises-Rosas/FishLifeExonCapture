#!/bin/sh

# make a directory to store alignments
mkdir Alignments/


# cat individual exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		contigs=$( grep -c ">" $directory/trinity*.$exon.*filtered.fa )
		if [  $contigs -eq 1  ]
		cat $directory/$exon.final_contigs.fasta >> Alignments/$exon.unaligned.fasta;
		fi;
		fi;
	done;
done < ../FishLifeExonCapture/ExonList.txt

# cat individual mitochondrial exons into files for alignment
while read -r exon;
do
	for directory in *;
	do
		if [  -d $directory  ];
		then
		contigs=$( grep -c ">" $directory/trinity*.$exon.*filtered.fa )
		if [  $contigs -eq 1  ]
		cat $directory/$exon.final_contigs.fasta >> Alignments/$exon.unaligned.fasta;
		fi;
		fi;
	done;
done < ../FishLifeExonCapture/MitochondrialExonList.txt
