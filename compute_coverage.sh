#!/bin/bash

HOME=`pwd`

main_dir=$HOME/scratch/
cd $main_dir

fn_array=()
name_array=()
WSnum_array=()

while IFS=, read -r col1 col2 col3; do
	fn_array+=($col1)
	name_array+=($col2)
	WSnum_array+=($col3)
done < $HOME/config/file_info.txt

coverage () {

	local fn=$1
	local name=$2
	local WSnum=$3

	wdir=$HOME/scratch/$name
	cd $wdir
	
	num=`cut -f2 $name.inserts.xk.uniq.reads.mapped | $HOME/scripts/addCols stdin`

	zcat $name.inserts | awk '{ print ">" $1 "\n" $1 }' | bowtie $HOME/reference/$WSnum/ce_WS230.rna.knownRNA -p 4 -v 0 --best --strata -a -f - > $name.inserts.knownRNA.bowtie

	cut -f1 $name.inserts.knownRNA.bowtie | uniq > $name.inserts.id

	$HOME/scripts/weedLines $name.inserts.id $name.inserts $name.inserts.xk

	awk '{ print ">" $1 "\n" $1 }' $name.inserts.xk | bowtie $HOME/reference/$WSnum/ce.hairpin -p 4 -v 0 --best --strata -a -f - > $name.inserts.xk.hairpin.bowtie
	cut -f1 $name.inserts.xk.hairpin.bowtie | uniq > $name.inserts.xk.id

	$HOME/scripts/weedLines $name.inserts.xk.id $name.inserts.xk $name.inserts.xkxh

	cat $name.inserts.xkxh | awk '{if (length($1) >= 21 && length($1) <= 23){print $0; } }' | awk '{if(substr($1,0,1)=="G"){print $0; } }' | awk '{ print ">" $1 "\n" $1 }' | bowtie $HOME/reference/$WSnum/c_elegans.WS230.genomic -p 4 -v 0 --best --strata -a --un $name.inserts.v0.22G.unMap -f - > tmp
	$HOME/scripts/bowtie_bed.pl tmp > tmp2
	bowtie $HOME/reference/$WSnum/ce_WS230.coding_transcript.juncs -p 4 -v 0 --best --strata -a --un $name.inserts.v0.22G.unMap.unJunc -f $name.inserts.v0.22G.unMap > tmp3
	$HOME/scripts/junc_bowtie_bed.pl tmp3 > tmp4

	cat tmp2 tmp4 | sort -k1,1 -k2,2n > $name.inserts.v0.22G.bed
	intersectBed -a $name.inserts.v0.22G.bed   -b $HOME/reference/$WSnum/ce_WS230.coding_transcript.exon.merge.gene.bed -wo -sorted > bed.tmp
	cat bed.tmp | awk '{if ($6 != $12){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $12;}}' | sort -k1,1 -k2,2n > bed.tmp.anti

	norm=`awk -v nTag=$num 'BEGIN{print 1000000/nTag}'`

	bedtools genomecov -i bed.tmp.anti -g $HOME/reference/$WSnum/c_elegans.WS230.genomic.chromInfo -scale $norm -bg | awk -F t '{ OFS=t; { print $1,$2,$3; }  }'> $name.inserts.v0.22G.norm.anti.bedgraph

	cat $name.inserts.xkxh | awk '{if (length($1) >= 25 && length($1) <= 27){print $0; } }' | awk '{if(substr($1,0,1)=="G"){print $0; } }' | awk '{ print ">" $1 "\n" $1 }' | bowtie $HOME/reference/$WSnum/c_elegans.WS230.genomic -p 4 -v 0 --best --strata -a --un $name.inserts.v0.26G.unMap -f - > tmp
	$HOME/scripts/bowtie_bed.pl tmp > tmp2
	bowtie $HOME/reference/$WSnum/ce_WS230.coding_transcript.juncs -p 4 -v 0 --best --strata -a --un $name.inserts.v0.26G.unMap.unJunc -f $name.inserts.v0.26G.unMap > tmp3
	$HOME/scripts/junc_bowtie_bed.pl tmp3 > tmp4

	cat tmp2 tmp4 | sort -k1,1 -k2,2n > $name.inserts.v0.26G.bed
	intersectBed -a $name.inserts.v0.26G.bed   -b $HOME/reference/$WSnum/ce_WS230.coding_transcript.exon.merge.gene.bed -wo -sorted > bed.tmp
	cat bed.tmp | awk '{if ($6 != $12){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $12;}}' | sort -k1,1 -k2,2n > bed.tmp.anti

	bedtools genomecov -i bed.tmp.anti -g $HOME/reference/$WSnum/c_elegans.WS230.genomic.chromInfo -scale $norm -bg | awk -F t '{ OFS=t; { print $1,$2,$3; }  }'> $name.inserts.v0.26G.norm.anti.bedgraph

	ou_path=$HOME/results
	mkdir $ou_path

	mkdir ${ou_path}/bedgraph
	cp *.bedgraph ${ou_path}/bedgraph

	intersectBed -a $name.inserts.v0.22G.bed   -b $HOME/reference/WS230/ce_WS230.coding_transcript.exon.merge.gene.bed -wo -sorted > bed.tmp
	cat bed.tmp | awk '{if ($6 != $12){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $12;}}' | sort -k1,1 -k2,2n > bed.tmp.anti

	bedtools genomecov -i bed.tmp.anti -g $HOME/reference/WS230/c_elegans.WS230.genomic.chromInfo -scale $norm -d | awk -F t '{ OFS=t; { print $1,$2,$3; }  }'> $name.inserts.v0.22G.norm.anti.depth

	intersectBed -a $name.inserts.v0.26G.bed   -b $HOME/reference/WS230/ce_WS230.coding_transcript.exon.merge.gene.bed -wo -sorted > bed.tmp
	cat bed.tmp | awk '{if ($6 != $12){print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $12;}}' | sort -k1,1 -k2,2n > bed.tmp.anti

	bedtools genomecov -i bed.tmp.anti -g $HOME/reference/WS230/c_elegans.WS230.genomic.chromInfo -scale $norm -d | awk -F t '{ OFS=t; { print $1,$2,$3; }  }'> $name.inserts.v0.26G.norm.anti.depth

	mkdir $HOME/cov/$name
	cp *.depth $HOME/cov/$name/
	cd $HOME/cov/$name

	cp $HOME/scripts/coverage.metagene.sum.rpm.R ./
	cp $HOME/scripts/coverage.metagene.merge.rpm.R ./
	cp $HOME/reference/WS230/ce_WS230.coding_transcript.exon.merge.gene.bed ./

	for i in chrI chrII chrIII chrIV chrV chrX; do
		cat $name.inserts.v0.22G.norm.anti.depth | awk -v ch="$i" '{if ($1 == ch){print $0; } }' > $name.inserts.v0.22G.norm.anti.depth.${i}

		cat $name.inserts.v0.26G.norm.anti.depth | awk -v ch="$i" '{if ($1 == ch){print $0; } }' > $name.inserts.v0.26G.norm.anti.depth.${i}

	done

	wait

	for i in chrI chrII chrIII chrIV chrV chrX; do
		
		for x in 22G 26G; do
			Rscript coverage.metagene.sum.rpm.R $name $x $i
		done

		wait

	done

	wait

	for x in 22G 26G; do
		Rscript coverage.metagene.merge.rpm.R $name $x
	done
	wait

	rm *chr*
	rm *.depth

}

for ((i=0;i<=$((${#fn_array[@]}-1));i++)); do
	coverage ${fn_array[i]} ${name_array[i]} ${WSnum_array[i]}
done
wait

