#!/bin/bash

HOME=`pwd`

main_dir=$HOME/scratch/
mkdir $main_dir
cd $main_dir

fn_array=()
name_array=()
WSnum_array=()

while IFS=, read -r col1 col2 col3; do
	fn_array+=($col1)
	name_array+=($col2)
	WSnum_array+=($col3)
done < $HOME/config/file_info.txt

pipeline () {

	local fn=$1
	local name=$2
	local WSnum=$3

	wdir=$HOME/scratch/$name
	mkdir $wdir
	cd $wdir

	cp $HOME/scripts/pipeline.ce.quantify.tu.R ./
	cp $HOME/reference/geneIDs.WS230 ./

	# trim barcode and 3' linker
	zcat $fn | $HOME/scripts/Extract_insert_6mer_1mm.pl - AGATCGGAAGAGCACACGTCT | cut -f1 | awk '{ if(length($1)>16) print}' > $name.inserts

	$HOME/scripts/inserts2uniqreads.pl $name.inserts 20 70 > $name.inserts.uniq.reads

	# map to known rna and remove them
	awk '{ print ">" $1 "\n" $1 }' $name.inserts.uniq.reads | bowtie $HOME/reference/$WSnum/ce_WS230.rna.knownRNA -p 4 -v 0 --best --strata -a -f - > $name.inserts.knownRNA.bowtie
	cut -f1 $name.inserts.knownRNA.bowtie | uniq > $name.inserts.id

	$HOME/scripts/weedLines $name.inserts.id $name.inserts.uniq.reads $name.inserts.xk.uniq.reads
	$HOME/scripts/weedLines -invert $name.inserts.id $name.inserts.uniq.reads $name.inserts.knownRNA.uniq.reads 

	$HOME/scripts/bowtie_all2ntm_mod.pl $name.inserts.uniq.reads $name.inserts.knownRNA.bowtie > $name.inserts.knownRNA.uniq.reads.v0.ntm

	# map to hairpin sequences and remove them
	awk '{ print ">" $1 "\n" $1 }' $name.inserts.xk.uniq.reads | bowtie $HOME/reference/$WSnum/ce.hairpin -p 4 -v 0 --best --strata -a -f - > $name.inserts.xk.hairpin.bowtie
	cut -f1 $name.inserts.xk.hairpin.bowtie | uniq > $name.inserts.xk.id

	$HOME/scripts/weedLines $name.inserts.xk.id $name.inserts.xk.uniq.reads $name.inserts.xkxh.uniq.reads
	$HOME/scripts/weedLines -invert $name.inserts.xk.id $name.inserts.xk.uniq.reads $name.inserts.xk.hairpin.uniq.reads 

	$HOME/scripts/bowtie_all2ntm_mod.pl $name.inserts.xk.uniq.reads $name.inserts.xk.hairpin.bowtie > $name.inserts.xk.hairpin.uniq.reads.v0.ntm

	# map to genome
	awk '{ print ">" $1 "\n" $1 }' $name.inserts.xk.uniq.reads | bowtie $HOME/reference/$WSnum/c_elegans.WS230.genomic -p 4 -v 0 --best --strata -a --un $name.inserts.xk.uniq.reads.v0.unMap -f - > $name.inserts.xk.uniq.reads.v0.bowtie
	awk '{ print ">" $1 "\n" $1 }' $name.inserts.xkxh.uniq.reads | bowtie $HOME/reference/$WSnum/c_elegans.WS230.genomic -p 4 -v 0 --best --strata -a --un $name.inserts.xkxh.uniq.reads.v0.unMap -f - > $name.inserts.xkxh.uniq.reads.v0.bowtie
	 
	# map to splice junctions
	bowtie $HOME/reference/$WSnum/ce_WS230.coding_transcript.juncs -p 4 -v 0 --best --strata -a --un $name.inserts.xk.uniq.reads.v0.unMap.unJuncs -f $name.inserts.xk.uniq.reads.v0.unMap > $name.inserts.xk.uniq.reads.juncs.v0.bowtie
	bowtie $HOME/reference/$WSnum/ce_WS230.coding_transcript.juncs -p 4 -v 0 --best --strata -a --un $name.inserts.xkxh.uniq.reads.v0.unMap.unJuncs -f $name.inserts.xkxh.uniq.reads.v0.unMap > $name.inserts.xkxh.uniq.reads.juncs.v0.bowtie

	# convert all bowtie files to ntm
	$HOME/scripts/bowtie_all2ntm_mod.pl $name.inserts.xk.uniq.reads   $name.inserts.xk.uniq.reads.v0.bowtie   > $name.inserts.xk.uniq.reads.v0.ntm
	$HOME/scripts/bowtie_all2ntm_mod.pl $name.inserts.xkxh.uniq.reads $name.inserts.xkxh.uniq.reads.v0.bowtie > $name.inserts.xkxh.uniq.reads.v0.ntm

	$HOME/scripts/junc_bowtie_all2ntm_mod.pl $name.inserts.xk.uniq.reads   $name.inserts.xk.uniq.reads.juncs.v0.bowtie   > $name.inserts.xk.uniq.reads.juncs.v0.ntm
	$HOME/scripts/junc_bowtie_all2ntm_mod.pl $name.inserts.xkxh.uniq.reads $name.inserts.xkxh.uniq.reads.juncs.v0.bowtie > $name.inserts.xkxh.uniq.reads.juncs.v0.ntm

	cat $name.inserts.xk.uniq.reads.v0.ntm   $name.inserts.xk.uniq.reads.juncs.v0.ntm   | cut -f4 | $HOME/scripts/weedLines stdin $name.inserts.xk.uniq.reads -invert $name.inserts.xk.uniq.reads.mapped
	cat $name.inserts.xkxh.uniq.reads.v0.ntm $name.inserts.xkxh.uniq.reads.juncs.v0.ntm | cut -f4 | $HOME/scripts/weedLines stdin $name.inserts.xkxh.uniq.reads -invert $name.inserts.xkxh.uniq.reads.mapped
	 
	# variables for R quantification
	num=`cut -f2 $name.inserts.xk.uniq.reads.mapped | $HOME/scripts/addCols stdin`
	num2=$num   ## `addCols $name.inserts.xkxh.uniq.reads.mapped | awk '{ print $2 }'`
	 
	# quantify the 21U (type-1)
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.ntm -b $HOME/reference/$WSnum/ce_WS230.ncrna.21urn.bed -f 1 -r -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.21u $num cosmid WS230

	# quantify the 21U (type-2)
	cat $HOME/reference/$WSnum/ce_WS230.ncrna.21urn.type2.bed | awk -F "	" '{ OFS=","; print $4,$4,$4; }' | sort -k1,1 -u >tmp.geneID.21uType2
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.ntm -b $HOME/reference/$WSnum/ce_WS230.ncrna.21urn.type2.bed -wo -f 0.90 -r | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.21u.type2 $num cosmid WS230 tmp.geneID.21uType2
	 
	## quantify the hairpin
	intersectBed -a $name.inserts.xk.uniq.reads.v0.ntm -b $HOME/reference/$WSnum/ce_WS230.ncrna.hairpin.mirbase.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp
	Rscript pipeline.ce.quantify.tu.R tmp $name.inserts.xk.uniq.reads.v0.hairpin $num cosmid WS230

	## quantify the miRNA
	intersectBed -a $name.inserts.xk.uniq.reads.v0.ntm -b $HOME/reference/$WSnum/ce_WS230.ncrna.mirna.mirbase.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp
	Rscript pipeline.ce.quantify.tu.R tmp $name.inserts.xk.uniq.reads.v0.mirna $num cosmid WS230

	## quantify the mRNA
	cat $name.inserts.xkxh.uniq.reads.v0.ntm $name.inserts.xkxh.uniq.reads.juncs.v0.ntm | intersectBed -a stdin -b $HOME/reference/$WSnum/ce_WS230.coding_transcript.exon.merge.gene.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.mrna $num2 WBGene WS230

	## quantify the pseudogene
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.ntm -b $HOME/reference/$WSnum/ce_WS230.pseudogene.exon.merge.gene.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.pseudogene.exon $num2 WBGene WS230

	## quantify the repeat
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.ntm -b $HOME/reference/$WSnum/ce_WS230.rmsk.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.rmsk $num2 WBGene WS230

	## quantify the 22G mRNA
	cat $name.inserts.xkxh.uniq.reads.v0.ntm | awk -F "\t" '{ OFS="\t"; if(length($4) >= 21 && length($4) <= 23){print $0; } }' | awk '{if(substr($4,0,1)=="G"){print $0; } }' >$name.inserts.xkxh.uniq.reads.v0.22G.ntm
	cat $name.inserts.xkxh.uniq.reads.juncs.v0.ntm | awk -F "\t" '{ OFS="\t"; if(length($4) >= 21 && length($4) <= 23){print $0; } }' | awk '{if(substr($4,0,1)=="G"){print $0; } }' >$name.inserts.xkxh.uniq.reads.juncs.v0.22G.ntm

	cat $name.inserts.xkxh.uniq.reads.v0.22G.ntm $name.inserts.xkxh.uniq.reads.juncs.v0.22G.ntm | intersectBed -a stdin -b $HOME/reference/$WSnum/ce_WS230.coding_transcript.exon.merge.gene.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript /gpfs/data/hlee-lab/lib/pipeline.ce.quantify.gfp.R tmp2 $name.inserts.xkxh.uniq.reads.v0.22G.mrna $num2 WBGene WS230

	## quantify the 22G pseudogene
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.22G.ntm -b $HOME/reference/$WSnum/ce_WS230.pseudogene.exon.merge.gene.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.22G.pseudogene.exon $num2 WBGene WS230

	## quantify the 22G repeat
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.22G.ntm -b $HOME/reference/$WSnum/ce_WS230.rmsk.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.22G.rmsk $num2 WBGene WS230

	## quantify the 26G mRNA
	cat $name.inserts.xkxh.uniq.reads.v0.ntm | awk -F "\t" '{ OFS="\t"; if(length($4) >= 25 && length($4) <= 27){print $0; } }' | awk '{if(substr($4,0,1)=="G"){print $0; } }' >$name.inserts.xkxh.uniq.reads.v0.26G.ntm
	cat $name.inserts.xkxh.uniq.reads.juncs.v0.ntm | awk -F "\t" '{ OFS="\t"; if(length($4) >= 25 && length($4) <= 27){print $0; } }' | awk '{if(substr($4,0,1)=="G"){print $0; } }' >$name.inserts.xkxh.uniq.reads.juncs.v0.26G.ntm

	cat $name.inserts.xkxh.uniq.reads.v0.26G.ntm $name.inserts.xkxh.uniq.reads.juncs.v0.26G.ntm | intersectBed -a stdin -b $HOME/reference/$WSnum/ce_WS230.coding_transcript.exon.merge.gene.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.26G.mrna $num2 WBGene WS230


	## quantify the 26G pseudogene
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.26G.ntm -b $HOME/reference/$WSnum/ce_WS230.pseudogene.exon.merge.gene.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.26G.pseudogene.exon $num2 WBGene WS230

	## quantify the 26G repeat
	intersectBed -a $name.inserts.xkxh.uniq.reads.v0.26G.ntm -b $HOME/reference/$WSnum/ce_WS230.rmsk.bed -wo | $HOME/scripts/intersectBed2Best.pl > tmp2
	Rscript pipeline.ce.quantify.tu.R tmp2 $name.inserts.xkxh.uniq.reads.v0.26G.rmsk $num2 WBGene WS230


	ou_path=$HOME/results
	mkdir $ou_path

	gzip $name.inserts
	gzip $name.inserts.uniq.reads
	gzip $name.inserts.xk.uniq.reads
	gzip $name.inserts.xkxh.uniq.reads
	gzip $name.inserts.knownRNA.uniq.reads
	gzip $name.inserts.hairpin.uniq.reads
	gzip *.unMap*	

	mkdir $ou_path/$name
	mkdir $ou_path/$name/rpm

	cp $wdir/$name.inserts.gz     $ou_path/$name
	cp $wdir/$name.inserts.uniq.reads.gz     $ou_path/$name
	cp $wdir/$name.inserts.xk.uniq.reads.gz  $ou_path/$name
	cp $wdir/$name.inserts.xkxh.uniq.reads.gz  $ou_path/$name
	cp $wdir/$name.inserts.knownRNA.uniq.reads.gz  $ou_path/$name
	cp $wdir/$name.inserts.xk.hairpin.uniq.reads.gz   $ou_path/$name
	cp $wdir/$name.inserts.xk.uniq.reads.mapped $ou_path/$name
	cp $wdir/$name.inserts.xkxh.uniq.reads.mapped $ou_path/$name
	cp $wdir/*v0*.ntm   $ou_path/$name
	cp $wdir/*.ppm   $ou_path/$name/rpm
	cp $wdir/*.unMap* $ou_path/$name


}

for ((i=0;i<=$((${#fn_array[@]}-1));i++)); do
	pipeline ${fn_array[i]} ${name_array[i]} ${WSnum_array[i]}
done
wait

mkdir $HOME/results/master

# mirna master file
Rscript merge.mirna.R ${name_array[@]} v0.mirna.sense

# 21U master file
Rscript merge.21rmsk.R ${name_array[@]} v0.21u.sense
Rscript merge.21rmsk.R ${name_array[@]} v0.21u.type2.sense

# 22G master file
Rscript merge.R ${name_array[@]} v0.22G.mrna.anti
Rscript merge.21rmsk.R ${name_array[@]} v0.22G.rmsk.anti

# 26G master file
Rscript merge.R ${name_array[@]} v0.26G.mrna.anti

Rscript merge.R ${name_array[@]} v0.mrna.sense
Rscript merge.R ${name_array[@]} v0.mrna.anti

