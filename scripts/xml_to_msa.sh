#!/usr/bin/bash
#### This is a bash script for running src/download_genebank.py in command line.
#### The xml is the file of blastn output generated in BLAST site by manual querying.
#### Usage:
####     xml_to_msa.sh A.xml B.xml C.xml ...
#### Output:
####     tempdir/job`date +%Y%m%d`_{A,B,C,...}-<qseq>-genebank.fasta
####     tempdir/job`date +%Y%m%d`_{A,B,C,...}-<qseq>-genebank.aln

echo "virtual env path: $VIRTUAL_ENV"
[ -n "$VIRTUAL_ENV" ] && source "$VIRTUAL_ENV"/bin/activate
[ -z "$JOBID" ] && JOBID="job`date +%Y%m%d`"


help(){
    grep -P "^####" $0 | sed -e "s/####//"
}

get_qseq(){
    in_xml=$1
    grep Hsp_qseq ${in_xml} | head -n1 | sed -e "s/<Hsp_qseq>//g" -e "s/<\/Hsp_qseq>//g" -e "s/ //g" | cut -c 1-30
}

preparing(){
    infile=$@
    for ff in $infile; do
        tarfn=`basename $ff`
        prefix=${tarfn%.xml}
        seq=`get_qseq $ff`
        tardir=tempdir/${JOBID}_${prefix}
        tarfile=${JOBID}_${prefix}-${seq}.xml
        mkdir -p $tardir
        mv $ff ${tardir}/${tarfile}
    done
}

preprocessing(){
## preprocessing
## CREATE_VIEW 
	job_dirs=`find tempdir/${JOBID}_* -type d `
	for ddd in $job_dirs; do
		cd $ddd
			for xml in *.xml; do
				cat $xml | grep -v "CREATE_VIEW" > ${xml}_correction
				mv $xml ${xml}.ori
				mv ${xml}_correction $xml
			done
		cd -
	done
}

run_step1(){
	flank_len=$1
	jobids=`ls -d tempdir/${JOBID}_*` # A directory name itself is the <jobid> argument of `download_genebank.py`
    for jobid in $jobids; do
		[ -f "$jobid" ] && continue # 디렉토리가 아니면 무시
        inxml=`ls ${jobid}/*.xml`
		qseqs=`get_qseq $inxml`
		for qseq in $qseqs; do 
			PYTHONPATH=. python src/download_genebank.py `basename $jobid` $qseq --al_ratio 0.8 --id_ratio 0.9 --max_updown_len $flank_len --final_strand 1
		done
	done
}

### MSA
## MSA 전처리
run_step2(){
	# 이번에는 특별히 없음
	echo "pass run_step2"
}
## MSA 를 한다.
run_step3(){
	for ori_fa in tempdir/${JOBID}_*-genebank.fasta; do
			PYTHONPATH=. python src/multiple_sequence_alignment.py --method clustalo $ori_fa
			PYTHONPATH=. python src/multiple_sequence_alignment.py --method mafft $ori_fa
	done
}

if [ $# == 0 ]; then
    help
    exit 0
fi

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    preparing $@
    preprocessing
	run_step1 500
	if [  $? -ne 0 ]; then 
		exit $?
	fi
	run_step3
	if [  $? -ne 0 ]; then 
		exit $?
	fi
	echo $JOBID
else
	echo "Done: import libraries"
fi

