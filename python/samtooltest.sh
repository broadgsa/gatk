#!/bin/tcsh

setenv refroot ~/work/refmut/
setenv refname test2
setenv ref "${refroot}/${refname}"
setenv cov 100
setenv sites 1
setenv x '' # '-x 6696745:6708910'

# testing samtools
./SimulateReads.py --ref ${ref}.fasta --coverage ${cov} -n ${sites} -l 10 -t None ${x}

# running samtools
samtools import ${refname}__mut.None_10.bp_nsites.${sites}_cov.${cov}.ref ${refname}__mut.None_10.bp_nsites.${sites}_cov.${cov}.sam samtest.bam

samtools sort samtest.bam samtestSorted

samtools index samtestSorted.bam

samtools pileup -f ${ref}.fasta samtestSorted.bam
