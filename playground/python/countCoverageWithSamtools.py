import farm_commands
import os.path
import sys

for line in open(sys.argv[1]):
    fastb = line.strip()
    head, fastb_filename = os.path.split(fastb)
    filebase = os.path.splitext(fastb_filename)[0]
    fasta = filebase + '.fasta'

    # convert the fasta
    if not os.path.exists(fasta):
        cmd = "Fastb2Fasta IN="+fastb+" OUT="+fasta
        farm_commands.cmd(cmd)

    qualb = os.path.join(head, filebase + '.new.qualb')
    quala = filebase + '.quala'
    if not os.path.exists(quala):
        cmd = "Qualb2Quala IN="+qualb+" OUT="+quala
        farm_commands.cmd(cmd)
                         
    fastq = filebase + '.fastq'
    if not os.path.exists(fastq):
        cmd = "FastaQuals2Fastq.py "+fasta+" "+quala+ " "+fastq
        farm_commands.cmd(cmd)

    filteredFastq = filebase + '.filtered.fastq'
    if not os.path.exists(filteredFastq):
        cmd = "/seq/dirseq/maq-0.7.1/maq catfilter "+fastq+" > "+filteredFastq
        farm_commands.cmd(cmd)

    filteredFasta = filebase + '.filtered.fasta'
    print 'Looping'
    if not os.path.exists(filteredFasta):
        out = open(filteredFasta,'w')
        iter = open(filteredFastq).__iter__();
        for line in iter:
            if line[0] == '@':
                print >> out, '>%s' % line[1:].strip()
                print >> out, iter.next().strip()
        
    sam = filebase + '.sam'
    cmd = "bwahuman samse 32 " + fastq + " " + sam
    print cmd

samtools view tcga-freeze3-tumor.rev_2.bam chr1:15,000,000-15,002,000 | wc --lines