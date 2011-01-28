echo "63025520" | awk '{ for(i = 0; i < $1; i += 100000) {print "20:" i+1 "-" (i+100000 < $1 ? i+100000 : $1)}}' > whole_genome_chunked.chr20.hg19.intervals
