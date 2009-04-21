package org.broadinstitute.sting.playground.gatk.walkers;

import java.io.File;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;

public class ReadFilterWalker extends ReadWalker<Integer,Integer> {
	@Argument(fullName="output_file", shortName="O",doc="SAM or BAM file to write filtered reads into (will be overwritten if exists)",required=true ) public String output;
	@Argument(fullName="max_read_length",doc="Discard reads with length greater than the specified value",required=false) public Integer max_len;


	private SAMFileWriter writer = null;
	
	@Override
    public boolean filter(LocusContext context, SAMRecord read) {
		if ( read.getReadLength() > max_len ) return false;
		return true;
	}

	@Override
	public Integer map(LocusContext context, SAMRecord read) {
		if ( writer == null )	writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(read.getHeader(), read.getHeader().getSortOrder() != SAMFileHeader.SortOrder.unsorted, new File(output));
		writer.addAlignment(read);
		return 1;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		return sum+value;
	}

	@Override
	public Integer reduceInit() {
		return 0;
	}
	
	public void onTraversalDone(Integer nReads) {
		super.onTraversalDone(nReads);
		 out.println(nReads +" reads passed the filter and were written into output file "+output);
	}

}
