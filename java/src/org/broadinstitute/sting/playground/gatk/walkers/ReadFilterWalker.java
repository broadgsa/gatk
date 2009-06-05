package org.broadinstitute.sting.playground.gatk.walkers;

import java.io.File;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;

public class ReadFilterWalker extends ReadWalker<Integer,Integer> {
	@Argument(fullName="output_file", shortName="O",doc="SAM or BAM file to write filtered reads into (will be overwritten if exists)",required=true ) public String output;
	@Argument(fullName="max_read_length",doc="Discard reads with length greater than the specified value",required=false) public Integer max_len;

	private SAMFileWriter writer = null;
	
    public void initialize() {
        SAMFileHeader header = this.getToolkit().getEngine().getSAMHeader();
        writer = Utils.createSAMFileWriterWithCompression(header, header.getSortOrder() != SAMFileHeader.SortOrder.unsorted, output, getToolkit().getBAMCompression());
    }

	@Override
    public boolean filter(char[] ref, SAMRecord read) {
		return read.getReadLength() <= max_len;
	}

	@Override
	public Integer map(char[] ref, SAMRecord read) {
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
		writer.close();
        out.println(nReads +" reads passed the filter and were written into output file "+output);
	}

}
