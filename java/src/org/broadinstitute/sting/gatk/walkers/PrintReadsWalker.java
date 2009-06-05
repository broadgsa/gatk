package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;

import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.io.File;
import java.util.Random;

public class PrintReadsWalker extends ReadWalker<SAMRecord, SAMFileWriter> {

    @Argument(fullName="outputBamFile", shortName="of", doc="Write output to this BAM filename instead of STDOUT", required=false)
    String outputBamFile = null;

    public SAMRecord map(char[] ref, SAMRecord read) {
        return read;
    }

    public SAMFileWriter reduceInit() {
        if ( outputBamFile != null ) { // ! outputBamFile.equals("") ) {
            SAMFileHeader header = this.getToolkit().getEngine().getSAMHeader();
            return Utils.createSAMFileWriterWithCompression(header, true, outputBamFile, getToolkit().getBAMCompression());
        }
        else {
            return null;
        }
    }

    public SAMFileWriter reduce(SAMRecord read, SAMFileWriter output) {
        if ( output != null ) {
            output.addAlignment(read);
        } else {
            out.println(read.format());
        }

        return output;
    }

    public void onTraversalDone(SAMFileWriter output) {
        if ( output != null ) {
            output.close();
        }
    }
}
