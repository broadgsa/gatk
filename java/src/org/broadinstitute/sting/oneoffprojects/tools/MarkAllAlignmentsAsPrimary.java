package org.broadinstitute.sting.oneoffprojects.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

import java.io.*;

public class MarkAllAlignmentsAsPrimary extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Mark all alignments as primary.";
    @Option(shortName="I", doc="Input file (bam or sam) to extract reads from. If not specified, reads from stdin.",
            optional=true) public File IN = null;
    @Option(shortName="O",doc="Output file (bam or sam).",
            optional=true) public File OUT = null;

    public static void main(final String[] argv) {
        System.exit(new MarkAllAlignmentsAsPrimary().instanceMain(argv));
    }

    protected int doWork() {
        SAMFileReader inReader = new SAMFileReader(IN);
        SAMFileWriter outWriter = (new SAMFileWriterFactory()).makeSAMOrBAMWriter(inReader.getFileHeader(), true, OUT);

        for (SAMRecord read : inReader ) {
            read.setNotPrimaryAlignmentFlag(false);

            outWriter.addAlignment(read);
        }

        inReader.close();
        outWriter.close();

        return 0;
    }
}
