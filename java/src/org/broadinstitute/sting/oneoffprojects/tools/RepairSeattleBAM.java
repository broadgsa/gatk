package org.broadinstitute.sting.oneoffprojects.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.samtools.*;

import java.io.*;

public class RepairSeattleBAM extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Fix read group info";
    @Option(shortName="I", doc="Input file (bam or sam) to extract reads from. If not specified, reads from stdin.",
            optional=true) public File IN = null;
    @Option(shortName="O",doc="Output file (bam or sam).",
            optional=true) public File OUT = null;
    @Option(shortName="S",doc="Sample.",
            optional=true) public String SAMPLE = null;

    public static void main(final String[] argv) {
        System.exit(new RepairSeattleBAM().instanceMain(argv));
    }

    protected int doWork() {
        SAMFileReader inReader = new SAMFileReader(IN);

        for (SAMReadGroupRecord rg : inReader.getFileHeader().getReadGroups()) {
            rg.setSample(SAMPLE);
        }

        SAMFileWriter outWriter = (new SAMFileWriterFactory()).makeSAMOrBAMWriter(inReader.getFileHeader(), true, OUT);

        for (SAMRecord read : inReader ) {
            //read.getReadGroup().setSample(SAMPLE);
            read.setAttribute("SM", SAMPLE);

            outWriter.addAlignment(read);
        }

        inReader.close();
        outWriter.close();

        return 0;
    }
}