package org.broadinstitute.sting.playground.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.*;

import java.io.File;
import java.util.Arrays;

/**
 * User: mdepristo
 *
 * Replaces read groups in a BAM file
 */
public class ReplaceReadGroups extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Creates a new read group, and assigns all reads from the I BAM file to this read group in the O BAM";
    @Option(shortName="I", doc="Input file (bam or sam).", optional=false)
    public File IN = null;
    @Option(shortName="O",doc="Output file (bam or sam).", optional=false)
    public File OUT = null;

    @Option(shortName="ID",doc="Read Group ID", optional=false)
    public String RGID = null;

    @Option(shortName="LB",doc="Read Group Library", optional=false)
    public String RGLB = null;

    @Option(shortName="PL",doc="Read Group platform", optional=false)
    public String RGPL = null;

    @Option(shortName="SM",doc="Read Group sample", optional=false)
    public String RGSM = null;

    private static final String RGFIELD = "RG"; // todo -- use binary tag that's private in picard

    // todo -- is it worth supporting these fields?
    // CN Name of sequencing center producing the read.
    // DS Description.
    // DT Date the run was produced (ISO8601 date or date/time).
    // PU Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identier.

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new ReplaceReadGroups().instanceMain(argv));
    }

    protected int doWork() {
        SAMFileReader inReader = new SAMFileReader(IN);

        // create the read group we'll be using
        SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
        rg.setLibrary(RGLB);
        rg.setPlatform(RGPL);
        rg.setSample(RGSM);
        System.out.printf("Created read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample());

        // create the new header and output file
        SAMFileHeader outHeader = inReader.getFileHeader().clone();
        outHeader.setReadGroups(Arrays.asList(rg));
        SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader, true, OUT) ;

        //
        // write the reads in contig order
        //
        for ( SAMRecord read : inReader ) {
            read.setAttribute(RGFIELD, rg.getId());
            outWriter.addAlignment(read);
        }

        // cleanup
        inReader.close();
        outWriter.close();
        return 0;
    }
}


