package org.broadinstitute.sting.playground.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.samtools.*;

import java.io.*;

import org.broadinstitute.sting.utils.fastq.FastqReader;
import org.broadinstitute.sting.utils.fastq.FastqRecord;

public class FastqToBam extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Extracts read sequences and qualities from the input fastq file and writes them into the output file in unaligned BAM format.";

    @Option(shortName="I1", doc="Input file (fastq.gz) to extract reads from (single-end fastq or, if paired, first end of the pair fastq).", optional=false) public File IN1 = null;
    @Option(shortName="I2", doc="Input file (fastq.gz) to extract reads from (if paired, second end of the pair fastq).", optional=true) public File IN2 = null;
    @Option(shortName="O",  doc="Output file (bam).", optional=false) public File OUT = null;
    @Option(shortName="RB", doc="Run barcode", optional=false) public String RUN_BARCODE;
    @Option(shortName="RG", doc="Read group name", optional=false) public String READ_GROUP_NAME;
    @Option(shortName="SM", doc="Sample name", optional=false) public String SAMPLE_NAME;
    @Option(shortName="V",  doc="Verbose mode", optional=true) public Boolean VERBOSE = false;

    public static void main(final String[] argv) {
        System.exit(new FastqToBam().instanceMain(argv));
    }

    private String getReadName(String fqrHeader) {
        String[] headerPieces = fqrHeader.split("\\s+");
        return headerPieces[0];
    }

    protected int doWork() {
        FastqReader end1 = new FastqReader(IN1);
        FastqReader end2 = (IN2 == null) ? null : new FastqReader(IN2);

        SAMReadGroupRecord srg = new SAMReadGroupRecord(READ_GROUP_NAME);
        srg.setSample(SAMPLE_NAME);

        SAMFileHeader sfh = new SAMFileHeader();
        sfh.addReadGroup(srg);
        sfh.setSortOrder(SAMFileHeader.SortOrder.queryname);

        SAMFileWriter sfw = (new SAMFileWriterFactory()).makeSAMOrBAMWriter(sfh, false, OUT);

        int readsSeen = 0;
        while (end1.hasNext() && (end2 == null || end2.hasNext())) {
            FastqRecord fqr1 = end1.next();
            FastqRecord fqr2 = (end2 == null) ? null : end2.next();

            String fqr1Name = getReadName(fqr1.getReadHeader());

            //if (fqr2 != null && !fqr1Name.equalsIgnoreCase(fqr2Name)) {
                //throw new StingException(String.format("In paired mode, but end 1 read name (%s) does not match end 2 read name (%s)", fqr1.getReadName(), fqr2.getReadName()));
            //}

            SAMRecord sr1 = new SAMRecord(sfh);
            sr1.setReadName(RUN_BARCODE + ":" + fqr1Name);
            sr1.setReadString(fqr1.getReadString());
            sr1.setBaseQualityString(fqr1.getBaseQualityString());
            sr1.setReadUmappedFlag(true);
            sr1.setReadPairedFlag(false);
            sr1.setAttribute("RG", READ_GROUP_NAME);

            SAMRecord sr2 = null;

            if (fqr2 != null) {
                sr1.setReadPairedFlag(true);
                sr1.setFirstOfPairFlag(true);
                sr1.setSecondOfPairFlag(false);
                sr1.setMateUnmappedFlag(true);

                String fqr2Name = getReadName(fqr2.getReadHeader());
                sr2 = new SAMRecord(sfh);
                sr2.setReadName(RUN_BARCODE + ":" + fqr2Name);
                sr2.setReadString(fqr2.getReadString());
                sr2.setBaseQualityString(fqr2.getBaseQualityString());
                sr2.setReadUmappedFlag(true);
                sr2.setReadPairedFlag(true);
                sr2.setAttribute("RG", READ_GROUP_NAME);
                sr2.setFirstOfPairFlag(false);
                sr2.setSecondOfPairFlag(true);
                sr2.setMateUnmappedFlag(true);
            }

            sfw.addAlignment(sr1);
            if (fqr2 != null) {
                sfw.addAlignment(sr2);
            }
            readsSeen++;

            if (VERBOSE) {
                System.out.println(sr1.format());
                if (fqr2 != null) { System.out.println(sr2.format()); }
            }
        }

        sfw.close();

        System.out.println(String.format("%s %s : processed %d reads\n", READ_GROUP_NAME, SAMPLE_NAME, readsSeen));

        return 0;
    }
}
