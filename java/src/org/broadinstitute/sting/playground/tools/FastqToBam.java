package org.broadinstitute.sting.playground.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.samtools.*;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.Iterator;

import org.broadinstitute.sting.utils.StingException;

class FastqRecord {
    private String seqHeader;
    private String seqLine;
    private String qualHeader;
    private String qualLine;

    private String accessionName;
    private String readName;
    private String runName;

    public FastqRecord(BufferedReader in) {
        try {
            if (in.ready()) {
                seqHeader = in.readLine();
                seqLine = in.readLine();
                qualHeader = in.readLine();
                qualLine = in.readLine();

                String[] seqHeaderPieces = seqHeader.split("\\s+");
                accessionName = seqHeaderPieces[0];
                readName = seqHeaderPieces[1];

                String[] readNamePieces = readName.split(":");
                runName = readNamePieces[0];
            }
        } catch (IOException e) {
            throw new StingException("Could not read from fastq file.");
        }
    }

    public String getReadName() { return readName; }
    public String getRunName() { return runName; }
    public String getReadString() { return seqLine; }
    public String getBaseQualityString() { return qualLine; }

    public String toString() {
        return String.format("%s %s : %s : %s", accessionName, readName, seqLine, qualLine);
    }
}

class FastqReader implements Iterator<FastqRecord>, Iterable<FastqRecord> {
    private BufferedReader in;
    private FastqRecord nextRecord;
    private String runName;

    public FastqReader(File file) {
        try {
            in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));

            nextRecord = new FastqRecord(in);
            runName = nextRecord.getRunName();
        } catch (IOException e) {
            throw new StingException("IO problem");
        }
    }

    public String getRunName() { return runName; }

    public boolean hasNext() { return nextRecord != null; }

    public FastqRecord next() {
        FastqRecord rec = nextRecord;

        try {
            if (in.ready()) {
                nextRecord = new FastqRecord(in);
            } else {
                nextRecord = null;
            }
        } catch (IOException e) {
            throw new StingException("IO problem");
        }
        
        return rec;
    }

    public void remove() { throw new UnsupportedOperationException("Unsupported operation"); }

    public Iterator<FastqRecord> iterator() { return this; }
}

public class FastqToBam extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Extracts read sequences and qualities from the input fastq file and writes them into the output file in unaligned BAM format.";

    @Option(shortName="I1", doc="Input file (fastq.gz) to extract reads from (single-end fastq or, if paired, first end of the pair fastq).", optional=false) public File IN1 = null;
    @Option(shortName="I2", doc="Input file (fastq.gz) to extract reads from (if paired, second end of the pair fastq).", optional=true) public File IN2 = null;
    @Option(shortName="O", doc="Output file (bam).", optional=false) public File OUT = null;
    @Option(shortName="RG", doc="Read group name", optional=false) public String READ_GROUP_NAME;
    @Option(shortName="SM", doc="Sample name", optional=false) public String SAMPLE_NAME;
    @Option(shortName="V", doc="Verbose mode", optional=true) public Boolean VERBOSE = false;

    public static void main(final String[] argv) {
        System.exit(new FastqToBam().instanceMain(argv));
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

            if (fqr2 != null && !fqr1.getReadName().equalsIgnoreCase(fqr2.getReadName())) {
                //throw new StingException(String.format("In paired mode, but end 1 read name (%s) does not match end 2 read name (%s)", fqr1.getReadName(), fqr2.getReadName()));
            }

            SAMRecord sr1 = new SAMRecord(sfh);
            sr1.setReadName(fqr1.getReadName());
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

                sr2 = new SAMRecord(sfh);

                sr2.setReadName(fqr2.getReadName());
                sr2.setReadString(fqr2.getReadString());
                sr2.setBaseQualityString(fqr2.getBaseQualityString());
                sr2.setReadUmappedFlag(true);
                sr2.setReadPairedFlag(true);
                sr2.setFirstOfPairFlag(false);
                sr2.setSecondOfPairFlag(false);
                sr2.setMateUnmappedFlag(true);
                sr2.setAttribute("RG", READ_GROUP_NAME);
            }

            sfw.addAlignment(sr1);
            if (fqr2 != null) {
                sfw.addAlignment(sr2);
            }
            readsSeen++;

            if (VERBOSE) {
                //System.out.println(fqr1);
                //if (fqr2 != null) { System.out.println(fqr2); }
                System.out.println(sr1.format());
                if (fqr2 != null) { System.out.println(sr2.format()); }
            }
        }

        sfw.close();

        System.out.println(String.format("%s %s : processed %d reads\n", READ_GROUP_NAME, SAMPLE_NAME, readsSeen));

        return 0;
    }
}
