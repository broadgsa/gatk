package org.broadinstitute.sting.playground.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.samtools.*;

import java.io.*;
import java.util.List;
import java.util.ArrayList;

public class CGUtilities extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Mark, document me";
    @Option(shortName="I", doc="Input file (bam or sam) to mark primary alignments in.  Must be sorted by name",
                optional=false) public File IN = null;
    @Option(shortName="O",doc="Output file (BAM). If not specified, output is printed to stdout.",
            optional=false) public File OUT = null;
    @Option(shortName="M",doc="Maximum number of reasd to process.",
            optional=true) public int MAX_READS = -1;

    public static void main(final String[] argv) {
        System.exit(new CGUtilities().instanceMain(argv));
    }

    protected int doWork() {
        InputStream ins;
        try {
            ins = new FileInputStream(IN);
        } catch ( FileNotFoundException ie ) {
            System.out.println("Failed to open input file "+IN+": "+ie.getCause());
            return 1;
        }

        SAMFileReader inReader = new SAMFileReader(ins);
        SAMFileHeader header = inReader.getFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUT);

        int i = 0;
        String currentReadName = null;
        List<SAMRecord> possibleAlignments = new ArrayList<SAMRecord>();

        for ( SAMRecord read : inReader ) {
            String readName = read.getReadName();
            //System.out.println("@" + readName);

            if ( currentReadName == null )
                currentReadName = readName;

            if ( ! currentReadName.equals(readName) ) {
                List<SAMRecord> firsts = subsetByPair(possibleAlignments, true);
                List<SAMRecord> seconds = subsetByPair(possibleAlignments, false);

                if ( firsts.size() + seconds.size() != possibleAlignments.size() ) {
                    throw new RuntimeException(String.format("Dropped reads %d + %d != %d at %s%n", firsts.size(), seconds.size(), possibleAlignments.size(), currentReadName));
                }

                if ( firsts.size() > 0 ) out.addAlignment(selectBestAlignment(firsts));
                if ( seconds.size() > 0 ) out.addAlignment(selectBestAlignment(seconds));
                possibleAlignments.clear();
                currentReadName = readName;
            }

            possibleAlignments.add(read);

            //out.addAlignment(read);
            if ( i++ > MAX_READS && MAX_READS != -1 ) break;
        }
        inReader.close();
        out.close();

        return 0;
    }

    protected SAMRecord selectBestAlignment(List<SAMRecord> possibleAlignments) {
        SAMRecord best = null;
        for ( SAMRecord possible : possibleAlignments ) {
            if ( best == null || isBetterAlignment(best, possible) )
                best = possible;
        }

        if ( best != null ) {
            best.setNotPrimaryAlignmentFlag(false); // we're the primary alignment!
        }

        return best;
    }

    protected List<SAMRecord> subsetByPair(List<SAMRecord> reads, boolean getFirst) {
        List<SAMRecord> sub = new ArrayList<SAMRecord>();

        for ( SAMRecord read : reads ) {
            boolean add = read.getFirstOfPairFlag() ? getFirst : ! getFirst;
            if ( ! read.getReadPairedFlag() )
                throw new RuntimeException("Read is unpaired! " + read.format());
            if ( add ) sub.add(read);
        }
        
        return sub;
    }

    protected boolean isBetterAlignment(SAMRecord champ, SAMRecord contender) {
        boolean replace = false;
        if ( contender.getMappingQuality() >= champ.getMappingQuality() ) {
            if ( contender.getMappingQuality() > champ.getMappingQuality() ) {
                replace = true;
            } else if ( contender.getReadLength() > champ.getReadLength() ) {
                replace = true;
            } else {
                ;
                //System.out.printf("Equivalent reads%n  %s%n  %s%n", champ.format(), contender.format());
            }
        }

        //if ( replace )
        //    System.out.printf("Better read %n  %s%n  %s%n", contender.format(), champ.format());

        //System.out.printf("Comparing %s vs. %s => contender is better? %s%n", champ, contender, replace);

        return replace;
    }
}