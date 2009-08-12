package org.broadinstitute.sting.playground.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

import java.io.File;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Aug 11, 2009
 * Time: 5:13:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class FilterReads extends CommandLineProgram {
        @Usage(programVersion="1.0") public String USAGE = "Filters reads: the output file will contain only reads satisfying all the selected criteria";
        @Option(shortName="I", doc="Input file (bam or sam) to extract reads from.",
                optional=false) public File IN = null;
        @Option(shortName="O",doc="Output file (bam or sam) to write extracted reads to.",
            optional=false) public File OUT = null;
        @Option(shortName="U", doc="Select only unmapped reads if true; only mapped reads if false; both if not specified.",
                optional=true) public Boolean UNMAPPED = null;
        @Option(shortName="MINQ", doc="Select only reads with minimum base quality across all bases at or above the specified value.",
            optional=true) public Integer MIN_QUAL = 0;
        @Option(shortName="AVQ", doc="Select only reads with average base quality at or above the specified value.",
            optional=true) public Double AVERAGE_QUAL = 0.0;
        @Option(shortName="MAPQ", doc="Select only reads with mapping quality at or above the specified value (does not affect unmapped reads, use 'U').",
            optional=true) public Integer MAPPING_QUAL = 0;
        @Option(shortName="MAXE",doc="Select only reads with edit distance from the reference at or below the specified value ('NM' tags must be present in the input file).",
            optional = true) public Integer MAX_ERRORS = INFINITY;
        @Option(shortName="MINE",doc="Select only reads with edit distance from the reference at or above the specified value ('NM' tags must be present in the input file).",
            optional = true) public Integer MIN_ERRORS = 0;

        private static int INFINITY = 1000000;
        UnmappedFilter uFilter;

        /** Required main method implementation. */
        public static void main(final String[] argv) {
            System.exit(new FilterReads().instanceMain(argv));
        }

        protected int doWork() {

            if ( UNMAPPED == null ) uFilter = UnmappedFilter.BOTH;
            else {
                if ( UNMAPPED.booleanValue() ) uFilter = UnmappedFilter.UNMAPPED;
                else uFilter = UnmappedFilter.MAPPED;
            }


            SAMFileReader inReader = new SAMFileReader(IN);

            SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(inReader.getFileHeader(), true, OUT) ;

            for ( SAMRecord read : inReader ) {
                switch ( uFilter ) {
                case UNMAPPED: if ( ! read.getReadUnmappedFlag() ) continue; break;
                case MAPPED:   if ( read.getReadUnmappedFlag() ) continue; break;
                }

                if ( ! read.getReadUnmappedFlag() ) {
                    // these filters are applicable only to mapped reads:
                    if ( read.getMappingQuality() < MAPPING_QUAL ) continue;
                    if ( MAX_ERRORS < INFINITY ) {
                        Object attr = read.getAttribute("NM");
                        if ( attr != null ) {
                            int nm = (Integer)attr;
                            if ( nm > MAX_ERRORS ) continue;
                            if ( nm < MIN_ERRORS ) continue;
                        }
                    }
                }


                if ( MIN_QUAL > 0 || AVERAGE_QUAL > 0 ) {
                    byte[] quals = read.getBaseQualities();
                    double av_q = 0.0;
                    boolean passed = true;
                    for ( int i = 0 ; i < quals.length ; i++ ) {
                        if ( quals[i] < MIN_QUAL ) {
                            passed = false;
                            break;
                        }
                        av_q += (double)quals[i];
                    }
                    if ( ! passed ) continue;
                    if ( av_q / read.getReadLength() < AVERAGE_QUAL ) continue;
                }

                outWriter.addAlignment(read);
            }

            inReader.close();
            outWriter.close();
            return 0;
        }

        enum UnmappedFilter {
            UNMAPPED, MAPPED, BOTH
        }

}


