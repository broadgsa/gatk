package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Pair;
import net.sf.samtools.*;
import edu.mit.broad.picard.reference.ReferenceSequence;

import java.util.ArrayList;
import java.util.Random;
import java.util.HashMap;
import java.io.File;

/**
 * ReadErrorRateWalker assesses the error rate per read position ('cycle') by comparing the
 * read to its home on the reference and noting the mismatch rate.  It ignores reads with
 * indels in them, treats high and low-quality references bases the same, and does not count
 * ambiguous bases as mismatches.  It's also thread-safe, so you can process a slew of reads
 * in short order.
 *
 * @author Kiran Garimella
 */
public class ReplaceQuals extends ReadWalker<SAMRecord, SAMFileWriter> {
    @Argument(shortName="inputQualsBAM",doc="BAM files containing qualities to be replaced",required=true)
    public String inputQualsBAM;

    @Argument(shortName="outputBAM", required=false, doc="output BAM file for reads with replaced quals")
    public String outputFilename = null;

    public int MAX_READS_TO_LOAD = -1;

    private HashMap<String, Pair<SAMRecord, SAMRecord>> readNameToPairs;
    private int READ_PRINT_MOD = 100000;
    private final boolean DEBUG = false;

    public void initialize() {
        readNameToPairs = new HashMap<String, Pair<SAMRecord, SAMRecord>>();

        SAMFileReader samReader = new SAMFileReader(new File(inputQualsBAM));
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        int nReads = 0;
        logger.info("Starting to read inputQualsBAM = " + inputQualsBAM);
        for ( SAMRecord read : samReader ) {
            //System.out.printf("READ is %s%n", read.format());

            final String name = read.getReadName();
            Pair<SAMRecord, SAMRecord> binding = readNameToPairs.containsKey(name) ? readNameToPairs.get(name) : new Pair<SAMRecord, SAMRecord>(null, null);
            if ( read.getFirstOfPairFlag() ) {
                binding.first = read;
            } else {
                binding.second = read;
            }
            readNameToPairs.put(name, binding);

            if ( ++nReads % READ_PRINT_MOD == 0 ) {
                logger.info(String.format("  Read %d reads so far...", nReads));
            }

            if ( nReads > MAX_READS_TO_LOAD && MAX_READS_TO_LOAD != -1 )
                break;
        }

        logger.info("Done reading input BAM");
    }

    /**
     *
     */
    public SAMRecord map(char[] ref, SAMRecord read) {
        final String name = read.getReadName();

        if ( readNameToPairs.containsKey(name) ) {
            Pair<SAMRecord, SAMRecord> binding = readNameToPairs.get(name);
            SAMRecord qRead = read.getFirstOfPairFlag() ? binding.first : binding.second;
            if (qRead != null) {

                if ( DEBUG ) {
                    System.out.printf("Replacing read %s quals with %s%n", read.getReadName(), qRead.getReadName());
                    System.out.printf("%s%n", read.getReadName());
                    System.out.printf("%s%n", qRead.getReadName());
                    System.out.printf("%s%n", read.getReadString());
                    System.out.printf("%s%n", qRead.getReadString());
                    System.out.printf("%s%n", read.getBaseQualityString());
                    System.out.printf("%s%n", qRead.getBaseQualityString());

                    //if (! read.getReadString().equals(qRead.getReadString()))
                    //    throw new RuntimeException(String.format("BUG: equating %s and %s but bases are different", read.getReadName(), qRead.getReadName()));
                }
                
                read.setBaseQualities(qRead.getBaseQualities());
            }
        }
        
        return read;
    }

    // -----------------------------------------------------------------------------------------------
    // Standard i/o reduce
    //
    public void onTraversalDone(SAMFileWriter output) {
        if ( output != null ) {
            output.close();
        }
    }

    public SAMFileWriter reduceInit() {
        if ( outputFilename != null ) { // ! outputBamFile.equals("") ) {
            SAMFileWriterFactory fact = new SAMFileWriterFactory();
            SAMFileHeader header = this.getToolkit().getEngine().getSAMHeader();
            return fact.makeBAMWriter(header, true, new File(outputFilename));
        }
        else {
            return null;
        }
    }

    /**
     *
     *
     */
    public SAMFileWriter reduce(SAMRecord read, SAMFileWriter output) {
        if ( output != null ) {
            output.addAlignment(read);
        } else {
            out.println(read.format());
        }

        return output;
    }
}