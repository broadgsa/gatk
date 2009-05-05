package org.broadinstitute.sting.gatk;

import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import static org.junit.Assert.fail;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * User: aaron
 * Date: Apr 28, 2009
 * Time: 5:21:29 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 28, 2009
 * <p/>
 * Class GenomeAnalysisTKTest
 * <p/>
 * A quick move of the unit test cases that were stuffed into GenomeAnalysisTK.
 */
public class GenomeAnalysisTKTest extends BaseTest {

    // our sequence dictionary to check
    private final File refFileName = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/Homo_sapiens_assembly17.fasta");

    // skip n number of chromesomes when we do a chromesome by chromesome jumping
    private final int skipChrome = 15;
    /**
     * This has been blindly moved out of the GenomeAnalysisTK.java where it was hanging on,
     * but the author makes very limited promises of any functionality 
     *
     */
    @Test
    public void testNewReferenceFeatures() {

        final FastaSequenceFile2 refFile = new FastaSequenceFile2(refFileName);
        GenomeLoc.setupRefContigOrdering(refFile);

        List<SAMSequenceRecord> refContigs = refFile.getSequenceDictionary().getSequences();

        /*
        for ( SAMSequenceRecord refContig: refContigs ) {
            System.out.printf("  Traversing from chr1 to %s would require jumping %d bytes%n",
                    refContig.getSequenceName(), refFile.getDistanceBetweenContigs("chr1", refContig.getSequenceName()));
        }
        */
        String lastContig = null;
        List<Double> timings = new ArrayList<Double>();
        for ( SAMSequenceRecord startContig : refFile.getSequenceDictionary().getSequences() ) {
            final String startContigName = startContig.getSequenceName();
            int skip = 1;
            for ( SAMSequenceRecord targetContig : refFile.getSequenceDictionary().getSequences() ) {
                if (skipChrome != skip) {
                    ++skip;
                    continue;
                } else {
                    skip = 1;
                }

                refFile.seekToContig(startContigName, true);
                logger.warn(String.format("Seeking: current=%s, target=%s", startContigName, targetContig.getSequenceName()));
                long lastTime = System.currentTimeMillis();
                final boolean success = refFile.seekToContig(targetContig.getSequenceName(), true);
                long curTime = System.currentTimeMillis();
                final double elapsed = (curTime - lastTime) / 1000.0;
                timings.add(elapsed);
                logger.warn(String.format("  -> Elapsed time %.2f, averaging %.2f sec / seek for %d seeks",
                        elapsed, Utils.averageDouble(timings), timings.size()));

                if ( ! success ) {
                    fail(String.format("Failured to seek to %s from %s", targetContig.getSequenceName(), lastContig ));

                }
                //System.exit(1);
            }
        }

        // code for randomly sampling the seeks
//        Random rnd = new Random();
//        String lastContig = null;
//        List<Double> timings = new ArrayList<Double>();
//        final int N_SAMPLES = 1000;
//        //try { refFile.seekToContig("chr3"); } catch ( IOException e ) {}
//        for ( int i = 0; i < N_SAMPLES; i++ ) {
//            final int nextIndex = rnd.nextInt(refContigs.size());
//            String nextContig = refFile.getSequenceDictionary().getSequence(nextIndex).getSequenceName();
//            //nextContig = "chr2";
//            try {
//                System.out.printf("Seeking: current=%s, target=%s%n", refFile.getContigName(), nextContig);
//                long lastTime = System.currentTimeMillis();
//                final boolean success = refFile.seekToContig(nextContig, true);
//                long curTime = System.currentTimeMillis();
//                final double elapsed = (curTime - lastTime) / 1000.0;
//                timings.add(elapsed);
//                System.out.printf("  -> Elapsed time %.2f, averaging %.2f sec / seek for %d seeks%n",
//                        elapsed, Utils.averageDouble(timings), timings.size());
//
//                if ( ! success ) {
//                    System.out.printf("Failured to seek to %s from %s%n", nextContig, lastContig );
//                }
//                //System.exit(1);
//            } catch ( IOException e ) {
//                System.out.printf("Failured to seek to %s from %s%n", nextContig, lastContig );
//                e.printStackTrace();
//            }
//
//            lastContig = nextContig;
//        }
//        System.exit(1);

/*
        final String targetChr = "chr10";
        try {
            refFile.seekToContig(targetChr);
        } catch ( IOException e ){
            System.out.printf("Failured to seek to %s%n", targetChr);
            e.printStackTrace();
        }
        System.exit(1);


        //List<Double> timings = new ArrayList<Double>();
        final long startTime = System.currentTimeMillis();
        long lastTime = System.currentTimeMillis();

        int i = 0;
        String prevNextContigName = null;
        logger.info(String.format("Walking reference sequence:%n"));
        for ( SAMSequenceRecord refContig: refContigs ) {
            long curTime = System.currentTimeMillis();
            ReferenceSequence contig = refFile.nextSequence();
            final double elapsed = (curTime - lastTime) / 1000.0;
            timings.add(elapsed);
            logger.info(String.format("%2d : expected %s contig, found %s with next of %s after %.2f seconds, average is %.2f", i,
                    refContig.getSequenceName(), contig.getName(), refFile.getNextContigName(), elapsed, Utils.averageDouble(timings)));
            if ( prevNextContigName != null && contig.getName() != null && ! prevNextContigName.equals(contig.getName()) )
                throw new RuntimeIOException(String.format("Unexpected contig ordering %s was expected next, but I found %s?",
                        prevNextContigName, contig.getName()));

            prevNextContigName = refFile.getNextContigName();
            lastTime = curTime;
            i++;

            logger.info(String.format("  Traversing from chr1 to %s would require jumping %d bytes",
                    contig.getName(), refFile.getDistanceBetweenContigs("chr1", contig.getName())));
        }
        */
    }

}
