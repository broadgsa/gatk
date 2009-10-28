package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.Iterator;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 23, 2009
 *
 * This walker is designed to work as the first pass in a two-pass processing step.
 * It does a by-locus traversal calculating a neighborhood quality score based on a number of factors:
 *  1.) Number of reads in neighborhood whose mate is mapped to a different chromosome.
 *  2.) Average mapping quality of reads in the neighborhood.
 *  3.) Average reference mismatch rate for all reads in the neighborhood.
 * The output file is a list of: (GenomeLoc QualityScore) for every locus
 *
 * This walker is designed to be used in conjunction with ReadQualityScoreWalker.
 */

/**
 * Example of what the output should look like:
 *
 * 1:10999919 25.781164
 * 1:10999920 30.321754
 * 1:10999921 30.321754
 * 1:10999922 30.321754
 * 1:10999923 30.005175
 * 1:10999924 29.82714
 * 1:10999925 29.901012
 * 1:10999926 24.971085
 * 1:10999927 24.634737
 * 1:10999928 21.552652
 * 1:10999929 21.95971
 * 1:10999930 21.95971
 * 1:10999931 20.272423
 * 1:10999932 18.20454
 * 1:10999933 18.20454
 * 1:10999934 18.20454
 */

public class NeighborhoodQualityWalker extends LocusWalker<Integer, Long> {

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        float neighborhoodQualityScore = 0.0f;
        int numReadsMismatchedMate = 0; // num of reads in this locus whose mate is mapped to different chromosome
        float percentMismatchedMate = 0.0f; // percentage of reads at this locus whose mate is mapped to diff chromosome
        float avgMappingQuality = 0.0f; // mean mapping quality for all reads at this locus
        float avgMismatchRate = 0.0f; // mean mismatch with reference rate over all reads at this locus

        int numValidMappingQuality = 0;
        long sumMappingQuality = 0L;
        float sumMismatchRate = 0.0f;
        int mappingQuality = 0; // preallocate for use in while loop below
        boolean isGoodPair = false; // preallocate for use in while loop below
        SAMRecord read = null; // preallocate for use in while loop below

        List<SAMRecord> reads = context.getReads();
        Iterator<SAMRecord> readsIter = reads.iterator();
        assert reads.size() > 0 : "This locus has no reads.";
        
        while( readsIter.hasNext() ) { // for each read in this context
            read = readsIter.next();

            // Only consider reads for this calculation whose mapping quality isn't 0 or 255
            mappingQuality = read.getMappingQuality();
            if ( mappingQuality > 0 && mappingQuality < 255 ) {

                // Generate sum of mapping quality for all reads at this locus
                sumMappingQuality += mappingQuality;
                numValidMappingQuality++;

                // Look to see if mate was mapped to different chromosome
                //isGoodPair = ( read.getReadPairedFlag() ? read.getProperPairFlag() : true );
                isGoodPair = ( !read.getReadPairedFlag() || read.getProperPairFlag() ); // optimized version of above line
                if ( !isGoodPair ) { numReadsMismatchedMate++; }

                // Generate sum number of mismatches for all reads at this locus
                if( read.getAttribute("NM") != null ) {
                    sumMismatchRate += ((float) Integer.parseInt(read.getAttribute("NM").toString())) / ((float) read.getReadLength());
                } else {
                    sumMismatchRate += 1.0f;
                }
            }
        }

        // Calculate averages from sums which accumulated during while loop above
        if ( numValidMappingQuality == 0 ) { numValidMappingQuality = 1; }
        percentMismatchedMate = ((float) numReadsMismatchedMate) / ((float) numValidMappingQuality);
        avgMappingQuality = sumMappingQuality / ((float) numValidMappingQuality);
        avgMismatchRate = sumMismatchRate / ((float) numValidMappingQuality);

        // Calculate the three metrics that go into a neighborhood quality score using exponential decay model
        // BUGBUG: some analysis is needed to determine reasonable rates and scale factors for the exponential functions
        float scoreMates = 40.0f * (float) Math.exp( -16.0f * (float) percentMismatchedMate );
                                    // exp decay with rate 16.0, scaled to Q=40 when mismatched mates is 0%
        float scoreMapping = 40.0f * (float) Math.exp( -0.02f * Math.max( 99.0f - avgMappingQuality, 0.0f ) );
                                    // exp decay with rate 0.02, scaled to Q=40 when avg map quality is 99
        float scoreMismatch = 40.0f * (float) Math.exp( -27.0f * avgMismatchRate );
                                    // exp decay with rate 27.0, scaled to Q=40 when avg mismatch rate in reads is 0

        // BUGBUG: some analysis is needed to determine reasonable weights for each metric
        neighborhoodQualityScore = 0.35f * scoreMates + 0.05f * scoreMapping + 0.6f * scoreMismatch;
        assert neighborhoodQualityScore >= 0.0f : "Neighborhood quality score must be nonnegative.";
        if( neighborhoodQualityScore < 1.0f ) { neighborhoodQualityScore = 1.0f; }

        // verbose debug printing lines
        logger.debug( context.getLocation() + " " + neighborhoodQualityScore );
        logger.debug( "mate mismatch% =\t" + percentMismatchedMate + " --> " + scoreMates );
        logger.debug( "mean mappingQ =\t" + avgMappingQuality + " --> " + scoreMapping );
        logger.debug( "mean mismatchRate =\t" + avgMismatchRate + " --> " + scoreMismatch );

        // This printout useful for making histograms of scores in Matlab
        //out.println( neighborhoodQualityScore + " " + scoreMates + " " + scoreMapping + " " + scoreMismatch);

        out.println( context.getLocation() + " " + neighborhoodQualityScore );
        return 0;
    }

    public Long reduceInit() {
        return 0L;
    }

    public Long reduce( Integer value, Long sum ) {
        return 0L; // nothing to do here
    }

    public void onTraversalDone( Long reduceResult ) {
    }

}