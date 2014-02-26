/*
* Copyright (c) 2012 The Broad Institute
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.sting.gatk.walkers.readutils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMFileWriter;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
/**
 * Utility tool to blindly strip base adaptors. Main application is for FASTQ/unaligned BAM pre-processing where libraries
 * have very short inserts, and hence a substantial part of the sequencing data will have adaptor sequence present.
 * <p>
 * By design, tool will only work for Illumina-like library constructs, where the typical library architecture is:
 * [Adaptor 1]-[Genomic Insert]-[Adaptor 2 (index/barcode)]
 * <p>
 * It is assumed that when data is paired, one read will span the forward strand and one read will span the reverse strand.
 * Hence, when specifying adaptors they should be specified as both forward and reverse-complement to make sure they're removed in all cases.
 * By design, as well, "circular" constructions where a read can have an insert, then adaptor, then more genomic insert, are not supported.
 * When an adaptor is detected, all bases downstream from it (i.e. in the 3' direction) will be removed.
 * Adaptor detection is carried out by looking for overlaps between forward and reverse reads in a pair.
 * If a sufficiently high overlap is found, the insert size is computed and if insert size < read lengths adaptor bases are removed from reads.
 *
 * Advantages over ReadClipper:
 * - No previous knowledge of adaptors or library structure is necessary
 *
 * Advantages over 3rd party tools like SeqPrep:
 * - Can do BAM streaming instead of having to convert to fastq
 * - No need to merge reads - merging reads can have some advantages, but complicates downstream processing and loses information that can be used,
 *   e.g. in variant calling
 * <p>
 *
 * <h2>Input</h2>
 * <p>
 * The input read data in BAM format. Read data MUST be in query name ordering as produced, for example with Picard's FastqToBam
 *
 * <h2>Output</h2>
 * <p>
 * A merged BAM file with unaligned reads
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T ReadAdaptorTrimmer \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -o trimmed_Reads.bam
 * </pre>
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.READ)
public class ReadAdaptorTrimmer extends ReadWalker<List<GATKSAMRecord>, SAMFileWriter> implements NanoSchedulable {
    @Output(doc="Write output to this BAM filename instead of STDOUT", required = false)
    SAMFileWriter out;

     /**
     * Only prints the first n reads of the file - for short testing
     */
     @Hidden
    @Argument(fullName = "number", shortName = "n", doc="Print the first n reads from the file, discarding the rest", required = false)
    int nReadsToPrint = -1;

    /**
     * Argument to control strictness of match between forward and reverse reads - by default, we require 15 matches between them to declare
     * an overlap.
     */
    @Advanced
    @Argument(fullName = "minMatches", shortName = "minMatches", doc="Minimum number of substring matches to detect pair overlaps", required = false)
    int minMatchesForOverlap = 15;


    /**
     * If true, this argument will make the walker discard unpaired reads instead of erroring out.
     */
    @Advanced
    @Argument(fullName = "removeUnpairedReads", shortName = "removeUnpairedReads", doc="Remove unpaired reads instead of erroring out", required = false)
    boolean cleanUnpairedReads = false;

     /**
     * private class members
     */
    private GATKSAMRecord firstReadInPair;
    private TrimStats trimStats = new TrimStats();

    static class TrimStats {
        long numReadsProcessed;
        long numReadsWithAdaptorTrimmed;
        long numUnpairedReadsFound;
    }

   /**
     * The reads filter function.
     *
     * @param ref  the reference bases that correspond to our read, if a reference was provided
     * @param read the read itself, as a GATKSAMRecord
     * @return true if the read passes the filter, false if it doesn't
     */
    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
         // check if we've reached the output limit
        if ( nReadsToPrint == 0 ) {
            return false;          // n == 0 means we've printed all we needed.
        }
        else if (nReadsToPrint > 0) {
            nReadsToPrint--;       // n > 0 means there are still reads to be printed.
        }
        return true;
    }
    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     *
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    public SAMFileWriter reduceInit() {
        return out;
    }

    public List<GATKSAMRecord> map( final ReferenceContext ref, final GATKSAMRecord readIn, final RefMetaDataTracker metaDataTracker ) {


        final List<GATKSAMRecord> readsToEmit = new ArrayList<GATKSAMRecord>();


        // cache first read in pair if flag set.
        if (readIn.getFirstOfPairFlag()) {
            firstReadInPair = GATKSAMRecord.emptyRead(readIn);
            firstReadInPair.setReadString(readIn.getReadString());
            firstReadInPair.setReadName(readIn.getReadName());
            firstReadInPair.setBaseQualities(readIn.getBaseQualities());
        }
        else {
            if (!readIn.getReadName().equals(firstReadInPair.getReadName()))  {
                if (cleanUnpairedReads) {
                    trimStats.numUnpairedReadsFound++;
                    return readsToEmit;
                }
                else // by default require that reads be completely paired
                    throw new IllegalStateException("Second read in pair must follow first read in pair: data not ordered?");
            }

            final int oldLength1 = firstReadInPair.getReadLength();
            final int oldLength2 = readIn.getReadLength();
            // try to strip any adaptor sequence in read pair
            final Integer result = trimReads(firstReadInPair, readIn, minMatchesForOverlap, logger);

            if (logger.isDebugEnabled()) {
                if (result == null)
                    logger.debug("No overlap found, insert size cannot be computed");
                else
                    logger.debug("Insert size estimate = " + result);

            }


            readsToEmit.add(firstReadInPair);
            readsToEmit.add(readIn);

            if (oldLength1 != firstReadInPair.getReadLength())
                trimStats.numReadsWithAdaptorTrimmed++;
            if (oldLength2 != readIn.getReadLength())
                trimStats.numReadsWithAdaptorTrimmed++;

         }


        trimStats.numReadsProcessed++;
        return readsToEmit;

    }

    /**
     * given a read and a output location, reduce by emitting the read
     *
     * @param readsToEmit   the read itself
     * @param output the output source
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public SAMFileWriter reduce( final List<GATKSAMRecord> readsToEmit, final SAMFileWriter output ) {
        for (final GATKSAMRecord read : readsToEmit)
             output.addAlignment(read);

        return output;
    }

    @Override
    public void onTraversalDone(SAMFileWriter output) {

        logger.info("Finished Trimming:");
        logger.info("Number of processed reads:                     "+ trimStats.numReadsProcessed);
        logger.info("Number of reads with adaptor sequence trimmed: "+ trimStats.numReadsWithAdaptorTrimmed);
        if (cleanUnpairedReads)
            logger.info("Number of unpaired reads thrown out: "+ trimStats.numUnpairedReadsFound);
    }


    /**
     *
     * Workhorse routines...
     *
     */
        /**
         * Core routine that does most underlying work for walker. Takes two reads and looks for overlaps in them.
         * An overlap is defined as a contiguous chunk of N bases that matches reverse-complement between reads.
         * Currently, the only insert structure that it will look for overlaps is as follows:
         * CASE 1: Insert shorter than read length:
         * 3' XXXXXXXXXXXXXXXX 5'            (second read)
         * 5'      YYYYYYYYYYYYYYYY 3'       (first read)
         *         ***********
         *
         * In this case, if X and Y are complements at the 11 positions marked by *, routine will do the following
         * iff minMatchesForOverlap <= 11:
         *  a) Cleave adaptor from end of second read (leftmost dangling part in diagram above)
         *  b) Cleave adaptor from end of first read (rightmost part in diagram).
         *
         * CASE 2: Insert size >= read length:
         * 3'             XXXXXXXXXXXXXXXX 5'           (second read)
         * 5'      YYYYYYYYYYYYYYYY 3'                  (first read)
         *                *********                        (overlap)
         *
         * In this case, no trimming is done and reads are left unchanged
         * @param first                      (I/O) First read in pair - read contents (bases/quals) can be modified if adaptor is detected
         * @param second                     (I/O) Second read in pair - read contents (bases/quals) can be modified if adaptor is detected
         * @param minMatchesForOverlap       Reads need to match in these # of bases to be joined
         * @return                           Offset between second and first read.
         *                                   If there's no detectable offset, return Null
         */
    @Requires({"first != null","second != null","minMatchesForOverlap>0"})
    protected static Integer trimReads(final GATKSAMRecord first,
                                       final GATKSAMRecord second,
                                       final int minMatchesForOverlap,
                                       final Logger logger) {

        final Integer insertSize = estimateInsertSize(first.getReadBases(), second.getReadBases(),
                minMatchesForOverlap, logger);

        if (insertSize == null)
            return insertSize;
        if (insertSize < first.getReadLength()) {
            // trim adaptor sequence from read
            first.setReadBases(Arrays.copyOfRange(first.getReadBases(),0,insertSize));
            first.setBaseQualities(Arrays.copyOfRange(first.getBaseQualities(),0,insertSize));
        }
        if (insertSize < second.getReadLength()) {
            // trim adaptor sequence from read
            second.setReadBases(Arrays.copyOfRange(second.getReadBases(),0,insertSize));
            second.setBaseQualities(Arrays.copyOfRange(second.getBaseQualities(),0,insertSize));
        }
        return insertSize;
    }

    /**
    * Brain-dead implementation of an aligner of two sequences, where it's assumed that there might be an overlap
    * from the first into the second. From this, an estimate of insert size is performed and returned
    * Assumes that reads come in reverse direction, so one of the base sequences needs to be reverse-complemented.]
    *
    * @param firstRead                           Bytes from first read
    * @param secondRead                          Bytes from second read (reverse direction)
    * @return                                  Estimated insert size based on offset between first and second read.
    *                                          If no overlap can be detected, return null
    */

    @Requires({"firstRead != null","secondRead != null","minMatches>0","firstRead.length == secondRead.length"})
    protected static Integer estimateInsertSize(final byte[] firstRead,
                                                                final byte[] secondRead,
                                                                final int minMatches,
                                                                final Logger logger) {
        final byte[] firstBases = firstRead;
        final byte[] secondBases = BaseUtils.simpleReverseComplement(secondRead);

        final Pair<Integer,Integer> overlaps = findOverlappingSequence(firstBases, secondBases);
        final int bestOffset = overlaps.first;
        final int maxScore = overlaps.second;
        if ( logger.isDebugEnabled()) {
            String sb="", s1 = new String(firstBases), s2 = new String(secondBases);
            for (int k=0; k < Math.abs(bestOffset); k++) sb+=" ";
            if (maxScore >= minMatches) {
                logger.debug(String.format("Match, Max Score = %d, best offset = %d\n",maxScore, bestOffset));
                if (bestOffset>0)
                    s2 = sb+s2;
                else
                    s1 = sb+s1;
            }
            else logger.debug("NoMatch:");
            logger.debug("R1:"+s1);
            logger.debug("R2:"+s2);


        }

        if (maxScore < minMatches)
            return null; // no overlap detected

        return bestOffset+secondRead.length;


    }


     /**
     * Tries to find overlapping sequence between two reads, and computes offset between them
      * For each possible offset, computes matching score, which is = MATCH_SCORE*Num_matches + MISMATCH_SCORE*num_mismatches
      * (like SW with infinite gap penalties).
     * @param first                              First read bytes
     * @param second                             Second read bytes
     * @return                                   Pair of integers (x,y). x = best offset between reads, y = corresponding score
     */
    @Requires({"first != null","second != null"})
    @Ensures("result != null")
    protected static Pair<Integer,Integer> findOverlappingSequence(final byte[] first,
                                                 final byte[] second) {
        final int MATCH_SCORE = 1;
        final int MISMATCH_SCORE = -1;
        // try every possible offset - O(N^2) algorithm

        // In case of following structure,
        //      111111111
        // 222222222
        // computed offset will be negative (=-5 in this case).
        // If however,
        //   111111111
        //      222222222
        // then offset will be positive (=3 in this case)
        int maxScore = 0, bestOffset =0;
        for (int offset = -second.length; offset < first.length; offset++) {
            int score = 0;
            // compute start index for each array
            int ind1 = (offset<0)?0:offset;
            int ind2 = (offset<0)?-offset:0;
            for (int k=0; k < Math.min(first.length, second.length) ; k++) {
                if (ind1 >= first.length)
                    break;
                if (ind2 >= second.length )
                    break;
                if (first[ind1] != 'N' && second[ind2] != 'N')  {
                    if (first[ind1] == second[ind2])
                        score += MATCH_SCORE;
                    else
                        score += MISMATCH_SCORE;
                }
                ind1++;
                ind2++;
            }
            if (score > maxScore) {
                maxScore = score;
                bestOffset = offset;
            }
        }
        return new Pair<Integer, Integer>(bestOffset,maxScore);
    }

}
