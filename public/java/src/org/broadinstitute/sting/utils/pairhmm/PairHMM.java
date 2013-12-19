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

package org.broadinstitute.sting.utils.pairhmm;

import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.variantcontext.Allele;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
/**
 * Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
 *
 * User: rpoplin
 * Date: 10/16/12
 */
public abstract class PairHMM {
    protected final static Logger logger = Logger.getLogger(PairHMM.class);

    protected boolean constantsAreInitialized = false;

    protected byte[] previousHaplotypeBases;
    protected int hapStartIndex;

    public enum HMM_IMPLEMENTATION {
        /* Very slow implementation which uses very accurate log10 sum functions. Only meant to be used as a reference test implementation */
        EXACT,
        /* PairHMM as implemented for the UnifiedGenotyper. Uses log10 sum functions accurate to only 1E-4 */
        ORIGINAL,
        /* Optimized version of the PairHMM which caches per-read computations and operations in real space to avoid costly sums of log10'ed likelihoods */
        LOGLESS_CACHING,
        /* Logless caching PairHMM that stores computations in 1D arrays instead of matrices, and which proceeds diagonally over the (read x haplotype) intersection matrix */
        ARRAY_LOGLESS
    }

    protected int maxHaplotypeLength, maxReadLength;
    protected int paddedMaxReadLength, paddedMaxHaplotypeLength;
    protected int paddedReadLength, paddedHaplotypeLength;
    protected boolean initialized = false;

    // only used for debugging purposes
    protected boolean doNotUseTristateCorrection = false;
    protected void doNotUseTristateCorrection() { doNotUseTristateCorrection = true; }

    /**
     * Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
     *
     * Note: Do not worry about padding, just provide the true max length of the read and haplotype. The HMM will take care of the padding.
     *
     * @param haplotypeMaxLength the max length of haplotypes we want to use with this PairHMM
     * @param readMaxLength the max length of reads we want to use with this PairHMM
     */
    public void initialize( final int readMaxLength, final int haplotypeMaxLength ) {
        if ( readMaxLength <= 0 ) throw new IllegalArgumentException("READ_MAX_LENGTH must be > 0 but got " + readMaxLength);
        if ( haplotypeMaxLength <= 0 ) throw new IllegalArgumentException("HAPLOTYPE_MAX_LENGTH must be > 0 but got " + haplotypeMaxLength);

        maxHaplotypeLength = haplotypeMaxLength;
        maxReadLength = readMaxLength;

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        paddedMaxReadLength = readMaxLength + 1;
        paddedMaxHaplotypeLength = haplotypeMaxLength + 1;

        previousHaplotypeBases = null;

        constantsAreInitialized = false;
        initialized = true;
    }

    protected int findMaxReadLength(final List<GATKSAMRecord> reads) {
        int listMaxReadLength = 0;
        for(GATKSAMRecord read : reads){
            final int readLength = read.getReadLength();
            if( readLength > listMaxReadLength ) { listMaxReadLength = readLength; }
        }
        return listMaxReadLength;
    }

    protected int findMaxHaplotypeLength(final Map<Allele, Haplotype> haplotypeMap) {
        int listMaxHaplotypeLength = 0;
        for( final Allele a: haplotypeMap.keySet() ) {
            final Haplotype h = haplotypeMap.get(a);
            final int haplotypeLength = h.getBases().length;
            if( haplotypeLength > listMaxHaplotypeLength ) { listMaxHaplotypeLength = haplotypeLength; }
        }
        return listMaxHaplotypeLength;
    }

    /**
     *  Given a list of reads and haplotypes, for every read compute the total probability of said read arising from
     *  each haplotype given base substitution, insertion, and deletion probabilities.
     *
     * @param reads the list of reads
     * @param alleleHaplotypeMap the list of haplotypes
     * @param GCPArrayMap Each read is associated with an array containing the gap continuation penalties for use in the model. Length of each GCP-array must match that of its read.
     * @return a PerReadAlleleLikelihoodMap containing each read, haplotype-allele, and the log10 probability of
     *          said read coming from the said haplotype under the provided error model
     */
    public PerReadAlleleLikelihoodMap computeLikelihoods(final List<GATKSAMRecord> reads, final Map<Allele, Haplotype> alleleHaplotypeMap, final Map<GATKSAMRecord, byte[]> GCPArrayMap) {

        // (re)initialize the pairHMM only if necessary
        final int readMaxLength = findMaxReadLength(reads);
        final int haplotypeMaxLength = findMaxHaplotypeLength(alleleHaplotypeMap);
        if (!initialized || readMaxLength > maxReadLength || haplotypeMaxLength > maxHaplotypeLength) { initialize(readMaxLength, haplotypeMaxLength); }

        final PerReadAlleleLikelihoodMap likelihoodMap = new PerReadAlleleLikelihoodMap();
        for(GATKSAMRecord read : reads){
            final byte[] readBases = read.getReadBases();
            final byte[] readQuals = read.getBaseQualities();
            final byte[] readInsQuals = read.getBaseInsertionQualities();
            final byte[] readDelQuals = read.getBaseDeletionQualities();
            final byte[] overallGCP = GCPArrayMap.get(read);

            // peak at the next haplotype in the list (necessary to get nextHaplotypeBases, which is required for caching in the array implementation)
            byte[] currentHaplotypeBases = null;
            boolean isFirstHaplotype = true;
            Allele currentAllele = null;
            double log10l;
            for (final Allele allele : alleleHaplotypeMap.keySet()){
                final Haplotype haplotype = alleleHaplotypeMap.get(allele);
                final byte[] nextHaplotypeBases = haplotype.getBases();
                if (currentHaplotypeBases != null) {
                     log10l = computeReadLikelihoodGivenHaplotypeLog10(currentHaplotypeBases,
                            readBases, readQuals, readInsQuals, readDelQuals, overallGCP, isFirstHaplotype, nextHaplotypeBases);
                    likelihoodMap.add(read, currentAllele, log10l);
                }
                // update the current haplotype
                currentHaplotypeBases = nextHaplotypeBases;
                currentAllele = allele;
            }
            // process the final haplotype
            if (currentHaplotypeBases != null) {

                // there is no next haplotype, so pass null for nextHaplotypeBases.
                log10l = computeReadLikelihoodGivenHaplotypeLog10(currentHaplotypeBases,
                        readBases, readQuals, readInsQuals, readDelQuals, overallGCP, isFirstHaplotype, null);
                likelihoodMap.add(read, currentAllele, log10l);
            }
        }
        return likelihoodMap;
    }

    /**
     * Compute the total probability of read arising from haplotypeBases given base substitution, insertion, and deletion
     * probabilities.
     *
     * Note on using hapStartIndex.  This allows you to compute the exact true likelihood of a full haplotypes
     * given a read, assuming that the previous calculation read over a full haplotype, recaching the read values,
     * starting only at the place where the new haplotype bases and the previous haplotype bases different.  This
     * index is 0-based, and can be computed with findFirstPositionWhereHaplotypesDiffer given the two haplotypes.
     * Note that this assumes that the read and all associated quals values are the same.
     *
     * @param haplotypeBases the full sequence (in standard SAM encoding) of the haplotype, must be >= than read bases in length
     * @param readBases the bases (in standard encoding) of the read, must be <= haplotype bases in length
     * @param readQuals the phred-scaled per base substitution quality scores of read.  Must be the same length as readBases
     * @param insertionGOP the phred-scaled per base insertion quality scores of read.  Must be the same length as readBases
     * @param deletionGOP the phred-scaled per base deletion quality scores of read.  Must be the same length as readBases
     * @param overallGCP the phred-scaled gap continuation penalties scores of read.  Must be the same length as readBases
     * @param recacheReadValues if false, we don't recalculate any cached results, assuming that readBases and its associated
     *                          parameters are the same, and only the haplotype bases are changing underneath us
     * @return the log10 probability of read coming from the haplotype under the provided error model
     */
    protected final double computeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                  final byte[] readBases,
                                                                  final byte[] readQuals,
                                                                  final byte[] insertionGOP,
                                                                  final byte[] deletionGOP,
                                                                  final byte[] overallGCP,
                                                                  final boolean recacheReadValues,
                                                                  final byte[] nextHaploytpeBases) {

        if ( ! initialized ) throw new IllegalStateException("Must call initialize before calling computeReadLikelihoodGivenHaplotypeLog10");
        if ( haplotypeBases == null ) throw new IllegalArgumentException("haplotypeBases cannot be null");
        if ( haplotypeBases.length > maxHaplotypeLength ) throw new IllegalArgumentException("Haplotype bases is too long, got " + haplotypeBases.length + " but max is " + maxHaplotypeLength);
        if ( readBases == null ) throw new IllegalArgumentException("readBases cannot be null");
        if ( readBases.length > maxReadLength ) throw new IllegalArgumentException("readBases is too long, got " + readBases.length + " but max is " + maxReadLength);
        if ( readQuals.length != readBases.length ) throw new IllegalArgumentException("Read bases and read quals aren't the same size: " + readBases.length + " vs " + readQuals.length);
        if ( insertionGOP.length != readBases.length ) throw new IllegalArgumentException("Read bases and read insertion quals aren't the same size: " + readBases.length + " vs " + insertionGOP.length);
        if ( deletionGOP.length != readBases.length ) throw new IllegalArgumentException("Read bases and read deletion quals aren't the same size: " + readBases.length + " vs " + deletionGOP.length);
        if ( overallGCP.length != readBases.length ) throw new IllegalArgumentException("Read bases and overall GCP aren't the same size: " + readBases.length + " vs " + overallGCP.length);

        paddedReadLength = readBases.length + 1;
        paddedHaplotypeLength = haplotypeBases.length + 1;

        hapStartIndex =  (recacheReadValues) ? 0 : hapStartIndex;

        // Pre-compute the difference between the current haplotype and the next one to be run
        // Looking ahead is necessary for the ArrayLoglessPairHMM implementation
        final int nextHapStartIndex =  (nextHaploytpeBases == null || haplotypeBases.length != nextHaploytpeBases.length) ? 0 : findFirstPositionWhereHaplotypesDiffer(haplotypeBases, nextHaploytpeBases);

        double result = subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues, nextHapStartIndex);

        if ( ! MathUtils.goodLog10Probability(result) )
            throw new IllegalStateException("PairHMM Log Probability cannot be greater than 0: " + String.format("haplotype: %s, read: %s, result: %f", Arrays.toString(haplotypeBases), Arrays.toString(readBases), result));

        // Warning: Careful if using the PairHMM in parallel! (this update has to be taken care of).
        // Warning: This assumes no downstream modification of the haplotype bases (saves us from copying the array). It is okay for the haplotype caller and the Unified Genotyper.
        previousHaplotypeBases = haplotypeBases;

        // For the next iteration, the hapStartIndex for the next haploytpe becomes the index for the current haplotype
        // The array implementation has to look ahead to the next haplotype to store caching info. It cannot do this if nextHapStart is before hapStart
        hapStartIndex = (nextHapStartIndex < hapStartIndex) ? 0: nextHapStartIndex;

        return result;
    }

    /**
     * To be overloaded by subclasses to actually do calculation for #computeReadLikelihoodGivenHaplotypeLog10
     */
    @Requires({"readBases.length == readQuals.length", "readBases.length == insertionGOP.length", "readBases.length == deletionGOP.length",
            "readBases.length == overallGCP.length", "matchMatrix!=null", "insertionMatrix!=null", "deletionMatrix!=null"})
    protected abstract double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                           final byte[] readBases,
                                                                           final byte[] readQuals,
                                                                           final byte[] insertionGOP,
                                                                           final byte[] deletionGOP,
                                                                           final byte[] overallGCP,
                                                                           final int hapStartIndex,
                                                                           final boolean recacheReadValues,
                                                                           final int nextHapStartIndex);

    /**
     * Compute the first position at which two haplotypes differ
     *
     * If the haplotypes are exact copies of each other, returns the min length of the two haplotypes.
     *
     * @param haplotype1 the first haplotype1
     * @param haplotype2 the second haplotype1
     * @return the index of the first position in haplotype1 and haplotype2 where the byte isn't the same
     */
    public static int findFirstPositionWhereHaplotypesDiffer(final byte[] haplotype1, final byte[] haplotype2) {
        if ( haplotype1 == null || haplotype1.length == 0 ) throw new IllegalArgumentException("Haplotype1 is bad " + Arrays.toString(haplotype1));
        if ( haplotype2 == null || haplotype2.length == 0 ) throw new IllegalArgumentException("Haplotype2 is bad " + Arrays.toString(haplotype2));

        for( int iii = 0; iii < haplotype1.length && iii < haplotype2.length; iii++ ) {
            if( haplotype1[iii] != haplotype2[iii] ) {
                return iii;
            }
        }

        return Math.min(haplotype1.length, haplotype2.length);
    }
}
