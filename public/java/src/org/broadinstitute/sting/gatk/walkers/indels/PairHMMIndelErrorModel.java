/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.PairHMM;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;


public class PairHMMIndelErrorModel {
    public static final int BASE_QUAL_THRESHOLD = 20;

    private boolean DEBUG = false;
    private boolean bandedLikelihoods = false;

    private static final int MAX_CACHED_QUAL = 127;

    private static final double baseMatchArray[];
    private static final double baseMismatchArray[];

    private final static double LOG_ONE_HALF;

    private static final int START_HRUN_GAP_IDX = 4;
    private static final int MAX_HRUN_GAP_IDX = 20;

    private static final byte MIN_GAP_OPEN_PENALTY = 30;
    private static final byte MIN_GAP_CONT_PENALTY = 10;
    private static final byte GAP_PENALTY_HRUN_STEP = 1; // each increase in hrun decreases gap penalty by this.

    private final byte[] GAP_OPEN_PROB_TABLE;
    private final byte[] GAP_CONT_PROB_TABLE;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////

    static {
        LOG_ONE_HALF= -Math.log10(2.0);

        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10, -k/10.);


            baseMatchArray[k] =  Math.log10(1-baseProb);
            baseMismatchArray[k] = Math.log10(baseProb);
        }
    }

    public PairHMMIndelErrorModel(byte indelGOP, byte indelGCP, boolean deb, boolean bandedLikelihoods) {
        this.DEBUG = deb;
        this.bandedLikelihoods = bandedLikelihoods;

        // fill gap penalty table, affine naive model:
        this.GAP_CONT_PROB_TABLE = new byte[MAX_HRUN_GAP_IDX];
        this.GAP_OPEN_PROB_TABLE = new byte[MAX_HRUN_GAP_IDX];


        for (int i = 0; i < START_HRUN_GAP_IDX; i++) {
            GAP_OPEN_PROB_TABLE[i] = indelGOP;
            GAP_CONT_PROB_TABLE[i] = indelGCP;
        }

        double step = GAP_PENALTY_HRUN_STEP/10.0;

        // initialize gop and gcp to their default values
        byte gop = indelGOP;
        byte gcp = indelGCP;

        // all of the following is computed in QUal-space
        for (int i=START_HRUN_GAP_IDX; i < MAX_HRUN_GAP_IDX; i++) {
            gop -= GAP_PENALTY_HRUN_STEP;
            if (gop < MIN_GAP_OPEN_PENALTY)
                gop = MIN_GAP_OPEN_PENALTY;

            gcp -= step;
            if(gcp < MIN_GAP_CONT_PENALTY)
                gcp = MIN_GAP_CONT_PENALTY;
            GAP_OPEN_PROB_TABLE[i] = gop;
            GAP_CONT_PROB_TABLE[i] = gcp;
        }

    }

    static private void getContextHomopolymerLength(final byte[] refBytes, final int[] hrunArray) {
        // compute forward hrun length, example:
        // AGGTGACCCCCCTGAGAG
        // 001000012345000000
        hrunArray[0] = 0;
        int[] hforward = new int[hrunArray.length];
        int[] hreverse = new int[hrunArray.length];

        for (int i = 1; i < refBytes.length; i++) {
            if (refBytes[i] == refBytes[i-1])
                hforward[i] = hforward[i-1]+1;
            else
                hforward[i] = 0;
        }

        // do similar thing for reverse length, example:
        // AGGTGACCCCCCTGAGAG
        // 021000543210000000
        // and then accumulate with forward values.
        // Total:
        // AGGTGACCCCCCTGAGAG
        // 022000555555000000
        for (int i=refBytes.length-1; i > 0; i--) {
            if (refBytes[i-1] == refBytes[i])
                hreverse[i-1] += hreverse[i]+1;
        }

        for (int i = 1; i < refBytes.length; i++)
            hrunArray[i] = hforward[i]+hreverse[i];
    }


    private void fillGapProbabilities(final int[] hrunProfile,
                                      final byte[] contextLogGapOpenProbabilities,
                                      final byte[] contextLogGapContinuationProbabilities) {
        // fill based on lookup table
        for (int i = 0; i < hrunProfile.length; i++) {
            if (hrunProfile[i] >= MAX_HRUN_GAP_IDX) {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
            }
            else {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[hrunProfile[i]];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[hrunProfile[i]];
            }
        }
    }


    public synchronized double[] computeDiploidReadHaplotypeLikelihoods(ReadBackedPileup pileup, LinkedHashMap<Allele, Haplotype> haplotypeMap, ReferenceContext ref, int eventLength, HashMap<PileupElement, LinkedHashMap<Allele, Double>> indelLikelihoodMap){
        final int numHaplotypes = haplotypeMap.size();

        final int readCounts[] = new int[pileup.getNumberOfElements()];
        final double[][] readLikelihoods = computeGeneralReadHaplotypeLikelihoods(pileup, haplotypeMap, ref, eventLength, indelLikelihoodMap, readCounts);
        return getDiploidHaplotypeLikelihoods(numHaplotypes, readCounts, readLikelihoods);
        
    }

    public synchronized double[][] computeGeneralReadHaplotypeLikelihoods(final ReadBackedPileup pileup, 
                                                                          final LinkedHashMap<Allele, Haplotype> haplotypeMap, 
                                                                          final ReferenceContext ref,
                                                                          final int eventLength, 
                                                                          final HashMap<PileupElement, LinkedHashMap<Allele, Double>> indelLikelihoodMap,
                                                                          final int[] readCounts) {
        final double readLikelihoods[][] = new double[pileup.getNumberOfElements()][haplotypeMap.size()];
        final PairHMM pairHMM = new PairHMM(bandedLikelihoods);

        int readIdx=0;
        for (PileupElement p: pileup) {
            // > 1 when the read is a consensus read representing multiple independent observations
            readCounts[readIdx] = p.getRepresentativeCount();

            // check if we've already computed likelihoods for this pileup element (i.e. for this read at this location)
            if (indelLikelihoodMap.containsKey(p)) {
                HashMap<Allele,Double> el = indelLikelihoodMap.get(p);
                int j=0;
                for (Allele a: haplotypeMap.keySet()) {
                    readLikelihoods[readIdx][j++] = el.get(a);
                }
            }
            else {
                final int refWindowStart = ref.getWindow().getStart();
                final int refWindowStop = ref.getWindow().getStop();

                if (DEBUG) {
                    System.out.format("Read Name:%s, aln start:%d aln stop:%d orig cigar:%s\n",p.getRead().getReadName(), p.getRead().getAlignmentStart(), p.getRead().getAlignmentEnd(), p.getRead().getCigarString());
                }

                GATKSAMRecord read = ReadClipper.hardClipAdaptorSequence(p.getRead());

                if (!read.isEmpty() && (read.getSoftEnd() > refWindowStop && read.getSoftStart() < refWindowStop))
                    read = ReadClipper.hardClipByReferenceCoordinatesRightTail(read, ref.getWindow().getStop());

                if (!read.isEmpty() && (read.getSoftStart() < refWindowStart && read.getSoftEnd() > refWindowStart))
                    read = ReadClipper.hardClipByReferenceCoordinatesLeftTail (read, ref.getWindow().getStart());

                if (read.isEmpty())
                    continue;

                // hard-clip low quality ends - this may introduce extra H elements in CIGAR string
                read = ReadClipper.hardClipLowQualEnds(read, (byte) BASE_QUAL_THRESHOLD );

                if (read.isEmpty())
                    continue;

                // get bases of candidate haplotypes that overlap with reads
                final int trailingBases = 3;
                final long readStart = read.getSoftStart();
                final long readEnd = read.getSoftEnd();

                // see if we want to use soft clipped bases. Aligners may soft clip all bases at insertions because they don't match,
                // but they're actually consistent with the insertion!
                // Rule: if a read starts in interval [eventStart-eventLength,eventStart+1] and we are at an insertion, we'll use all soft clipped bases at the beginning.
                // Conversely, if a read ends at [eventStart,eventStart+eventLength] we'll use all soft clipped bases in the end of the read.
                final long eventStartPos = ref.getLocus().getStart();

                // compute total number of clipped bases (soft or hard clipped) and only use them if necessary
                final boolean softClips = useSoftClippedBases(read, eventStartPos, eventLength);
                final int numStartSoftClippedBases = softClips ? read.getAlignmentStart()- read.getSoftStart() : 0;
                final int numEndSoftClippedBases = softClips ? read.getSoftEnd()- read.getAlignmentEnd() : 0 ;
                final byte [] unclippedReadBases = read.getReadBases();
                final byte [] unclippedReadQuals = read.getBaseQualities();
                final int extraOffset = Math.abs(eventLength);

                /**
                 * Compute genomic locations that candidate haplotypes will span.
                 * Read start and stop locations (variables readStart and readEnd) are the original unclipped positions from SAMRecord,
                 * adjusted by hard clips from Cigar string and by qual-based soft-clipping performed above.
                 * We will propose haplotypes that overlap the read with some padding.
                 * True read start = readStart + numStartSoftClippedBases - ReadUtils.getFirstInsertionOffset(read)
                 * Last term is because if a read starts with an insertion then these bases are not accounted for in readStart.
                 * trailingBases is a padding constant(=3) and we additionally add abs(eventLength) to both sides of read to be able to
                 * differentiate context between two haplotypes
                 */
                long startLocationInRefForHaplotypes = Math.max(readStart + numStartSoftClippedBases - trailingBases - ReadUtils.getFirstInsertionOffset(read)-extraOffset, 0);
                long stopLocationInRefForHaplotypes =  readEnd -numEndSoftClippedBases  + trailingBases + ReadUtils.getLastInsertionOffset(read)+extraOffset;

                if (DEBUG)
                    System.out.format("orig Start:%d orig stop: %d\n", startLocationInRefForHaplotypes, stopLocationInRefForHaplotypes);

                int readLength = read.getReadLength()-numStartSoftClippedBases-numEndSoftClippedBases;

                if (startLocationInRefForHaplotypes < ref.getWindow().getStart()) {
                    startLocationInRefForHaplotypes = ref.getWindow().getStart();                                       // read starts before haplotype: read will have to be cut numStartSoftClippedBases += ref.getWindow().getStart() - startLocationInRefForHaplotypes;
                }

                if (stopLocationInRefForHaplotypes > ref.getWindow().getStop()) {
                    stopLocationInRefForHaplotypes = ref.getWindow().getStop();                                         // check also if end of read will go beyond reference context
                }

                if (stopLocationInRefForHaplotypes <= startLocationInRefForHaplotypes + readLength) {
                    stopLocationInRefForHaplotypes = startLocationInRefForHaplotypes + readLength-1;                    // if there's an insertion in the read, the read stop position will be less than start + read legnth, but we want to compute likelihoods in the whole region that a read might overlap
                }

                // ok, we now figured out the total number of clipped bases on both ends.
                // Figure out where we want to place the haplotype to score read against

                if (DEBUG)
                    System.out.format("numStartSoftClippedBases: %d numEndSoftClippedBases: %d WinStart:%d WinStop:%d start: %d stop: %d readLength: %d\n",
                            numStartSoftClippedBases, numEndSoftClippedBases, ref.getWindow().getStart(), ref.getWindow().getStop(), startLocationInRefForHaplotypes, stopLocationInRefForHaplotypes, read.getReadLength());

                LinkedHashMap<Allele,Double> readEl = new LinkedHashMap<Allele,Double>();

                /**
                 * Check if we'll end up with an empty read once all clipping is done
                 */
                if (numStartSoftClippedBases + numEndSoftClippedBases >= unclippedReadBases.length) {
                    int j=0;
                    for (Allele a: haplotypeMap.keySet()) {
                        readEl.put(a,0.0);
                        readLikelihoods[readIdx][j++] = 0.0;
                    }
                }
                else {
                    final byte[] readBases = Arrays.copyOfRange(unclippedReadBases,numStartSoftClippedBases, unclippedReadBases.length-numEndSoftClippedBases);
                    final byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals,numStartSoftClippedBases, unclippedReadBases.length-numEndSoftClippedBases);
                    int j=0;

                    // initialize path metric and traceback memories for likelihood computation
                    double[][] matchMetricArray = null, XMetricArray = null, YMetricArray = null;
                    byte[] previousHaplotypeSeen = null;
                    final byte[] contextLogGapOpenProbabilities = new byte[readBases.length];
                    final byte[] contextLogGapContinuationProbabilities  = new byte[readBases.length];

                    // get homopolymer length profile for current haplotype
                    final int[] hrunProfile = new int[readBases.length];
                    getContextHomopolymerLength(readBases,hrunProfile);
                    fillGapProbabilities(hrunProfile, contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities);


                    for (Allele a: haplotypeMap.keySet()) {

                        Haplotype haplotype = haplotypeMap.get(a);

                        if (stopLocationInRefForHaplotypes > haplotype.getStopPosition())
                            stopLocationInRefForHaplotypes = haplotype.getStopPosition();

                        if (startLocationInRefForHaplotypes < haplotype.getStartPosition())
                            startLocationInRefForHaplotypes = haplotype.getStartPosition();

                        final long indStart = startLocationInRefForHaplotypes - haplotype.getStartPosition();
                        final long indStop =  stopLocationInRefForHaplotypes - haplotype.getStartPosition();

                        double readLikelihood;
                        if (DEBUG)
                            System.out.format("indStart: %d indStop: %d WinStart:%d WinStop:%d start: %d stop: %d readLength: %d C:%s\n",
                                    indStart, indStop, ref.getWindow().getStart(), ref.getWindow().getStop(), startLocationInRefForHaplotypes, stopLocationInRefForHaplotypes, read.getReadLength(), read.getCigar().toString());



                            final byte[] haplotypeBases = Arrays.copyOfRange(haplotype.getBases(),
                                    (int)indStart, (int)indStop);

                            final int X_METRIC_LENGTH = readBases.length+2;
                            final int Y_METRIC_LENGTH = haplotypeBases.length+2;

                            if (matchMetricArray == null) {
                                //no need to reallocate arrays for each new haplotype, as length won't change
                                matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                                XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                                YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];


                                PairHMM.initializeArrays(matchMetricArray, XMetricArray, YMetricArray, X_METRIC_LENGTH);
                            }

                            int startIndexInHaplotype = 0;
                            if (previousHaplotypeSeen != null)
                                startIndexInHaplotype = computeFirstDifferingPosition(haplotypeBases, previousHaplotypeSeen);
                            previousHaplotypeSeen = haplotypeBases.clone();

                            readLikelihood = pairHMM.computeReadLikelihoodGivenHaplotype(haplotypeBases, readBases, readQuals,
                                    contextLogGapOpenProbabilities, contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities,
                                    startIndexInHaplotype, matchMetricArray, XMetricArray, YMetricArray);


                            if (DEBUG) {
                                System.out.println("H:"+new String(haplotypeBases));
                                System.out.println("R:"+new String(readBases));
                                System.out.format("L:%4.2f\n",readLikelihood);
                                System.out.format("StPos:%d\n", startIndexInHaplotype);
                            }
                        readEl.put(a,readLikelihood);
                        readLikelihoods[readIdx][j++] = readLikelihood;
                    }
                }
                indelLikelihoodMap.put(p,readEl);
            }
            readIdx++;
        }

        if (DEBUG) {
            System.out.println("\nLikelihood summary");
            for (readIdx=0; readIdx < pileup.getNumberOfElements(); readIdx++) {
                System.out.format("Read Index: %d ",readIdx);
                for (int i=0; i < readLikelihoods[readIdx].length; i++)
                    System.out.format("L%d: %f ",i,readLikelihoods[readIdx][i]);
                System.out.println();
            }

        }

        return readLikelihoods;
    }

    private boolean useSoftClippedBases(GATKSAMRecord read, long eventStartPos, int eventLength) {
        return !((read.getAlignmentStart() >= eventStartPos-eventLength && read.getAlignmentStart() <= eventStartPos+1) || (read.getAlignmentEnd() >= eventStartPos && read.getAlignmentEnd() <= eventStartPos + eventLength));
    }

    private int computeFirstDifferingPosition(byte[] b1, byte[] b2) {
        if (b1.length != b2.length)
            return 0; // sanity check

        for (int i=0; i < b1.length; i++ ){
            if ( b1[i]!= b2[i] )
                return i;
        }
        return b1.length;
    }

    private static double[] getDiploidHaplotypeLikelihoods(final int numHaplotypes, final int readCounts[], final double readLikelihoods[][]) {
        final double[][] haplotypeLikehoodMatrix = new double[numHaplotypes][numHaplotypes];

        // todo: MAD 09/26/11 -- I'm almost certain this calculation can be simplified to just a single loop without the intermediate NxN matrix
        for (int i=0; i < numHaplotypes; i++) {
            for (int j=i; j < numHaplotypes; j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                // L(Hi, Hj) = sum_reads ( Pr(R|Hi)/2 + Pr(R|Hj)/2)
                //readLikelihoods[k][j] has log10(Pr(R_k) | H[j] )
                for (int readIdx = 0; readIdx < readLikelihoods.length; readIdx++) {
                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    if (Double.isInfinite(readLikelihoods[readIdx][i]) && Double.isInfinite(readLikelihoods[readIdx][j]))
                        continue;
                    final double li = readLikelihoods[readIdx][i];
                    final double lj = readLikelihoods[readIdx][j];
                    final int readCount = readCounts[readIdx];
                    haplotypeLikehoodMatrix[i][j] += readCount * (MathUtils.approximateLog10SumLog10(li, lj) + LOG_ONE_HALF);
                }
            }
        }

        final double[] genotypeLikelihoods = new double[numHaplotypes*(numHaplotypes+1)/2];
        int k=0;
        for (int j=0; j < numHaplotypes; j++) {
            for (int i=0; i <= j; i++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
            }
        }

        // renormalize so that max element is zero.
        return MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
    }
}
