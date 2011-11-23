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

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
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
    private final static double END_GAP_COST;

    private static final int START_HRUN_GAP_IDX = 4;
    private static final int MAX_HRUN_GAP_IDX = 20;

    private static final double MIN_GAP_OPEN_PENALTY = 30.0;
    private static final double MIN_GAP_CONT_PENALTY = 10.0;
    private static final double GAP_PENALTY_HRUN_STEP = 1.0; // each increase in hrun decreases gap penalty by this.

    private final double[] GAP_OPEN_PROB_TABLE;
    private final double[] GAP_CONT_PROB_TABLE;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////

    static {
        LOG_ONE_HALF= -Math.log10(2.0);
        END_GAP_COST = LOG_ONE_HALF;

        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10, -k/10.);


            baseMatchArray[k] =  Math.log10(1-baseProb);
            baseMismatchArray[k] = Math.log10(baseProb);
        }
    }

    public PairHMMIndelErrorModel(double indelGOP, double indelGCP, boolean deb, boolean bandedLikelihoods) {
        this.DEBUG = deb;
        this.bandedLikelihoods = bandedLikelihoods;

        // fill gap penalty table, affine naive model:
        this.GAP_CONT_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];
        this.GAP_OPEN_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];

        double gop = -indelGOP/10.0;
        double gcp = -indelGCP/10.0;

        for (int i = 0; i < START_HRUN_GAP_IDX; i++) {
            GAP_OPEN_PROB_TABLE[i] = gop;
            GAP_CONT_PROB_TABLE[i] = gcp;
        }

        double step = GAP_PENALTY_HRUN_STEP/10.0;

        double maxGOP = -MIN_GAP_OPEN_PENALTY/10.0;  // phred to log prob
        double maxGCP = -MIN_GAP_CONT_PENALTY/10.0;  // phred to log prob

        for (int i=START_HRUN_GAP_IDX; i < MAX_HRUN_GAP_IDX; i++) {
            gop += step;
            if (gop > maxGOP)
                gop = maxGOP;

            gcp += step;
            if(gcp > maxGCP)
                gcp = maxGCP;
            GAP_OPEN_PROB_TABLE[i] = gop;
            GAP_CONT_PROB_TABLE[i] = gcp;
        }

    }

    static private void getContextHomopolymerLength(final byte[] refBytes, int[] hrunArray) {
        // compute forward hrun length, example:
        // AGGTGACCCCCCTGAGAG
        // 001000012345000000
        int runCount = 0;
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


    private void updateCell(final int indI, final int indJ, final int X_METRIC_LENGTH, final int Y_METRIC_LENGTH, byte[] readBases, byte[] readQuals, byte[] haplotypeBases,
                            double[] currentGOP, double[] currentGCP,  double[][] matchMetricArray,  double[][] XMetricArray,  double[][] YMetricArray) {
        if (indI > 0 && indJ > 0) {
            final int im1 = indI -1;
            final int jm1 = indJ - 1;
            // update current point
            final byte x = readBases[im1];
            final byte y = haplotypeBases[jm1];
            final byte qual = readQuals[im1] < 1 ? 1 : (readQuals[im1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[im1]);

            final double pBaseRead =  (x == y)? baseMatchArray[(int)qual]:baseMismatchArray[(int)qual];

            matchMetricArray[indI][indJ] = MathUtils.softMax(matchMetricArray[im1][jm1] + pBaseRead, XMetricArray[im1][jm1] + pBaseRead,
                    YMetricArray[im1][jm1] + pBaseRead);

            final double c1 = indJ == Y_METRIC_LENGTH-1 ? END_GAP_COST : currentGOP[jm1];
            final double d1 = indJ == Y_METRIC_LENGTH-1 ? END_GAP_COST : currentGCP[jm1];

            XMetricArray[indI][indJ] = MathUtils.softMax(matchMetricArray[im1][indJ] + c1, XMetricArray[im1][indJ] + d1);

            // update Y array
            final double c2 = indI == X_METRIC_LENGTH-1 ? END_GAP_COST : currentGOP[jm1];
            final double d2 = indI == X_METRIC_LENGTH-1 ? END_GAP_COST : currentGCP[jm1];
            YMetricArray[indI][indJ] = MathUtils.softMax(matchMetricArray[indI][jm1] + c2, YMetricArray[indI][jm1] + d2);
        }
    }

    private double computeReadLikelihoodGivenHaplotypeAffineGaps(byte[] haplotypeBases, byte[] readBases, byte[] readQuals,
                                                                 double[] currentGOP, double[] currentGCP, int indToStart,
                                                                 double[][] matchMetricArray, double[][] XMetricArray, double[][] YMetricArray) {

        final int X_METRIC_LENGTH = readBases.length+1;
        final int Y_METRIC_LENGTH = haplotypeBases.length+1;

        if (indToStart == 0) {
            // default initialization for all arrays

            for (int i=0; i < X_METRIC_LENGTH; i++) {
                Arrays.fill(matchMetricArray[i],Double.NEGATIVE_INFINITY);
                Arrays.fill(YMetricArray[i],Double.NEGATIVE_INFINITY);
                Arrays.fill(XMetricArray[i],Double.NEGATIVE_INFINITY);
            }

            for (int i=1; i < X_METRIC_LENGTH; i++) {
                //initialize first column
                XMetricArray[i][0]      = END_GAP_COST*(i);
            }

            for (int j=1; j < Y_METRIC_LENGTH; j++) {
                // initialize first row
                YMetricArray[0][j]      = END_GAP_COST*(j);
            }
            matchMetricArray[0][0]= END_GAP_COST;//Double.NEGATIVE_INFINITY;
            XMetricArray[0][0]=  YMetricArray[0][0] = 0;
        }


        if (bandedLikelihoods) {
            final double DIAG_TOL = 20; // means that max - min element in diags have to be > this number for banding to take effect.

            final int numDiags = X_METRIC_LENGTH +  Y_METRIC_LENGTH -1;
            final int elemsInDiag = Math.min(X_METRIC_LENGTH, Y_METRIC_LENGTH);

            int idxWithMaxElement = 0;

            for (int  diag=indToStart; diag <  numDiags; diag++) {
                // compute default I and J start positions at edge of diagonals
                int indI = 0;
                int indJ = diag;
                if (diag >= Y_METRIC_LENGTH ) {
                    indI = diag-(Y_METRIC_LENGTH-1);
                    indJ = Y_METRIC_LENGTH-1;
                }

                // first pass: from max element to edge
                int idxLow =  idxWithMaxElement;

                // reset diag max value before starting
                double maxElementInDiag = Double.NEGATIVE_INFINITY;
                // set indI, indJ to correct values
                indI += idxLow;
                indJ -= idxLow;
                if (indI >= X_METRIC_LENGTH || indJ < 0) {
                    idxLow--;
                    indI--;
                    indJ++;
                }


                for (int el = idxLow; el < elemsInDiag; el++) {
                    updateCell(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, readBases, readQuals, haplotypeBases,
                            currentGOP, currentGCP,  matchMetricArray,  XMetricArray, YMetricArray);
                    // update max in diagonal
                    final double bestMetric = MathUtils.max(matchMetricArray[indI][indJ], XMetricArray[indI][indJ], YMetricArray[indI][indJ]);

                    // check if we've fallen off diagonal value by threshold
                    if (bestMetric > maxElementInDiag) {
                        maxElementInDiag = bestMetric;
                        idxWithMaxElement = el;
                    }
                    else if (bestMetric < maxElementInDiag - DIAG_TOL && idxWithMaxElement > 0)
                        break; // done w/current diagonal

                    indI++;
                    if (indI >=X_METRIC_LENGTH )
                        break;
                    indJ--;
                    if (indJ <= 0)
                        break;
                }
                if (idxLow > 0) {
                    // now do second part in opposite direction
                    indI = 0;
                    indJ = diag;
                    if (diag >= Y_METRIC_LENGTH ) {
                        indI = diag-(Y_METRIC_LENGTH-1);
                        indJ = Y_METRIC_LENGTH-1;
                    }

                    indI += idxLow-1;
                    indJ -= idxLow-1;
                    for (int el = idxLow-1; el >= 0; el--) {

                        updateCell(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, readBases, readQuals, haplotypeBases,
                                currentGOP, currentGCP,  matchMetricArray,  XMetricArray, YMetricArray);
                        // update max in diagonal
                        final double bestMetric = MathUtils.max(matchMetricArray[indI][indJ], XMetricArray[indI][indJ], YMetricArray[indI][indJ]);

                        // check if we've fallen off diagonal value by threshold
                        if (bestMetric > maxElementInDiag) {
                            maxElementInDiag = bestMetric;
                            idxWithMaxElement = el;
                        }
                        else if (bestMetric < maxElementInDiag - DIAG_TOL)
                            break; // done w/current diagonal

                        indJ++;
                        if (indJ >= Y_METRIC_LENGTH )
                            break;
                        indI--;
                        if (indI <= 0)
                            break;
                    }
                }
                // if (DEBUG)
                //     System.out.format("Max:%4.1f el:%d\n",maxElementInDiag,  idxWithMaxElement);
            }
        }
        else {
            // simplified rectangular version of update loop
            for (int indI=1; indI < X_METRIC_LENGTH; indI++) {
                for (int indJ=indToStart+1; indJ < Y_METRIC_LENGTH; indJ++) {
                    updateCell(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, readBases, readQuals, haplotypeBases,
                            currentGOP, currentGCP,  matchMetricArray,  XMetricArray, YMetricArray);

                }
            }
        }



        final int bestI = X_METRIC_LENGTH - 1, bestJ = Y_METRIC_LENGTH - 1;
        final double bestMetric = MathUtils.softMax(matchMetricArray[bestI][bestJ],
                XMetricArray[bestI][bestJ],
                YMetricArray[bestI][bestJ]);

        /*
        if (DEBUG) {
            PrintStream outx, outy, outm, outs;
            double[][] sumMetrics = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
            try {
                outx = new PrintStream("datax.txt");
                outy = new PrintStream("datay.txt");
                outm = new PrintStream("datam.txt");
                outs = new PrintStream("datas.txt");
                double metrics[] = new double[3];
                for (int indI=0; indI < X_METRIC_LENGTH; indI++) {
                    for (int indJ=0; indJ < Y_METRIC_LENGTH; indJ++) {
                        metrics[0] = matchMetricArray[indI][indJ];
                        metrics[1] = XMetricArray[indI][indJ];
                        metrics[2] = YMetricArray[indI][indJ];
                        //sumMetrics[indI][indJ] = MathUtils.softMax(metrics);
                        outx.format("%4.1f ", metrics[1]);
                        outy.format("%4.1f ", metrics[2]);
                        outm.format("%4.1f ", metrics[0]);
                        outs.format("%4.1f ", MathUtils.softMax(metrics));
                    }
                    outx.println();  outm.println();outy.println(); outs.println();
                }
                outm.close(); outx.close(); outy.close();
            } catch (java.io.IOException e) { throw new UserException("bla");}
        }
        */

        return bestMetric;

    }

    private void fillGapProbabilities(int[] hrunProfile,
                                      double[] contextLogGapOpenProbabilities, double[] contextLogGapContinuationProbabilities) {
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
    public synchronized double[] computeReadHaplotypeLikelihoods(ReadBackedPileup pileup, LinkedHashMap<Allele,Haplotype> haplotypeMap,
                                                                 ReferenceContext ref, int eventLength,
                                                                 HashMap<PileupElement, LinkedHashMap<Allele,Double>> indelLikelihoodMap){

        int numHaplotypes = haplotypeMap.size();
        final double readLikelihoods[][] = new double[pileup.getNumberOfElements()][numHaplotypes];
        final int readCounts[] = new int[pileup.getNumberOfElements()];
        int readIdx=0;

        LinkedHashMap<Allele,double[]> gapOpenProbabilityMap = new LinkedHashMap<Allele,double[]>();
        LinkedHashMap<Allele,double[]> gapContProbabilityMap = new LinkedHashMap<Allele,double[]>();

        // will context dependent probabilities based on homopolymer run. Probabilities are filled based on total complete haplotypes.
        // todo -- refactor into separate function
        for (Allele a: haplotypeMap.keySet()) {
            Haplotype haplotype = haplotypeMap.get(a);
            byte[] haplotypeBases = haplotype.getBasesAsBytes();
            double[] contextLogGapOpenProbabilities = new double[haplotypeBases.length];
            double[] contextLogGapContinuationProbabilities = new double[haplotypeBases.length];

            // get homopolymer length profile for current haplotype
            int[] hrunProfile = new int[haplotypeBases.length];
            getContextHomopolymerLength(haplotypeBases,hrunProfile);
            fillGapProbabilities(hrunProfile, contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities);

            gapOpenProbabilityMap.put(a,contextLogGapOpenProbabilities);
            gapContProbabilityMap.put(a,contextLogGapContinuationProbabilities);

        }

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
                //System.out.format("%d %s\n",p.getRead().getAlignmentStart(), p.getRead().getClass().getName());
                SAMRecord read = ReadUtils.hardClipAdaptorSequence(p.getRead());
                if (read == null)
                    continue;

                if(ReadUtils.is454Read(read)) {
                    continue;
                }

                // get bases of candidate haplotypes that overlap with reads
                final int trailingBases = 3;

                long readStart = read.getUnclippedStart();
                long readEnd = read.getUnclippedEnd();

                int numStartSoftClippedBases, numEndSoftClippedBases;

                // see if we want to use soft clipped bases. Aligners may soft clip all bases at insertions because they don't match,
                // but they're actually consistent with the insertion!
                // Rule: if a read starts in interval [eventStart-eventLength,eventStart+1] and we are at an insertion, we'll use all soft clipped bases at the beginning.
                // Conversely, if a read ends at [eventStart,eventStart+eventLength] we'll use all soft clipped bases in the end of the read.
                long eventStartPos = ref.getLocus().getStart();

                // compute total number of clipped bases (soft or hard clipped)
                numStartSoftClippedBases = read.getAlignmentStart()- read.getUnclippedStart();
                numEndSoftClippedBases = read.getUnclippedEnd()- read.getAlignmentEnd();

                // check for hard clips (never consider these bases):
                Cigar c = read.getCigar();
                CigarElement first = c.getCigarElement(0);
                CigarElement last = c.getCigarElement(c.numCigarElements()-1);
                int numStartHardClippedBases = 0, numEndHardClippedBases = 0;

                if (first.getOperator() == CigarOperator.H) {
                    numStartHardClippedBases = first.getLength();
                }

                if (last.getOperator() == CigarOperator.H) {
                    numEndHardClippedBases = last.getLength();
                }

                // correct for hard clips
                numStartSoftClippedBases -= numStartHardClippedBases;
                numEndSoftClippedBases -= numEndHardClippedBases;
                readStart += numStartHardClippedBases;
                readEnd -= numEndHardClippedBases;

                // remove soft clips if necessary
                if ((read.getAlignmentStart()>=eventStartPos-eventLength && read.getAlignmentStart() <= eventStartPos+1) ||
                        (read.getAlignmentEnd() >= eventStartPos && read.getAlignmentEnd() <= eventStartPos + eventLength)) {
                    numStartSoftClippedBases = 0;
                    numEndSoftClippedBases = 0;
                }



                byte[] unclippedReadBases, unclippedReadQuals;

                int numStartClippedBases = numStartSoftClippedBases;
                int numEndClippedBases = numEndSoftClippedBases;
                unclippedReadBases = read.getReadBases();
                unclippedReadQuals = read.getBaseQualities();

                // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
                // and may leave a string of Q2 bases still hanging off the reads.
                for (int i=numStartSoftClippedBases; i < unclippedReadBases.length; i++) {
                    if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                        numStartClippedBases++;
                    else
                        break;

                }
                for (int i=unclippedReadBases.length-numEndSoftClippedBases-1; i >= 0; i-- ){
                    if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                        numEndClippedBases++;
                    else
                        break;
                }

                int extraOffset = Math.abs(eventLength);

                long start = Math.max(readStart + numStartClippedBases - trailingBases - ReadUtils.getFirstInsertionOffset(read)-extraOffset, 0);
                long stop =  readEnd -numEndClippedBases  + trailingBases + ReadUtils.getLastInsertionOffset(read)+extraOffset;

                // Variables start and stop are coordinates (inclusive) where we want to get the haplotype from.
                int readLength = read.getReadLength()-numStartSoftClippedBases-numEndSoftClippedBases;
                // check if start of read will be before start of reference context
                if (start < ref.getWindow().getStart())// read starts before haplotype: read will have to be cut
                    start = ref.getWindow().getStart();

                // check also if end of read will go beyond reference context
                if (stop > ref.getWindow().getStop())
                    stop = ref.getWindow().getStop();

                // if there's an insertion in the read, the read stop position will be less than start + read legnth,
                // but we want to compute likelihoods in the whole region that a read might overlap
                if (stop <= start + readLength) {
                    stop = start + readLength-1;
                }

                // ok, we now figured out total number of clipped bases on both ends.
                // Figure out where we want to place the haplotype to score read against
                /*
               if (DEBUG)
                   System.out.format("numStartClippedBases: %d numEndClippedBases: %d WinStart:%d WinStop:%d start: %d stop: %d readLength: %d\n",
                           numStartClippedBases, numEndClippedBases, ref.getWindow().getStart(), ref.getWindow().getStop(), start, stop, read.getReadLength());
                */


                LinkedHashMap<Allele,Double> readEl = new LinkedHashMap<Allele,Double>();

                if (numStartClippedBases + numEndClippedBases >= unclippedReadBases.length) {
                    int j=0;
                    for (Allele a: haplotypeMap.keySet()) {
                        readEl.put(a,0.0);
                        readLikelihoods[readIdx][j++] = 0.0;
                    }

                }
                else {
                    final byte[] readBases = Arrays.copyOfRange(unclippedReadBases,numStartClippedBases,
                            unclippedReadBases.length-numEndClippedBases);

                    final byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals,numStartClippedBases,
                            unclippedReadBases.length-numEndClippedBases);

                    int j=0;

                    // initialize path metric and traceback memories for likelihood computation
                    double[][] matchMetricArray = null, XMetricArray = null, YMetricArray = null;
                    byte[] previousHaplotypeSeen = null;
                    double[] previousGOP = null;
                    double[] previousGCP = null;
                    int startIdx;
                    for (Allele a: haplotypeMap.keySet()) {


                        Haplotype haplotype = haplotypeMap.get(a);
                        if (stop > haplotype.getStopPosition())
                            stop = haplotype.getStopPosition();

                        if (start < haplotype.getStartPosition())
                            start = haplotype.getStartPosition();

                        // cut haplotype bases
                        long indStart = start - haplotype.getStartPosition();
                        long indStop =  stop - haplotype.getStartPosition();
                        double readLikelihood;
                        if (DEBUG)
                            System.out.format("indStart: %d indStop: %d WinStart:%d WinStop:%d start: %d stop: %d readLength: %d C:%s\n",
                                    indStart, indStop, ref.getWindow().getStart(), ref.getWindow().getStop(), start, stop, read.getReadLength(), read.getCigar().toString());

                        if (indStart < 0 || indStop >= haplotype.getBasesAsBytes().length || indStart > indStop) {
                            // read spanned more than allowed reference context: we currently can't deal with this
                            readLikelihood =0;
                        } else
                        {
                            final byte[] haplotypeBases = Arrays.copyOfRange(haplotype.getBasesAsBytes(),
                                    (int)indStart, (int)indStop);

                            if (matchMetricArray == null) {
                                final int X_METRIC_LENGTH = readBases.length+1;
                                final int Y_METRIC_LENGTH = haplotypeBases.length+1;

                                matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                                XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                                YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                            }
                            final double[] currentContextGOP = Arrays.copyOfRange(gapOpenProbabilityMap.get(a), (int)indStart, (int)indStop);
                            final double[] currentContextGCP = Arrays.copyOfRange(gapContProbabilityMap.get(a), (int)indStart, (int)indStop);
                            if (previousHaplotypeSeen == null)
                                startIdx = 0;
                            else {
                                final int s1 = computeFirstDifferingPosition(haplotypeBases, previousHaplotypeSeen);
                                final int s2 = computeFirstDifferingPosition(currentContextGOP, previousGOP);
                                final int s3 = computeFirstDifferingPosition(currentContextGCP, previousGCP);
                                startIdx = Math.min(Math.min(s1, s2), s3);
                            }
                            previousHaplotypeSeen = haplotypeBases.clone();
                            previousGOP = currentContextGOP.clone();
                            previousGCP = currentContextGCP.clone();


                            readLikelihood = computeReadLikelihoodGivenHaplotypeAffineGaps(haplotypeBases, readBases, readQuals,
                                    currentContextGOP, currentContextGCP, startIdx, matchMetricArray, XMetricArray, YMetricArray);

                            if (DEBUG) {
                                System.out.println("H:"+new String(haplotypeBases));
                                System.out.println("R:"+new String(readBases));
                                System.out.format("L:%4.2f\n",readLikelihood);
                                System.out.format("StPos:%d\n", startIdx);
                            }
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

        return getHaplotypeLikelihoods(numHaplotypes, readCounts, readLikelihoods);
    }

    private int computeFirstDifferingPosition(byte[] b1, byte[] b2) {
        if (b1.length != b2.length)
            return 0; // sanity check

        for (int i=0; i < b1.length; i++ ){
            if ( b1[i]!= b2[i])
                return i;
        }
        return b1.length;
    }

    private int computeFirstDifferingPosition(double[] b1, double[] b2) {
        if (b1.length != b2.length)
            return 0; // sanity check

        for (int i=0; i < b1.length; i++ ){
            if ( b1[i]!= b2[i])
                return i;
        }
        return b1.length;
    }

    private final static double[] getHaplotypeLikelihoods(final int numHaplotypes, final int readCounts[], final double readLikelihoods[][]) {
        final double[][] haplotypeLikehoodMatrix = new double[numHaplotypes][numHaplotypes];

        // todo: MAD 09/26/11 -- I'm almost certain this calculation can be simplied to just a single loop without the intermediate NxN matrix
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
                    haplotypeLikehoodMatrix[i][j] += readCount * (MathUtils.softMax(li, lj) + LOG_ONE_HALF);
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

        // renormalize   so that max element is zero.
        return MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
    }
}
