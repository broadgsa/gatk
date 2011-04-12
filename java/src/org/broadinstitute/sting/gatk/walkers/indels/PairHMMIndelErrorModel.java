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
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.Haplotype;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.Arrays;
import java.util.List;


public class PairHMMIndelErrorModel {


    private final int BASE_QUAL_THRESHOLD = 10;


    private static final int MATCH_OFFSET = 0;
    private static final int X_OFFSET = 1;
    private static final int Y_OFFSET = 2;

    private static final int DIAG = 0;
    private static final int UP = 1;
    private static final int LEFT = 2;

    private static final int DIAG_GOTO_M = 0;
    private static final int DIAG_GOTO_X = 1;
    private static final int DIAG_GOTO_Y = 2;

    private static final int UP_GOTO_M = 4;
    private static final int UP_GOTO_X = 5;
    private static final int UP_GOTO_Y = 6;

    private static final int LEFT_GOTO_M = 8;
    private static final int LEFT_GOTO_X = 9;
    private static final int LEFT_GOTO_Y = 10;

    private static final int[] ACTIONS_M = {DIAG_GOTO_M, DIAG_GOTO_X, DIAG_GOTO_Y};
    private static final int[] ACTIONS_X = {UP_GOTO_M, UP_GOTO_X, UP_GOTO_Y};
    private static final int[] ACTIONS_Y = {LEFT_GOTO_M, LEFT_GOTO_X, LEFT_GOTO_Y};


    private final double logGapOpenProbability;
    private final double logGapContinuationProbability;

    private boolean DEBUG = false;

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


    private boolean doViterbi = false;

    private final boolean useAffineGapModel = true;
    private boolean doContextDependentPenalties = false;

    private final double[] GAP_OPEN_PROB_TABLE;
    private final double[] GAP_CONT_PROB_TABLE;


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

    public  PairHMMIndelErrorModel(double indelGOP, double indelGCP, boolean deb, boolean doCDP, boolean dovit) {
        this(indelGOP, indelGCP, deb, doCDP);
        this.doViterbi = dovit;
    }

    public  PairHMMIndelErrorModel(double indelGOP, double indelGCP, boolean deb, boolean doCDP) {


        this.logGapOpenProbability = -indelGOP/10.0; // QUAL to log prob
        this.logGapContinuationProbability = -indelGCP/10.0; // QUAL to log prob
        this.doContextDependentPenalties = doCDP;
        this.DEBUG = deb;


        // fill gap penalty table, affine naive model:
        this.GAP_CONT_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];
        this.GAP_OPEN_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];

        for (int i = 0; i < START_HRUN_GAP_IDX; i++) {
            GAP_OPEN_PROB_TABLE[i] = logGapOpenProbability;
            GAP_CONT_PROB_TABLE[i] = logGapContinuationProbability;
        }

        double gop = logGapOpenProbability;
        double gcp = logGapContinuationProbability;
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

    private double computeReadLikelihoodGivenHaplotype(byte[] haplotypeBases, byte[] readBases, byte[] readQuals) {
        final int X_METRIC_LENGTH = readBases.length+1;
        final int Y_METRIC_LENGTH = haplotypeBases.length+1;

        // initialize path metric and traceback memories for likelihood computation
        double[][] pathMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestMetricArray = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        pathMetricArray[0][0]= 0;//Double.NEGATIVE_INFINITY;

        for (int i=1; i < X_METRIC_LENGTH; i++) {
            pathMetricArray[i][0] = 0;
            bestMetricArray[i][0] = UP;
        }

        for (int j=1; j < Y_METRIC_LENGTH; j++) {
            pathMetricArray[0][j] = 0;//logGapOpenProbability + (j-1) * logGapContinuationProbability;
            bestMetricArray[0][j] = LEFT;
        }

        for (int indI=1; indI < X_METRIC_LENGTH; indI++) {
            for (int indJ=1; indJ < Y_METRIC_LENGTH; indJ++) {

                byte x = readBases[indI-1];
                byte y = haplotypeBases[indJ-1];
                byte qual = readQuals[indI-1];

                double bestMetric = 0.0;
                int bestMetricIdx = 0;

                // compute metric for match/mismatch
                // workaround for reads whose bases quality = 0,
                if (qual < 1)
                    qual = 1;

                if (qual > MAX_CACHED_QUAL)
                    qual = MAX_CACHED_QUAL;

                double pBaseRead =  (x == y)? baseMatchArray[(int)qual]:baseMismatchArray[(int)qual];
                double[] metrics = new double[3];

                metrics[DIAG] = pathMetricArray[indI-1][indJ-1] + pBaseRead;
                metrics[UP] = pathMetricArray[indI-1][indJ] + logGapOpenProbability;//(end?0.0:logGapOpenProbability);
                metrics[LEFT] = pathMetricArray[indI][indJ-1] + logGapOpenProbability;//(end?0.0:logGapOpenProbability);

                if (doViterbi) {
                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(metrics);

                pathMetricArray[indI][indJ] = bestMetric;
                bestMetricArray[indI][indJ] = bestMetricIdx;

            }
        }


        double bestMetric=0.0;
        int bestMetricIdx=0,bestI=X_METRIC_LENGTH - 1, bestJ=Y_METRIC_LENGTH - 1;

        for (int i=0; i < X_METRIC_LENGTH; i ++ ) {
            int j= Y_METRIC_LENGTH-1;

            if (pathMetricArray[i][j] > bestMetric) {
                bestMetric = pathMetricArray[i][j];
                bestI = i;
                bestJ = j;
            }
        }
        for (int j=0; j < Y_METRIC_LENGTH; j++ ) {
            int i= X_METRIC_LENGTH-1;
            if (pathMetricArray[i][j] >= bestMetric) {
                bestMetric = pathMetricArray[i][j];
                bestI = i;
                bestJ = j;
            }
        }

        if (DEBUG && doViterbi) {

            String haplotypeString = new String (haplotypeBases);
            String readString = new String(readBases);


            int i = bestI;
            int j = bestJ;


            System.out.println("Simple NW");

            while (i >0 || j >0) {
                bestMetricIdx = bestMetricArray[i][j];
                System.out.print(bestMetricIdx);
                if (bestMetricIdx == UP) {
                    // insert gap in Y
                    haplotypeString = haplotypeString.substring(0,j)+"-"+haplotypeString.substring(j);
                    i--;
                } else if (bestMetricIdx == LEFT) {
                    readString = readString.substring(0,i)+"-"+readString.substring(i);
                    j--;
                }
                else {
                    i--; j--;
                }
            }




            System.out.println("\nAlignment: ");
            System.out.println("R:"+readString);
            System.out.println("H:"+haplotypeString);
            System.out.println();


        }
        if (DEBUG)
            System.out.format("Likelihood: %5.4f\n", bestMetric);

        return bestMetric;


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
    private Pair<Double,Double> getGapPenalties(final int indI, final int indJ, final int X_METRIC_LENGTH,
                                                final int Y_METRIC_LENGTH, final int tableToUpdate, double[] currentContextGOP, double[] currentContextGCP) {

        double c=0.0,d=0.0;

        if (doContextDependentPenalties) {
            switch (tableToUpdate) {
                case MATCH_OFFSET:

                    break;
                case X_OFFSET:
                    c = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: currentContextGOP[indJ-1]);
                    d = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: currentContextGCP[indJ-1]);

                    break;

                case Y_OFFSET:
                    c = (indI==X_METRIC_LENGTH-1? END_GAP_COST: currentContextGOP[indJ-1]);
                    d = (indI==X_METRIC_LENGTH-1? END_GAP_COST: currentContextGCP[indJ-1]);

                    break;
                default:
                    throw new StingException("BUG!! Invalid table offset");
            }
        }
        else {
            switch (tableToUpdate) {
                case MATCH_OFFSET:

                    break;
                case X_OFFSET:
                    c = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: logGapOpenProbability);
                    d = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: logGapContinuationProbability);

                    break;

                case Y_OFFSET:
                    c = (indI==X_METRIC_LENGTH-1? END_GAP_COST: logGapOpenProbability);
                    d = (indI==X_METRIC_LENGTH-1? END_GAP_COST: logGapContinuationProbability);

                    break;
                default:
                    throw new StingException("BUG!! Invalid table offset");
            }
        }
        return new Pair<Double,Double>(Double.valueOf(c),Double.valueOf(d));
    }

    private double computeReadLikelihoodGivenHaplotypeAffineGaps(byte[] haplotypeBases, byte[] readBases, byte[] readQuals,
                                                                 double[] currentGOP, double[] currentGCP) {

        final int X_METRIC_LENGTH = readBases.length+1;
        final int Y_METRIC_LENGTH = haplotypeBases.length+1;

        // initialize path metric and traceback memories for likelihood computation
        double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayM = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayX = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayY = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        matchMetricArray[0][0]= END_GAP_COST;//Double.NEGATIVE_INFINITY;

        for (int i=1; i < X_METRIC_LENGTH; i++) {
            //initialize first column
            matchMetricArray[i][0]  = Double.NEGATIVE_INFINITY;
            YMetricArray[i][0]      = Double.NEGATIVE_INFINITY;
            XMetricArray[i][0]      = END_GAP_COST*(i);//logGapOpenProbability + (i-1)*logGapContinuationProbability;

            bestActionArrayX[i][0] = bestActionArrayY[i][0] = bestActionArrayM[i][0] = UP_GOTO_X;
        }

        for (int j=1; j < Y_METRIC_LENGTH; j++) {
            // initialize first row
            matchMetricArray[0][j]  = Double.NEGATIVE_INFINITY;
            XMetricArray[0][j]      = Double.NEGATIVE_INFINITY;
            YMetricArray[0][j]      = END_GAP_COST*(j);//logGapOpenProbability + (j-1) * logGapContinuationProbability;

            bestActionArrayY[0][j] = bestActionArrayM[0][j] = bestActionArrayX[0][j] = LEFT_GOTO_Y;
        }

        for (int indI=1; indI < X_METRIC_LENGTH; indI++) {
            for (int indJ=1; indJ < Y_METRIC_LENGTH; indJ++) {

                byte x = readBases[indI-1];
                byte y = haplotypeBases[indJ-1];
                byte qual = readQuals[indI-1];

                double bestMetric = 0.0;
                int bestMetricIdx = 0;

                // compute metric for match/mismatch
                // workaround for reads whose bases quality = 0,
                if (qual < 1)
                    qual = 1;

                if (qual > MAX_CACHED_QUAL)
                    qual = MAX_CACHED_QUAL;

                double pBaseRead =  (x == y)? baseMatchArray[(int)qual]:baseMismatchArray[(int)qual];


                double[] metrics = new double[3];

                // update match array
                metrics[MATCH_OFFSET] = matchMetricArray[indI-1][indJ-1] + pBaseRead;
                metrics[X_OFFSET] = XMetricArray[indI-1][indJ-1] + pBaseRead;
                metrics[Y_OFFSET] = YMetricArray[indI-1][indJ-1] + pBaseRead;

                if (doViterbi) {
                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(metrics);

                matchMetricArray[indI][indJ] = bestMetric;
                bestActionArrayM[indI][indJ] = ACTIONS_M[bestMetricIdx];

                // update X array
                 // State X(i,j): X(1:i) aligned to a gap in Y(1:j).
                // When in last column of X, ie X(1:i) aligned to full Y, we don't want to penalize gaps
                Pair<Double,Double> p = getGapPenalties(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, X_OFFSET, currentGOP, currentGCP);

                metrics[MATCH_OFFSET] = matchMetricArray[indI-1][indJ] + p.first;
                metrics[X_OFFSET] = XMetricArray[indI-1][indJ] + p.second;
                metrics[Y_OFFSET] = Double.NEGATIVE_INFINITY; //YMetricArray[indI-1][indJ] + logGapOpenProbability;

                if (doViterbi) {
                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(metrics);

                XMetricArray[indI][indJ] = bestMetric;
                bestActionArrayX[indI][indJ] = ACTIONS_X[bestMetricIdx];

                // update Y array
                p = getGapPenalties(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, Y_OFFSET, currentGOP, currentGCP);

                metrics[MATCH_OFFSET] = matchMetricArray[indI][indJ-1] + p.first;
                metrics[X_OFFSET] = Double.NEGATIVE_INFINITY; //XMetricArray[indI][indJ-1] + logGapOpenProbability;
                metrics[Y_OFFSET] = YMetricArray[indI][indJ-1] + p.second;

                if (doViterbi) {
                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(metrics);

                YMetricArray[indI][indJ] = bestMetric;
                bestActionArrayY[indI][indJ] = ACTIONS_Y[bestMetricIdx];



            }
        }

        double bestMetric;
        double metrics[] = new double[3];
        int bestTable=0, bestI=X_METRIC_LENGTH - 1, bestJ=Y_METRIC_LENGTH - 1;
        metrics[MATCH_OFFSET] = matchMetricArray[bestI][bestJ];
        metrics[X_OFFSET] = XMetricArray[bestI][bestJ];
        metrics[Y_OFFSET] = YMetricArray[bestI][bestJ];
        if (doViterbi) {
            bestTable = MathUtils.maxElementIndex(metrics);
            bestMetric = metrics[bestTable];
        }
        else
            bestMetric = MathUtils.softMax(metrics);

        // Do traceback (needed only for debugging!)
        if (DEBUG && doViterbi) {

            int bestAction;
            int i = bestI;
            int j = bestJ;


            System.out.println("Affine gap NW");


            String haplotypeString = new String (haplotypeBases);
            String readString = new String(readBases);


            while (i >0 || j >0) {
                if (bestTable == X_OFFSET) {
                    // insert gap in Y
                    haplotypeString = haplotypeString.substring(0,j)+"-"+haplotypeString.substring(j);
                    bestAction = bestActionArrayX[i][j];
                }
                else if (bestTable == Y_OFFSET) {
                    readString = readString.substring(0,i)+"-"+readString.substring(i);
                    bestAction = bestActionArrayY[i][j];

                }
                else {
                    bestAction = bestActionArrayM[i][j];
                }
                System.out.print(bestAction);


                // bestAction contains action to take at next step
                // encoding of bestAction: upper 2 bits = direction, lower 2 bits = next table

                // bestTable and nextDirection for next step
                bestTable = bestAction & 0x3;
                int nextDirection = bestAction >> 2;
                if (nextDirection == UP) {
                    i--;
                } else if (nextDirection == LEFT) {
                    j--;
                } else { //  if (nextDirection == DIAG)
                    i--; j--;
                }

             }




            System.out.println("\nAlignment: ");
            System.out.println("R:"+readString);
            System.out.println("H:"+haplotypeString);
            System.out.println();


        }
        if (DEBUG)
            System.out.format("Likelihood: %5.4f\n", bestMetric);

        return bestMetric;

    }

    private void fillGapProbabilities(int hIndex, int[] hrunProfile,
                                      double[][] contextLogGapOpenProbabilities, double[][] contextLogGapContinuationProbabilities) {
        // fill based on lookup table
        for (int i = 0; i < hrunProfile.length; i++) {
            if (hrunProfile[i] >= MAX_HRUN_GAP_IDX) {
                contextLogGapOpenProbabilities[hIndex][i] = GAP_OPEN_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
                contextLogGapContinuationProbabilities[hIndex][i] = GAP_CONT_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
            }
            else {
                contextLogGapOpenProbabilities[hIndex][i] = GAP_OPEN_PROB_TABLE[hrunProfile[i]];
                contextLogGapContinuationProbabilities[hIndex][i] = GAP_CONT_PROB_TABLE[hrunProfile[i]];
            }                
        }
    }
    public synchronized double[][] computeReadHaplotypeLikelihoods(ReadBackedPileup pileup, List<Haplotype> haplotypesInVC,
                                                                   ReferenceContext ref, int haplotypeSize, int eventLength){
        double[][] haplotypeLikehoodMatrix = new double[haplotypesInVC.size()][haplotypesInVC.size()];
        double readLikelihoods[][] = new double[pileup.getReads().size()][haplotypesInVC.size()];
        int readIdx=0;

        double[][] contextLogGapOpenProbabilities = null;
        double[][] contextLogGapContinuationProbabilities = null;


        if (DEBUG) {
            System.out.println("Reference bases:");
            System.out.println(new String(ref.getBases()));
        }

        if (doContextDependentPenalties)   {
            // will context dependent probabilities based on homopolymet run. Probabilities are filled based on total complete haplotypes.

            for (int j=0; j < haplotypesInVC.size(); j++) {
                Haplotype haplotype = haplotypesInVC.get(j);
                byte[] haplotypeBases = haplotype.getBasesAsBytes();
                if (contextLogGapOpenProbabilities == null) {
                    contextLogGapOpenProbabilities = new double[haplotypesInVC.size()][haplotypeBases.length];
                    contextLogGapContinuationProbabilities = new double[haplotypesInVC.size()][haplotypeBases.length];
                }
                // get homopolymer length profile for current haplotype
                int[] hrunProfile = new int[haplotypeBases.length];
                getContextHomopolymerLength(haplotypeBases,hrunProfile);
                if (DEBUG) {
                    System.out.println("Haplotype bases:");
                    System.out.println(new String(haplotypeBases));
                    for (int i=0; i < hrunProfile.length; i++)
                        System.out.format("%d",hrunProfile[i]);
                    System.out.println();
                }
                fillGapProbabilities(j, hrunProfile, contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities);

            }
        }
        for (SAMRecord read : pileup.getReads()) {
            if(ReadUtils.is454Read(read)) {
                continue;
            }

            // for each read/haplotype combination, compute likelihoods, ie -10*log10(Pr(R | Hi))
            // = sum_j(-10*log10(Pr(R_j | Hi) since reads are assumed to be independent
            if (DEBUG)
                System.out.format("\n\nStarting read:%s S:%d US:%d E:%d UE:%d C:%s\n",read.getReadName(),
                        read.getAlignmentStart(),
                        read.getUnclippedStart(), read.getAlignmentEnd(), read.getUnclippedEnd(),
                        read.getCigarString());


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

            // default: discard soft-clipped bases
            numStartSoftClippedBases = read.getAlignmentStart()- read.getUnclippedStart();
            numEndSoftClippedBases = read.getUnclippedEnd()- read.getAlignmentEnd();

            /*if (eventLength > 0) */
            {
                if ((read.getAlignmentStart()>=eventStartPos-eventLength && read.getAlignmentStart() <= eventStartPos+1) ||
                    (read.getAlignmentEnd() >= eventStartPos && read.getAlignmentEnd() <= eventStartPos + eventLength)) {
                    numStartSoftClippedBases = 0;
                    numEndSoftClippedBases = 0;
                }
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
            if (DEBUG)
                System.out.format("numStartClippedBases: %d numEndClippedBases: %d WinStart:%d WinStop:%d start: %d stop: %d readLength: %d\n",
                        numStartClippedBases, numEndClippedBases, ref.getWindow().getStart(), ref.getWindow().getStop(), start, stop, read.getReadLength());



            if (numStartClippedBases + numEndClippedBases >= unclippedReadBases.length) {
                if (DEBUG)
                    System.out.println("BAD READ!!");

                for (int j=0; j < haplotypesInVC.size(); j++) {
                    readLikelihoods[readIdx][j]= 0;
                }
            }
            else {
                byte[] readBases = Arrays.copyOfRange(unclippedReadBases,numStartClippedBases,
                        unclippedReadBases.length-numEndClippedBases);

                byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals,numStartClippedBases,
                        unclippedReadBases.length-numEndClippedBases);


                if (DEBUG) {
                    System.out.println("Read bases:");
                    System.out.println(new String(readBases));
                }


                // start and stop have indices into

                for (int j=0; j < haplotypesInVC.size(); j++) {
                    Haplotype haplotype = haplotypesInVC.get(j);

                    if (stop > haplotype.getStopPosition())
                        stop = haplotype.getStopPosition();

                    if (start < haplotype.getStartPosition())
                        start = haplotype.getStartPosition();

                    // cut haplotype bases
                    long indStart = start - haplotype.getStartPosition();
                    long indStop =  stop - haplotype.getStartPosition();

                    byte[] haplotypeBases = Arrays.copyOfRange(haplotype.getBasesAsBytes(),
                            (int)indStart, (int)indStop);

                    if (DEBUG) {
                        System.out.println("Haplotype to test:");
                        System.out.println(new String(haplotypeBases));
                    }

                    if (useAffineGapModel) {

                        double[] currentContextGOP = null;
                        double[] currentContextGCP = null;

                        if (doContextDependentPenalties) {
                            currentContextGOP = Arrays.copyOfRange(contextLogGapOpenProbabilities[j], (int)indStart, (int)indStop);
                            currentContextGCP = Arrays.copyOfRange(contextLogGapContinuationProbabilities[j], (int)indStart, (int)indStop);
                        }
                        readLikelihoods[readIdx][j]= computeReadLikelihoodGivenHaplotypeAffineGaps(haplotypeBases, readBases, readQuals, currentContextGOP, currentContextGCP);

                    }
                    else
                        readLikelihoods[readIdx][j]= computeReadLikelihoodGivenHaplotype(haplotypeBases, readBases, readQuals);


                }
            }

            readIdx++;
        }

        if (DEBUG) {
            System.out.println("\nLikelihood summary");
            for (readIdx=0; readIdx < pileup.getReads().size(); readIdx++) {
                System.out.format("Read Index: %d ",readIdx);
                for (int i=0; i < readLikelihoods[readIdx].length; i++)
                    System.out.format("L%d: %f ",i,readLikelihoods[readIdx][i]);
                System.out.println();
            }

        }
        for (int i=0; i < haplotypesInVC.size(); i++) {
            for (int j=i; j < haplotypesInVC.size(); j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                // L(Hi, Hj) = sum_reads ( Pr(R|Hi)/2 + Pr(R|Hj)/2)
                //readLikelihoods[k][j] has log10(Pr(R_k) | H[j] )
                double[] readLikelihood = new double[2]; // diploid sample
                for (readIdx=0; readIdx < pileup.getReads().size(); readIdx++) {

                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    if (Double.isInfinite(readLikelihoods[readIdx][i]) || Double.isInfinite(readLikelihoods[readIdx][j]))
                        continue;
                    haplotypeLikehoodMatrix[i][j] += ( MathUtils.softMax(readLikelihoods[readIdx][i],
                            readLikelihoods[readIdx][j]) + LOG_ONE_HALF);

                }


            }
        }

        return haplotypeLikehoodMatrix;

    }

    public static double[] getHaplotypeLikelihoods(double[][] haplotypeLikehoodMatrix) {
        int hSize = haplotypeLikehoodMatrix.length;
        double[] genotypeLikelihoods = new double[hSize*(hSize+1)/2];

        int k=0;
        double maxElement = Double.NEGATIVE_INFINITY;
        for (int i=0; i < hSize; i++) {
            for (int j=i; j < hSize; j++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
                if (haplotypeLikehoodMatrix[i][j] > maxElement)
                    maxElement = haplotypeLikehoodMatrix[i][j];
            }
        }

        // renormalize
        for (int i=0; i < genotypeLikelihoods.length; i++)
            genotypeLikelihoods[i] -= maxElement;

        return genotypeLikelihoods;

    }

 
}
