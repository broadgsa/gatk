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

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.Haplotype;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: Aug 24, 2010
 * Time: 1:31:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class HaplotypeIndelErrorModel {

    private final int maxReadDeletionLength; // maximum length of deletion on a read
    private final double noDeletionProbability; // alpha for geometric distribution for deletion length
    private final int haplotypeSize;
    private final int maxReadLength;
    private final int PATH_METRIC_TABLE_LENGTH;
    private final int RIGHT_ALIGN_INDEX;
    private final int LEFT_ALIGN_INDEX;

    private static final double INFINITE = 10000000000000.0;

    private double deletionErrorProbabilities[];

    private double pathMetricArray[];
    private int bestStateIndexArray[][];

    private final double logOneMinusInsertionStartProbability;
    private final double logInsertionStartProbability;
    private final double logInsertionEndProbability;
    private final double logOneMinusInsertionEndProbability;

    public  HaplotypeIndelErrorModel(int mrdl, double insStart, double insEnd, double alpha, int haplotypeSize, int maxReadLength) {
        this.maxReadDeletionLength = mrdl;
        this.noDeletionProbability = 1-alpha;
        this.haplotypeSize = haplotypeSize;
        this.maxReadLength = maxReadLength;

        PATH_METRIC_TABLE_LENGTH = haplotypeSize+2;
        RIGHT_ALIGN_INDEX = PATH_METRIC_TABLE_LENGTH-1;
        LEFT_ALIGN_INDEX = 0;

        logOneMinusInsertionStartProbability = probToQual(1-insStart);
        logInsertionStartProbability = probToQual(insStart);
        logInsertionEndProbability = probToQual(insEnd);
        logOneMinusInsertionEndProbability = probToQual(1-insEnd);


        // initialize path metric and traceback memories for Viterbi computation
        pathMetricArray = new double[2*PATH_METRIC_TABLE_LENGTH];
        bestStateIndexArray = new int[maxReadLength][2*(PATH_METRIC_TABLE_LENGTH)];

        for (int k=0; k < 2*PATH_METRIC_TABLE_LENGTH; k++)
            pathMetricArray[k] = 0;

        // fill in probability for read deletion of length k = C*exp(k-1)
        double prob = 1.0;
        double sumProb = 0.0;

        deletionErrorProbabilities = new double[maxReadDeletionLength+1];

        deletionErrorProbabilities[1] = noDeletionProbability;
        for (int k=2; k <= maxReadDeletionLength; k++) {
            deletionErrorProbabilities[k] = prob;
            sumProb = sumProb + prob;
            prob = prob*Math.exp(-1);
        }

        // now get everything in log domain, set normalizing constant so that probabilities sum to one
        deletionErrorProbabilities[1] = probToQual(deletionErrorProbabilities[1]);
        for (int k=2; k <= maxReadDeletionLength; k++) {
            deletionErrorProbabilities[k] = probToQual((1-noDeletionProbability)*deletionErrorProbabilities[k]/sumProb);
        }



    }

    public static double probToQual(double prob) {
        // TODO: see if I can use QualityUtils version, right now I don't want rounding or byte conversion
        return -10.0*Math.log10(prob);
    }

    public double computeReadLikelihoodGivenHaplotype(Haplotype haplotype, SAMRecord read) {

        final long readStartPosition = read.getAlignmentStart();
        final long haplotypeStartPosition = haplotype.getStartPosition();
        final long readEndPosition = read.getUnclippedEnd();
        final long haplotypeEndPosition = haplotype.getStartPosition() + haplotypeSize-1;
        final byte[] readBases = read.getReadBases();
        final byte[] readQuals = read.getBaseQualities();

        double[] previousPathMetricArray = new double[2*PATH_METRIC_TABLE_LENGTH];
        int readStartIdx, initialIndexInHaplotype;

  /*      List<AlignmentBlock> b = read.getAlignmentBlocks();
        // TODO temp hack, dont take realigned reads for now
        if (b.size()>1) {
            rrc++;
            //System.out.format("found realigned read %d\n", rrc);
                return 0.0;
        
        }
    */
        // case 1: read started before haplotype start position, we need to iterate only over remainder of read
        if (readStartPosition < haplotypeStartPosition) {
            if (readEndPosition<= haplotypeEndPosition)
                readStartIdx = (int)(haplotypeStartPosition - (readEndPosition - read.getReadBases().length)-1);
            else
                readStartIdx = (int)(haplotypeStartPosition- readStartPosition);
            initialIndexInHaplotype = 0; // implicit above

        }
        else {
            readStartIdx = 0;
            initialIndexInHaplotype = (int)(readStartPosition-haplotypeStartPosition);
        }


        // initialize path metric array to hard-code certainty that initial position is aligned
        for (int k=0; k < 2*PATH_METRIC_TABLE_LENGTH; k++)
            pathMetricArray[k] = INFINITE;

        pathMetricArray[initialIndexInHaplotype] = 0;

        /*
        System.out.println(read.getReadName());

        System.out.println(haplotypeStartPosition);
         System.out.println(readStartPosition);
        System.out.println(readStartIdx);
        System.out.println(initialIndexInHaplotype);
          */
        
        // Update path metric computations based on branch metric (Add/Compare/Select operations)
        // do forward direction first, ie from anchor to end of read
        // outer loop
        for (int indR=readStartIdx; indR < readBases.length; indR++) {
            byte readBase = readBases[indR];
            byte readQual = readQuals[indR];

            System.arraycopy(pathMetricArray, 0, previousPathMetricArray, 0, 2*PATH_METRIC_TABLE_LENGTH);



            for (int indX=LEFT_ALIGN_INDEX; indX <= RIGHT_ALIGN_INDEX; indX++) {


                byte haplotypeBase;
                if (indX > LEFT_ALIGN_INDEX && indX <= haplotypeSize)
                    haplotypeBase = haplotype.getBasesAsBytes()[indX-1];
                else
                    haplotypeBase = readBase;

                updatePathMetrics(haplotypeBase, indX, indR, readBase, readQual, previousPathMetricArray, true);
            }
        }

        double[] forwardPathMetricArray = new double[2*PATH_METRIC_TABLE_LENGTH];
        System.arraycopy(pathMetricArray, 0, forwardPathMetricArray, 0, 2*PATH_METRIC_TABLE_LENGTH);
        double forwardMetric = MathUtils.arrayMin(pathMetricArray);

        // reinitialize path metric array to known position at start of read
        // initialize path metric array to hard-code certainty that initial position is aligned
        for (int k=0; k < 2*PATH_METRIC_TABLE_LENGTH; k++)
            pathMetricArray[k] = INFINITE;

        pathMetricArray[initialIndexInHaplotype] = 0;
        
        // do now backward direction (from anchor to start of read)
        // outer loop
        for (int indR=readStartIdx-1; indR >= 0; indR--) {
            byte readBase = readBases[indR];
            byte readQual = readQuals[indR];

            System.arraycopy(pathMetricArray, 0, previousPathMetricArray, 0, 2*PATH_METRIC_TABLE_LENGTH);

            for (int indX=LEFT_ALIGN_INDEX; indX <= RIGHT_ALIGN_INDEX; indX++) {
                byte haplotypeBase;
                if (indX > 0 && indX <= haplotypeSize)
                    haplotypeBase = haplotype.getBasesAsBytes()[indX-1];
                else
                    haplotypeBase = readBase;

                updatePathMetrics(haplotypeBase, indX, indR, readBase, readQual, previousPathMetricArray, false);

            }
        }


        // for debugging only: compute backtracking to find optimal route through trellis. Since I'm only interested
        // in log-likelihood of best state, this isn't really necessary.
/*
        int[] bestIndexArray = new int[readBases.length+1];
         

        int bestIndex = MathUtils.minElementIndex(forwardPathMetricArray);
        bestIndexArray[readBases.length] = bestIndex;

        for (int k=readBases.length-1; k>=0; k--) {
            bestIndex = bestStateIndexArray[k][bestIndex];
            bestIndexArray[k] = bestIndex;
        }
        

        System.out.println(read.getReadName());
        System.out.print("Haplotype:");

        for (int k=0; k <haplotype.getBasesAsBytes().length; k++) {
            System.out.format("%c ", haplotype.getBasesAsBytes()[k]);
        }
        System.out.println();

        System.out.print("Read bases: ");
        for (int k=0; k <readBases.length; k++) {
            System.out.format("%c ", readBases[k]);
        }
        System.out.println();

        System.out.print("Alignment: ");
        for (int k=0; k <readBases.length; k++) {
            System.out.format("%d ", bestIndexArray[k]);
        }
        System.out.println();

        */
        // now just take optimum along all path metrics: that's the log likelihood of best alignment
        double backwardMetric = MathUtils.arrayMin(pathMetricArray);
        return forwardMetric + backwardMetric;

    }

    private void updatePathMetrics(byte haplotypeBase, int indX, int indR, byte readBase, byte readQual,
                                   double[] previousPathMetricArray, boolean isForward) {

        double bmetric;

        // Pr(Rb', Xb',Ib'| Xb,Ib=0) is populated starting on Xb , =
        // Pr(Rb'|Xb',Ib') Pr(Xb',Ib'|Xb,Ib) =
        // Pr(Rb'|Xb',Ib') Pr(Xb'|Ib',Xb,Ib) Pr(Ib'|Ib)

        // start with case Ib=0, Ib'=0: no insertion
        double bestMetric = INFINITE;
        int bestMetricIndex = 0;
        double pBaseRead = getProbabilityOfReadBaseGivenXandI(haplotypeBase, readBase, readQual, indX, 0);

        if (isForward) {
            if (indX==LEFT_ALIGN_INDEX){
                // there can't be a transition into Xb=0 
                // special case: only one incoming path
                // record best path metric
//                pathMetricArray[indX] = INFINITE;
//                bestStateIndexArray[indR][indX] = 0;
//                bestMetric = INFINITE;
//                bestMetricIndex = 0;
            }
            else {
                for (int indXold = indX-1; indXold >= indX-this.maxReadDeletionLength; indXold--){

                    if (indXold < 0)
                        break;

                    // fetch path metric and add branch metric
                    bmetric = previousPathMetricArray[indXold] + deletionErrorProbabilities[indX-indXold] +
                            logOneMinusInsertionStartProbability + pBaseRead;

                    if (bmetric < bestMetric) {
                        bestMetric = bmetric;
                        bestMetricIndex = indXold;
                    }
                }
                if (indX == RIGHT_ALIGN_INDEX)  {
                    // special treatment for R->R case  (read aligns to the right of haplotype) when Xold = X = R
                    bmetric = previousPathMetricArray[indX] + pBaseRead;

                    if (bmetric < bestMetric) {
                        bestMetric = bmetric;
                        bestMetricIndex = indX;
                    }


                }
            }            
        }
        else {
            for (int indXold = indX+1; indXold <=this.maxReadDeletionLength + indX; indXold++){
                if (indXold > this.haplotypeSize+1)
                    break;

                // fetch path metric and add branch metric
                bmetric = previousPathMetricArray[indXold] + deletionErrorProbabilities[indXold-indX] +
                        logOneMinusInsertionStartProbability + pBaseRead;

                if (bmetric < bestMetric) {
                    bestMetric = bmetric;
                    bestMetricIndex = indXold;
                }
                if (indX == LEFT_ALIGN_INDEX)  {
                    // special treatment for R->R case  (read aligns to the right of haplotype) when Xold = X = R
                    bmetric = previousPathMetricArray[indX] + pBaseRead;

                    if (bmetric < bestMetric) {
                        bestMetric = bmetric;
                        bestMetricIndex = indX;
                    }


                }
            }

        }
        // now Pr(Xb',Ib'=0| Xb,Ib=1): an insertion ended. Only case if Xb'=Xb+1 (walk diag down in path array).

        if (indX > LEFT_ALIGN_INDEX && indX < RIGHT_ALIGN_INDEX)  {
            if (isForward)
                bmetric = previousPathMetricArray[indX-1+PATH_METRIC_TABLE_LENGTH] + logInsertionEndProbability + pBaseRead;
            else
                bmetric = previousPathMetricArray[indX+1+PATH_METRIC_TABLE_LENGTH] + logInsertionEndProbability + pBaseRead;

            if (bmetric < bestMetric) {
                bestMetric = bmetric;
                if (isForward)
                    bestMetricIndex = indX-1 + PATH_METRIC_TABLE_LENGTH;
                else
                    bestMetricIndex = indX+1 + PATH_METRIC_TABLE_LENGTH;

            }
        }
        // record best path metric
        pathMetricArray[indX] = bestMetric;
        bestStateIndexArray[indR][indX] = bestMetricIndex;

        // now Pr(Xb', Ib'=1 | Xb, Ib=0) : an insertion started: (walk right in path array)
        bestMetric = INFINITE;
        pBaseRead = getProbabilityOfReadBaseGivenXandI(haplotypeBase, readBase, readQual, indX, 1);

        bmetric = previousPathMetricArray[indX] + logInsertionStartProbability + pBaseRead;
        if (bmetric < bestMetric) {     //if superfluous, could clean up
            bestMetric = bmetric;
            bestMetricIndex = indX;
        }

        // final case: Pr(Xb', Ib'=1 | Xb, Ib=1): insertion continues, also walk right in path array
        if (indX > LEFT_ALIGN_INDEX && indX < RIGHT_ALIGN_INDEX)  {
            bmetric = previousPathMetricArray[indX+PATH_METRIC_TABLE_LENGTH] + logOneMinusInsertionEndProbability + pBaseRead;
            if (bmetric < bestMetric) {     //if superfluous, could clean up
                bestMetric = bmetric;
                bestMetricIndex = indX+PATH_METRIC_TABLE_LENGTH;
            }

            // record best path metric
            pathMetricArray[indX+PATH_METRIC_TABLE_LENGTH] = bestMetric;
            bestStateIndexArray[indR][indX+PATH_METRIC_TABLE_LENGTH] = bestMetricIndex;
        }
    }

    private double getProbabilityOfReadBaseGivenXandI(byte haplotypeBase, byte readBase, byte readQual, int indX, int indI) {


        if (indX == this.haplotypeSize  || indX == 0) {
            // X = L or R
            double baseProb = QualityUtils.qualToProb(readQual);
            return probToQual(baseProb);
        }

        if (haplotypeBase == readBase) {
            double baseProb = QualityUtils.qualToProb(readQual);
            return probToQual(baseProb);
        }
        else {
            return (double)(readQual);
        }
        // TODO - optimize for speed by avoiding so many log conversions, can use a LUT

    }

}
