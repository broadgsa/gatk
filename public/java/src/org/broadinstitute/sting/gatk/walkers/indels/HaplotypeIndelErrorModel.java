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

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;

public class HaplotypeIndelErrorModel {

    private final int maxReadDeletionLength; // maximum length of deletion on a read
    private final double noDeletionProbability; // alpha for geometric distribution for deletion length
    private final int haplotypeSize;
    
    private final int BASE_QUAL_THRESHOLD = 6;

    private final int PATH_METRIC_TABLE_LENGTH;
    private final int RIGHT_ALIGN_INDEX;
    private final int LEFT_ALIGN_INDEX;


    private double deletionErrorProbabilities[];

    private double pathMetricArray[][];
    private int bestStateIndexArray[][];

    private final double logOneMinusInsertionStartProbability;
    private final double logInsertionStartProbability;
    private final double logInsertionEndProbability;
    private final double logOneMinusInsertionEndProbability;

    private boolean DEBUG = false;
    private boolean doSimpleCalculationModel = false;

    private static final double QUAL_ONE_HALF = -10*Math.log10(0.5);

    private static final int MAX_CACHED_QUAL = 60;

    private static final double baseMatchArray[];
    private static final double baseMismatchArray[];

    static {
        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = QualityUtils.qualToProb((byte)k);


            baseMatchArray[k] =  probToQual(baseProb);
            baseMismatchArray[k] = (double)(k);
        }
    }

    public  HaplotypeIndelErrorModel(int mrdl, double insStart, double insEnd, double alpha, int haplotypeSize,
                                     boolean dosimple, boolean deb) {
        this(mrdl, insStart, insEnd, alpha, haplotypeSize);
        this.DEBUG = deb;
        this.doSimpleCalculationModel = dosimple;
    }
    public  HaplotypeIndelErrorModel(int mrdl, double insStart, double insEnd, double alpha, int haplotypeSize) {
        this.maxReadDeletionLength = mrdl;
        this.noDeletionProbability = 1-alpha;
        this.haplotypeSize = haplotypeSize;

        PATH_METRIC_TABLE_LENGTH = haplotypeSize+2;
        RIGHT_ALIGN_INDEX = PATH_METRIC_TABLE_LENGTH-1;
        LEFT_ALIGN_INDEX = 0;

        logOneMinusInsertionStartProbability = probToQual(1-insStart);
        logInsertionStartProbability = probToQual(insStart);
        logInsertionEndProbability = probToQual(insEnd);
        logOneMinusInsertionEndProbability = probToQual(1-insEnd);


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

        long numStartClippedBases = 0;
        long numEndClippedBases = 0;


        byte[] unclippedReadQuals = read.getBaseQualities();
        byte[] unclippedReadBases = read.getReadBases();


        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i=0; i < read.getReadLength(); i++) {
            if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                numStartClippedBases++;
            else
                break;

        }
        for (int i=read.getReadLength()-1; i >= 0; i-- ){
            if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                numEndClippedBases++;
            else
                break;
        }
        //System.out.format("numstart: %d numend: %d\n", numStartClippedBases, numEndClippedBases);
        if (numStartClippedBases + numEndClippedBases >= read.getReadLength()) {
            return 0;///Double.POSITIVE_INFINITY;
        }
        byte[] readBases = Arrays.copyOfRange(unclippedReadBases,(int)numStartClippedBases,
                (int)(read.getReadBases().length-numEndClippedBases));

        byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals,(int)numStartClippedBases,
                (int)(read.getReadBases().length-numEndClippedBases));


        int readLength = readBases.length;

        // initialize path metric and traceback memories for Viterbi computation
        pathMetricArray = new double[readLength+1][PATH_METRIC_TABLE_LENGTH];
        bestStateIndexArray = new int[readLength+1][PATH_METRIC_TABLE_LENGTH];

        for (int k=1; k < PATH_METRIC_TABLE_LENGTH; k++)
            pathMetricArray[0][k] = 0;

 /*

         if (doSimpleCalculationModel) {

            // No Viterbi algorithm - assume no sequencing indel artifacts,

            // so we can collapse computations and pr(read | haplotype) is just probability of observing overlap
            // of read with haplotype.
            int haplotypeIndex = initialIndexInHaplotype;
            double c =  0.0;//deletionErrorProbabilities[1] +logOneMinusInsertionStartProbability;
            // compute likelihood of portion of base to the left of the haplotype
            for (int indR=readStartIdx-1; indR >= 0; indR--) {
                byte readBase = readBases[indR];
                byte readQual = readQuals[indR];
                if (readQual <= 2)
                    continue;
                double pBaseRead = getProbabilityOfReadBaseGivenXandI((byte)0, readBase, readQual, LEFT_ALIGN_INDEX, 0);

                // pBaseRead has -10*log10(Prob(base[i]|haplotype[i])
                pRead += pBaseRead;

            }
            //System.out.format("\nSt: %d Pre-Likelihood:%f\n",readStartIdx, pRead);

            for (int indR=readStartIdx; indR < readBases.length; indR++) {
                byte readBase = readBases[indR];
                byte readQual = readQuals[indR];

                byte haplotypeBase;
                if (haplotypeIndex < RIGHT_ALIGN_INDEX)
                    haplotypeBase = haplotype.getBasesAsBytes()[haplotypeIndex];
                else
                    haplotypeBase = (byte)0; // dummy

                double pBaseRead = getProbabilityOfReadBaseGivenXandI(haplotypeBase, readBase, readQual, haplotypeIndex, 0);
                if (haplotypeBase != 0)
                    pBaseRead += c;

                // pBaseRead has -10*log10(Prob(base[i]|haplotype[i])
                if (readQual > 3)
                    pRead += pBaseRead;
                haplotypeIndex++;
                if (haplotypeIndex >= haplotype.getBasesAsBytes().length)
                    haplotypeIndex = RIGHT_ALIGN_INDEX;
                //System.out.format("H:%c R:%c RQ:%d HI:%d %4.5f %4.5f\n", haplotypeBase, readBase, (int)readQual, haplotypeIndex, pBaseRead, pRead);
             }
            //System.out.format("\nSt: %d Post-Likelihood:%f\n",readStartIdx, pRead);

            if (DEBUG) {
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
                System.out.format("\nLikelihood:%f\n",pRead);

            }

            if (read.getReadName().contains("106880")) {

                System.out.println("aca");

                System.out.println("Haplotype:");

                for (int k=initialIndexInHaplotype; k <haplotype.getBasesAsBytes().length; k++) {
                    System.out.format("%c ", haplotype.getBasesAsBytes()[k]);
                }
                System.out.println();

                System.out.println("Read bases: ");
                for (int k=readStartIdx; k <readBases.length; k++) {
                    System.out.format("%c ", readBases[k]);
                }

            }
            return pRead;

        }
                */


        // Update path metric computations based on branch metric (Add/Compare/Select operations)
        // do forward direction first, ie from anchor to end of read
        // outer loop
        for (int indR=0; indR < readLength; indR++) {
            byte readBase = readBases[indR];
            byte readQual = readQuals[indR];

             for (int indX=LEFT_ALIGN_INDEX; indX <= RIGHT_ALIGN_INDEX; indX++) {


                byte haplotypeBase;
                if (indX > LEFT_ALIGN_INDEX && indX < RIGHT_ALIGN_INDEX)
                    haplotypeBase = haplotype.getBasesAsBytes()[indX-1];
                else
                    haplotypeBase = readBase;

                updatePathMetrics(haplotypeBase, indX, indR, readBase, readQual);
            }
        }



        // for debugging only: compute backtracking to find optimal route through trellis. Since I'm only interested
        // in log-likelihood of best state, this isn't really necessary.
        double bestMetric = MathUtils.arrayMin(pathMetricArray[readLength]);

        if (DEBUG) {



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

            System.out.print("Read quals: ");
            for (int k=0; k <readQuals.length; k++) {
                System.out.format("%d ", (int)readQuals[k]);
            }
            System.out.println();

            // start from last position of read, go backwards to find optimal alignment
            int[] bestIndexArray = new int[readLength];
            int bestIndex = MathUtils.minElementIndex(pathMetricArray[readLength]);
            bestIndexArray[readLength-1] = bestIndex;

            for (int k=readLength-2; k>=0; k--) {
                bestIndex = bestStateIndexArray[k][bestIndex];
                bestIndexArray[k] = bestIndex;
            }
            
            System.out.print("Alignment: ");
            for (int k=0; k <readBases.length; k++) {
                System.out.format("%d ", bestIndexArray[k]);
            }
            System.out.println();
        }
        // now just take optimum along all path metrics: that's the log likelihood of best alignment
        if (DEBUG)
            System.out.format("Likelihood: %5.4f\n", bestMetric);
        return bestMetric;

    }

    private void updatePathMetrics(byte haplotypeBase, int indX, int indR, byte readBase, byte readQual) {

        double bmetric;

        double bestMetric = Double.POSITIVE_INFINITY;
        int bestMetricIndex = -1;

        // compute metric for match/mismatch
        double pBaseRead;

        // workaround for reads whose bases quality = 0,
        if (readQual < 1)
            readQual = 1;

        if (readQual > MAX_CACHED_QUAL)
            readQual = MAX_CACHED_QUAL;

        double pBaseMatch =  baseMatchArray[(int)readQual];
        double pBaseMismatch = baseMismatchArray[(int)readQual];

        if (haplotypeBase == readBase)
            pBaseRead =  pBaseMatch;
        else
            pBaseRead = pBaseMismatch;





        if (indX == LEFT_ALIGN_INDEX) {
            // special treatment for L->L case  (read aligns to the right/left of haplotype) when Xold = X = R
            bestMetric = pathMetricArray[indR][indX] + QUAL_ONE_HALF; //1/2 in log scale
            bestMetricIndex = indX;
        }
        else {
            for (int indXold = indX-1; indXold >= indX-this.maxReadDeletionLength; indXold--) {

                if (indXold < 0)
                    break;

                // fetch path metric and add branch metric
                bmetric = pathMetricArray[indR][indXold] + deletionErrorProbabilities[indX-indXold] + pBaseRead;

                if (indXold == indX-1) {
                    // special case for exact walk down diagonal: need to consider that an insertion may have ended
               //     bmetric += logInsertionEndProbability;
                } else {
                    bmetric += logOneMinusInsertionStartProbability;
                }

                if (bmetric < bestMetric) {
                    bestMetric = bmetric;
                    bestMetricIndex = indXold;
                }
            }
            // additional case: a walk right (ie indXold = indX). Can be because of an insertion in the middle of reads,
            // or we're aligned to right of read
            bmetric = pathMetricArray[indR][indX]+pBaseRead;
            if (indX < RIGHT_ALIGN_INDEX) {
                bmetric += logInsertionStartProbability + logOneMinusInsertionEndProbability;
            }
            else {
                // anything extra to do for R->R case?
                bmetric = pathMetricArray[indR][indX] + QUAL_ONE_HALF;
            }

            if (bmetric < bestMetric) {
                bestMetric = bmetric;
                bestMetricIndex = indX;
            }


        }

         // record best path metric
        pathMetricArray[indR+1][indX] = bestMetric;
        bestStateIndexArray[indR+1][indX] = bestMetricIndex;

    }

    public double[] computeReadHaplotypeLikelihoods(ReadBackedPileup pileup, HashMap<Allele,Haplotype> haplotypesInVC){
        double[][] haplotypeLikehoodMatrix = new double[haplotypesInVC.size()][haplotypesInVC.size()];
        double readLikelihoods[][] = new double[pileup.getReads().size()][haplotypesInVC.size()];
        int i=0;
        for (SAMRecord read : pileup.getReads()) {
            if(ReadUtils.is454Read(read)) {
                continue;
            }
            // for each read/haplotype combination, compute likelihoods, ie -10*log10(Pr(R | Hi))
            // = sum_j(-10*log10(Pr(R_j | Hi) since reads are assumed to be independent
            int j=0;
            for (Allele a: haplotypesInVC.keySet()) {
                readLikelihoods[i][j]= computeReadLikelihoodGivenHaplotype(haplotypesInVC.get(a), read);
                if (DEBUG) {
                    System.out.print(read.getReadName()+" ");

                    System.out.format("%d %d S:%d US:%d E:%d UE:%d C:%s %3.4f\n",i, j, read.getAlignmentStart(),
                            read.getUnclippedStart(), read.getAlignmentEnd(), read.getUnclippedEnd(),
                            read.getCigarString(), readLikelihoods[i][j]);
                }
                j++;
            }
            i++;
        }

        for (i=0; i < haplotypesInVC.size(); i++) {
            for (int j=i; j < haplotypesInVC.size(); j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                // L(Hi, Hj) = sum_reads ( Pr(R|Hi)/2 + Pr(R|Hj)/2)
                //readLikelihoods[k][j] has log10(Pr(R_k) | H[j] )
                double[] readLikelihood = new double[2]; // diploid sample
                for (int readIdx=0; readIdx < pileup.getReads().size(); readIdx++) {
                    readLikelihood[0] = -readLikelihoods[readIdx][i]/10;
                    readLikelihood[1] = -readLikelihoods[readIdx][j]/10;

                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+x0^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    // Second term is a constant added to both likelihoods so will be ignored
                    haplotypeLikehoodMatrix[i][j] += MathUtils.softMax(readLikelihood[0],
                            readLikelihood[1]);

                }


            }
        }

        return getHaplotypeLikelihoods(haplotypeLikehoodMatrix);

    }

    private double[] getHaplotypeLikelihoods(double[][] haplotypeLikehoodMatrix) {
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
