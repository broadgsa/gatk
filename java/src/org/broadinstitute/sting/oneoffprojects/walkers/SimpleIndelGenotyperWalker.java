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

package org.broadinstitute.sting.oneoffprojects.walkers;

import net.sf.samtools.SAMRecord;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.indels.HaplotypeIndelErrorModel;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.CircularArray;
import org.broadinstitute.sting.utils.genotype.Haplotype;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.*;

import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;


@Reference(window=@Window(start=-30,stop=200))
@Allows({DataSource.READS, DataSource.REFERENCE})
public class SimpleIndelGenotyperWalker extends RefWalker<Integer,Integer> {
    @Output
    PrintStream out;
    @Argument(fullName="maxReadDeletionLength",shortName="maxReadDeletionLength",doc="Max deletion length allowed when aligning reads to candidate haplotypes.",required=false)
    int maxReadDeletionLength = 3;

    @Argument(fullName="insertionStartProbability",shortName="insertionStartProbability",doc="Assess only sites with coverage at or below the specified value.",required=false)
    double insertionStartProbability = 0.01;

    @Argument(fullName="insertionEndProbability",shortName="insertionEndProbability",doc="Assess only sites with coverage at or below the specified value.",required=false)
    double insertionEndProbability = 0.5;

    @Argument(fullName="alphaDeletionProbability",shortName="alphaDeletionProbability",doc="Assess only sites with coverage at or below the specified value.",required=false)
    double alphaDeletionProbability = 0.01;


    @Argument(fullName="haplotypeSize",shortName="hsize",doc="Size of haplotypes to evaluate calls.",required=true)
    int HAPLOTYPE_SIZE = 40;


    @Argument(fullName="indelHeterozygozity",shortName="indhet",doc="Indel Heterozygosity (assumed 1/6500 from empirical human data)",required=false)
    double indelHeterozygosity = (double)1.0/6500.0;

    @Argument(fullName="sampleName",shortName="sample",doc="Sample name to evaluate genotypes in",required=true)
    String sampleName;

    @Argument(fullName="doSimpleCalculationModel",shortName="simple",doc="Use Simple Calculation Model for Pr(Reads | Haplotype)",required=false)
    boolean doSimple = false;

    @Argument(fullName="enableDebugOutput",shortName="debugout",doc="Output debug data",required=false)
    boolean DEBUG = false;

    @Argument(fullName="doTrio",shortName="trio",doc="Output 1KG CEU trio genotype data (sampleName ignored)",required=false)
    boolean doTrio = false;

    @Argument(fullName="useFlatPriors",shortName="flat",doc="If present, use flat priors on haplotypes",required=false)
    boolean useFlatPriors = false;

    @Override
    public boolean generateExtendedEvents() { return true; }

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }


    private HaplotypeIndelErrorModel model;

    private static final int MAX_READ_LENGTH = 1000; // TODO- make this dynamic




    @Override
    public void initialize() {
        model = new HaplotypeIndelErrorModel(maxReadDeletionLength, insertionStartProbability,
                insertionEndProbability, alphaDeletionProbability, HAPLOTYPE_SIZE, MAX_READ_LENGTH, doSimple, DEBUG);

    }



    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        if (!context.hasBasePileup())
            return 0;

        VariantContext vc = tracker.getVariantContext(ref, "indels", null, context.getLocation(), true);
        // ignore places where we don't have a variant
        if ( vc == null )
            return 0;


        if (!vc.isIndel())
            return 0;


        List<Haplotype> haplotypesInVC = Haplotype.makeHaplotypeListFromVariantContextAlleles( vc, ref, HAPLOTYPE_SIZE);

        // combine likelihoods for all possible haplotype pair (n*(n-1)/2 combinations)

        double[][] haplotypeLikehoodMatrix = new double[haplotypesInVC.size()][haplotypesInVC.size()];
        double bestLikelihood = 1e13;
        // for given allele haplotype, compute likelihood of read pileup given haplotype
        ReadBackedPileup pileup =  context.getBasePileup().getPileupWithoutMappingQualityZeroReads();


        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.
        // In general, we'll assume: even spread of indels throughout genome (not true, but simplifying assumption),
        // and memoryless spread (i.e. probability that an indel lies in an interval A is independent of probability of
        // another indel lying in interval B iff A and B don't overlap), then we can approximate inter-indel distances
        // by an exponential distribution of mean 1/theta (theta = heterozygozity), and the number of indels on an interval
        // of size L is Poisson-distributed with parameter lambda = theta*L.

        // Since typically, for small haplotype sizes and human heterozygozity, lambda will be <<1, we'll further approximate it
        // by assuming that only one indel can happen in a particular interval, with Pr(indel present) = lambda*exp(-lambda), and
        // pr(no indel) = 1-lambda*exp(-lambda) ~= exp(-lambda) for small lambda.

        // We also assume that a deletion is equally likely as an insertion (empirical observation, see e.g. Mills et al, Genome Research 2006)
        // and we assume the following frequency spectrum for indel sizes Pr(event Length = L)= K*abs(L)^(-1.89)*10^(-0.015*abs(L)),
        // taking positive L = insertions, negative L = deletions. K turns out to be about 1.5716 for probabilities to sum to one.
        // so -10*log10(Pr event Length = L) =-10*log10(K)+ 18.9*log10(abs(L)) + 0.15*abs(L).
        // Hence, Pr(observe event size = L in interval) ~ Pr(observe event L | event present) Pr (event present in interval)
        // and -10*log10(above) = -10*log10(K)+ 18.9*log10(abs(L)) + 0.15*abs(L) - 10*log10(theta*L), and we ignore terms that would be
        // added to ref hypothesis.
        // Equation above is prior model.
        int eventLength = vc.getReference().getBaseString().length() - vc.getAlternateAllele(0).getBaseString().length(); // assume only one alt allele for now
        if (eventLength<0)
            eventLength = - eventLength;

        if (!useFlatPriors) {

            double lambda = (double)HAPLOTYPE_SIZE * indelHeterozygosity;
            double altPrior = HaplotypeIndelErrorModel.probToQual(lambda)-HaplotypeIndelErrorModel.probToQual(eventLength)*1.89 + 0.15*eventLength
                    + HaplotypeIndelErrorModel.probToQual(1.5716)+ HaplotypeIndelErrorModel.probToQual(0.5);

            haplotypeLikehoodMatrix[0][1] = altPrior;
            haplotypeLikehoodMatrix[1][1] = 2*altPrior;
        }

        int bestIndexI =-1, bestIndexJ=-1;
        double callConfidence = 0.0;

        if (pileup.getReads().size() > 0) {
            double readLikelihoods[][] = new double[pileup.getReads().size()][haplotypesInVC.size()];
            int i=0;
            for (SAMRecord read : pileup.getReads()) {
                // for each read/haplotype combination, compute likelihoods, ie -10*log10(Pr(R | Hi))
                // = sum_j(-10*log10(Pr(R_j | Hi) since reads are assumed to be independent
                for (int j=0; j < haplotypesInVC.size(); j++) {
                    readLikelihoods[i][j]= model.computeReadLikelihoodGivenHaplotype(haplotypesInVC.get(j), read, vc, eventLength);
                    if (DEBUG) {
                        System.out.print(read.getReadName()+" ");

                        System.out.format("%d %d S:%d US:%d E:%d UE:%d C:%s %3.4f\n",i, j, read.getAlignmentStart(),
                                read.getUnclippedStart(), read.getAlignmentEnd(), read.getUnclippedEnd(),
                                read.getCigarString(), readLikelihoods[i][j]);
                    }

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

                        double probRGivenHPair = MathUtils.sumLog10(readLikelihood)/2;
                        haplotypeLikehoodMatrix[i][j] += HaplotypeIndelErrorModel.probToQual(probRGivenHPair);

                    }


                    if (haplotypeLikehoodMatrix[i][j] < bestLikelihood) {
                        bestIndexI = i;
                        bestIndexJ = j;
                        callConfidence = bestLikelihood - haplotypeLikehoodMatrix[i][j];
                        bestLikelihood = haplotypeLikehoodMatrix[i][j];
                    }
                }
            }
        }

        //we say that most likely genotype at site is index (i,j) maximizing likelihood matrix
        String type;
        if (vc.isDeletion())
            type = "DEL";
        else if (vc.isInsertion())
            type = "INS";
        else
            type = "OTH";

        type += vc.getIndelLengths().toString();


        out.format("%s %d %s ",vc.getChr(), vc.getStart(), type);

        if (doTrio) {
            Genotype originalGenotype = vc.getGenotype("NA12878");
            Genotype dadGenotype = vc.getGenotype("NA12891");
            Genotype momGenotype = vc.getGenotype("NA12892");

            String dadG, momG, oldG, newG;
            oldG = getGenotypeString(originalGenotype);
            dadG = getGenotypeString(dadGenotype);
            momG = getGenotypeString(momGenotype);

            int x = bestIndexI+bestIndexJ;
            if (x == 0)
                newG = "HOMREF";
            else if (x == 1)
                newG = "HET";
            else if (x == 2)
                newG = "HOMVAR";
            else if (x < 0)
                newG = "NOCALL";
            else
                newG = "OTHER";

            out.format("NewG %s OldG %s DadG %s MomG %s\n", newG, oldG, dadG, momG);
        }
        else {
            Genotype originalGenotype = vc.getGenotype(sampleName);
            String oldG, newG;
            oldG = getGenotypeString(originalGenotype);
            int x = bestIndexI+bestIndexJ;
            if (x == 0)
                newG = "HOMREF";
            else if (x == 1)
                newG = "HET";
            else if (x == 2)
                newG = "HOMVAR";
            else if (x < 0)
                newG = "NOCALL";
            else
                newG = "OTHER";

            out.format("NewG %s OldG %s\n", newG, oldG);

        }

/*
        if ( context.hasExtendedEventPileup() ) {

            // if we got indels at current position:

            ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup().getPileupWithoutMappingQualityZeroReads();
            if ( pileup.size() < MIN_COVERAGE ) return 0;

            if ( pileup.getNumberOfDeletions() + pileup.getNumberOfInsertions() > MAX_INDELS ) {
                // we got too many indel events. Maybe it's even a true event, and what we are looking for are
                // errors rather than true calls. Hence, we do not need these indels. We have to 1) discard
                // all remaining indels from the buffer: if they are still in the buffer, they are too close
                // to the current position; and 2) make sure that the next position at which we attempt to count again is
                // sufficiently far *after* the current position.
      //                      System.out.println("Non countable indel event at "+pileup.getLocation());
                countableIndelBuffer.clear();
                coverageBuffer.clear(); // we do not want to count observations (read bases) around non-countable indel as well
                skipToLoc = GenomeLocParser.createGenomeLoc(pileup.getLocation().getContigIndex(),pileup.getLocation().getStop()+pileup.getMaxDeletionLength()+MIN_DISTANCE+1);
 //                       System.out.println("Skip to "+skipToLoc);
            } else {
                // pileup does not contain too many indels, we need to store them in the buffer and count them later,
                // if a non-countable indel event(s) do not show up too soon:
                countableIndelBuffer.add(pileup);
            }
            return 0;
        }
         */


        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    private String getGenotypeString(Genotype genotype) {
        String oldG;
        if (genotype.isHomRef())
            oldG = "HOMREF";
        else if (genotype.isHet())
            oldG = "HET";
        else if (genotype.isHomVar())
            oldG = "HOMVAR";
        else
            oldG = "OTHER";

        return oldG;
    }
    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return value+sum;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void onTraversalDone(Integer result) {
    }

}

