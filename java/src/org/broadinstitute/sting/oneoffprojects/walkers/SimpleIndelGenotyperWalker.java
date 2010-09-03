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


@Reference(window=@Window(start=-10,stop=80))
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



    @Override
    public boolean generateExtendedEvents() { return true; }

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }


    private HaplotypeIndelErrorModel model;

    private static final int MAX_READ_LENGTH = 200; // TODO- make this dynamic




    @Override
    public void initialize() {
        model = new HaplotypeIndelErrorModel(maxReadDeletionLength, insertionStartProbability,
                insertionEndProbability, alphaDeletionProbability, HAPLOTYPE_SIZE, MAX_READ_LENGTH);

    }


    /*   private void countIndels(ReadBackedExtendedEventPileup p) {
         for ( ExtendedEventPileupElement pe : p.toExtendedIterable() ) {
             if ( ! pe.isIndel() ) continue;
             if ( pe.getEventLength() > MAX_LENGTH ) continue;
             if ( pe.isInsertion() ) insCounts[pe.getEventLength()-1]++;
             else delCounts[pe.getEventLength()-1]++;
         }
     }
    */
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


        int bestIndexI =0, bestIndexJ=0;
        for (int i=0; i < haplotypesInVC.size(); i++) {
            for (int j=i; j < haplotypesInVC.size(); j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                for (SAMRecord read : pileup.getReads()) {
                    // compute Pr(r | hi, hj) = 1/2*(Pr(r|hi) + Pr(r|hj)

                    double readLikelihood[] = new double[2];
                    readLikelihood[0]= -model.computeReadLikelihoodGivenHaplotype(haplotypesInVC.get(i), read)/10.0;
                    if (i != j)
                        readLikelihood[1] = -model.computeReadLikelihoodGivenHaplotype(haplotypesInVC.get(j), read)/10.0;
                    else
                        readLikelihood[1] = readLikelihood[0];

                    // likelihood is by convention -10*log10(pr(r|h))
                    // sumlog10 computes 10^x[0]+10^x[1]+...
                    double probRGivenHPair = MathUtils.sumLog10(readLikelihood)/2;
                    haplotypeLikehoodMatrix[i][j] += HaplotypeIndelErrorModel.probToQual(probRGivenHPair);

                }


                if (haplotypeLikehoodMatrix[i][j] < bestLikelihood) {
                    bestIndexI = i;
                    bestIndexJ = j;
                    bestLikelihood = haplotypeLikehoodMatrix[i][j];
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

        Genotype originalGenotype = vc.getGenotype("NA12878");

        String oldG, newG;
        if (originalGenotype.isHomRef())
            oldG = "HOMREF";
        else if (originalGenotype.isHet())
            oldG = "HET";
        else if (originalGenotype.isHomVar())
            oldG = "HOMVAR";
        else
            oldG = "OTHER";

        int x = bestIndexI+bestIndexJ;
        if (x == 0)
            newG = "HOMREF";
        else if (x == 1)
            newG = "HET";
        else if (x == 2)
            newG = "HOMVAR";
        else
            newG = "OTHER";

        out.format("NewG %s OldG %s\n", newG, oldG);


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

