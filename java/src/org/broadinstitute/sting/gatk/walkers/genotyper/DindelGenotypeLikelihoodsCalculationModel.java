/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.indels.HaplotypeIndelErrorModel;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.Haplotype;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broad.tribble.util.variantcontext.Allele;

import java.util.*;

public class DindelGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {
    private final int maxReadDeletionLength = 3;
    private final double insertionStartProbability = 1e-3;
    private final double insertionEndProbability = 0.5;
    private final double alphaDeletionProbability = 1e-3;
    private final int HAPLOTYPE_SIZE = 80;

    // todo - the following  need to be exposed for command line argument control
    private final double indelHeterozygosity = 1.0/8000;
    boolean useFlatPriors = true;
    boolean DEBUGOUT = false;

    private HaplotypeIndelErrorModel model;
    private ArrayList<Integer> sitesVisited;

    protected DindelGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
        model = new HaplotypeIndelErrorModel(maxReadDeletionLength, insertionStartProbability,
                insertionEndProbability, alphaDeletionProbability, HAPLOTYPE_SIZE, false, DEBUGOUT);
        sitesVisited = new  ArrayList<Integer>();
    }

    public Allele getLikelihoods(RefMetaDataTracker tracker,
                                 ReferenceContext ref,
                                 Map<String, StratifiedAlignmentContext> contexts,
                                 StratifiedAlignmentContext.StratifiedContextType contextType,
                                 GenotypePriors priors,
                                 Map<String, BiallelicGenotypeLikelihoods> GLs,
                                 Allele alternateAlleleToUse) {  // TODO: use this instead of reading the 'indels' ROD from the tracker below                               

        if ( tracker == null )
            return null;

        VariantContext vc = null;

        for( final VariantContext vc_input : tracker.getVariantContexts(ref, "indels", null, ref.getLocus(), false, false) ) {
            if( vc_input != null && vc_input.isIndel() && ref.getLocus().getStart() == vc_input.getStart()) {
                vc = vc_input;
                break;
            }
        }
        // ignore places where we don't have a variant
        if ( vc == null )
            return null;


        if (!vc.isIndel())
            return null;

        boolean visitedBefore = false;
        synchronized (this) {

            if (sitesVisited.contains(new Integer(vc.getStart())) &&
                    contextType.equals(StratifiedAlignmentContext.StratifiedContextType.COMPLETE))
                 visitedBefore = true;
            else {
                sitesVisited.add(new Integer(vc.getStart()));
            }
        }

        if (visitedBefore)
                return null;

        // protect against having an indel too close to the edge of a contig
        if (vc.getStart() <= HAPLOTYPE_SIZE)
            return null;
        
        if ( !(priors instanceof DiploidIndelGenotypePriors) )
             throw new StingException("Only diploid-based Indel priors are supported in the DINDEL GL model");


        int eventLength = vc.getReference().getBaseString().length() - vc.getAlternateAllele(0).getBaseString().length();
        // assume only one alt allele for now
        if (eventLength<0)
            eventLength = - eventLength;

        int currentHaplotypeSize = HAPLOTYPE_SIZE;

        // int numSamples = getNSamples(contexts);
        List<Haplotype> haplotypesInVC = Haplotype.makeHaplotypeListFromVariantContextAlleles( vc, ref, currentHaplotypeSize);
        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.
        // initialize the GenotypeLikelihoods
        GLs.clear();

        double[][] haplotypeLikehoodMatrix;

        if (useFlatPriors) {
            priors = new DiploidIndelGenotypePriors();
        }
        else
            priors = new DiploidIndelGenotypePriors(indelHeterozygosity,eventLength,currentHaplotypeSize);

        //double[] priorLikelihoods = priors.getPriors();

        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            AlignmentContext context = sample.getValue().getContext(contextType);

            ReadBackedPileup pileup = null;
            if (context.hasExtendedEventPileup())
                pileup = context.getExtendedEventPileup();
            else if (context.hasBasePileup())
                pileup = context.getBasePileup();

            if (pileup != null ) {
                haplotypeLikehoodMatrix = model.computeReadHaplotypeLikelihoods( pileup, haplotypesInVC, vc, eventLength);


                double[] genotypeLikelihoods = HaplotypeIndelErrorModel.getHaplotypeLikelihoods( haplotypeLikehoodMatrix);

                GLs.put(sample.getKey(), new BiallelicGenotypeLikelihoods(sample.getKey(),
                        vc.getReference(),
                        vc.getAlternateAllele(0),
                        genotypeLikelihoods[0],
                        genotypeLikelihoods[1],
                        genotypeLikelihoods[2],
                        getFilteredDepth(pileup)));
                //System.out.format("%4.2f %4.2f %4.2f\n",genotypeLikelihoods[0],genotypeLikelihoods[1], genotypeLikelihoods[2]);
            }
        }

        return vc.getReference();
    }
}