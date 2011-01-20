/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.contexts;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.*;

import java.util.*;

/**
 * Useful class for storing different AlignmentContexts
 * User: ebanks
 * Modified: chartl (split by read group)
 */
public class StratifiedAlignmentContext<RBP extends ReadBackedPileup> implements HasGenomeLocation {

    // Definitions:
    //   COMPLETE = full alignment context
    //   FORWARD  = reads on forward strand
    //   REVERSE  = reads on forward strand
    //
    public enum StratifiedContextType { COMPLETE, FORWARD, REVERSE }

    private GenomeLoc loc;
    private RBP basePileup = null;

    //
    // accessors
    //
    public GenomeLoc getLocation() { return loc; }

    public StratifiedAlignmentContext(GenomeLoc loc) {
        this(loc,null);
    }

    public StratifiedAlignmentContext(GenomeLoc loc, RBP pileup) {
        this.loc = loc;
        this.basePileup = pileup;
    }

    public AlignmentContext getContext(StratifiedContextType type) {
        switch(type) {
            case COMPLETE:
                return new AlignmentContext(loc,basePileup);
            case FORWARD:
                return new AlignmentContext(loc,basePileup.getPositiveStrandPileup());
            case REVERSE:
                return new AlignmentContext(loc,basePileup.getNegativeStrandPileup());
            default:
                throw new ReviewedStingException("Unable to get alignment context for type = " + type);
        }
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample.
     *
     * @param pileup                the original pileup
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static <RBP extends ReadBackedPileup,PE extends PileupElement> Map<Sample, StratifiedAlignmentContext> splitContextBySample(RBP pileup) {
        return splitContextBySample(pileup, null);
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample.
     *
     * @param pileup                the original pileup
     * @param assumedSingleSample   if not null, any read without a readgroup will be given this sample name
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static <RBP extends ReadBackedPileup> Map<Sample, StratifiedAlignmentContext> splitContextBySample(RBP pileup, Sample assumedSingleSample) {

        GenomeLoc loc = pileup.getLocation();
        HashMap<Sample, StratifiedAlignmentContext> contexts = new HashMap<Sample, StratifiedAlignmentContext>();

        for(Sample sample: pileup.getSamples()) {
            RBP pileupBySample = (RBP)pileup.getPileupForSample(sample);

            // Don't add empty pileups to the split context.
            if(pileupBySample.size() == 0)
                continue;

            if(sample != null)
                contexts.put(sample,new StratifiedAlignmentContext<RBP>(loc,pileupBySample));
            else {
                if(assumedSingleSample == null) {
                    throw new UserException.ReadMissingReadGroup(pileupBySample.iterator().next().getRead());
                }
                contexts.put(assumedSingleSample,new StratifiedAlignmentContext<RBP>(loc,pileupBySample));
            }
        }

        return contexts;
    }



    public static <RBP extends ReadBackedPileup,PE extends PileupElement> Map<String, StratifiedAlignmentContext> splitContextBySampleName(RBP pileup) {
        return splitContextBySampleName(pileup, null);
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample, but referencd by sample name instead
     * of sample object.
     *
     * @param pileup                the original pileup
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static <RBP extends ReadBackedPileup> Map<String, StratifiedAlignmentContext> splitContextBySampleName(RBP pileup, String assumedSingleSample) {

        GenomeLoc loc = pileup.getLocation();
        HashMap<String, StratifiedAlignmentContext> contexts = new HashMap<String, StratifiedAlignmentContext>();

        for(String sample: pileup.getSampleNames()) {
            RBP pileupBySample = (RBP)pileup.getPileupForSampleName(sample);

            // Don't add empty pileups to the split context.
            if(pileupBySample.size() == 0)
                continue;

            if(sample != null)
                contexts.put(sample,new StratifiedAlignmentContext<RBP>(loc,pileupBySample));
            else {
                if(assumedSingleSample == null) {
                    throw new UserException.ReadMissingReadGroup(pileupBySample.iterator().next().getRead());
                }
                contexts.put(assumedSingleSample,new StratifiedAlignmentContext<RBP>(loc,pileupBySample));
            }
        }

        return contexts;
    }


    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per read group.
     *
     * @param pileup                the original pileup
     * @return a Map of sample name to StratifiedAlignmentContext
     * TODO - support for collapsing or assuming read groups if they are missing
     *
     **/
    public static <RBP extends ReadBackedPileup> Map<String,StratifiedAlignmentContext<RBP>> splitContextByReadGroup(RBP pileup) {
        HashMap<String,StratifiedAlignmentContext<RBP>> contexts = new HashMap<String,StratifiedAlignmentContext<RBP>>();
        for(String readGroupId: pileup.getReadGroups())
            contexts.put(readGroupId,new StratifiedAlignmentContext<RBP>(pileup.getLocation(),(RBP)pileup.getPileupForReadGroup(readGroupId)));
        return contexts;
    }

    public static AlignmentContext joinContexts(Collection<StratifiedAlignmentContext> contexts) {

        // validation
        GenomeLoc loc = contexts.iterator().next().getLocation();
        boolean isExtended = contexts.iterator().next().basePileup instanceof ReadBackedExtendedEventPileup;
        for(StratifiedAlignmentContext context: contexts) {
            if(!loc.equals(context.getLocation()))
                throw new ReviewedStingException("Illegal attempt to join contexts from different genomic locations");
            if(isExtended != (context.basePileup instanceof ReadBackedExtendedEventPileup))
                throw new ReviewedStingException("Illegal attempt to join simple and extended contexts");
        }

        AlignmentContext jointContext;
        if(isExtended) {
            List<ExtendedEventPileupElement> pe = new ArrayList<ExtendedEventPileupElement>();
            for(StratifiedAlignmentContext context: contexts) {
                for(PileupElement pileupElement: context.basePileup)
                    pe.add((ExtendedEventPileupElement)pileupElement);
            }
            jointContext = new AlignmentContext(loc, new ReadBackedExtendedEventPileupImpl(loc,pe));
        }
        else {
            List<PileupElement> pe = new ArrayList<PileupElement>();
            for(StratifiedAlignmentContext context: contexts) {
                for(PileupElement pileupElement: context.basePileup)
                    pe.add(pileupElement);
            }
            jointContext = new AlignmentContext(loc, new ReadBackedPileupImpl(loc,pe));
        }

        return jointContext;
    }
}