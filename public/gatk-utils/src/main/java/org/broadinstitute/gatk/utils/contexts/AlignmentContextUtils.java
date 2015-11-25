/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.contexts;

import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.*;

import java.util.*;

/**
 * Useful utilities for storing different AlignmentContexts
 * User: ebanks
 */
public class AlignmentContextUtils {

    // Definitions:
    //   COMPLETE = full alignment context
    //   FORWARD  = reads on forward strand
    //   REVERSE  = reads on forward strand
    //
    public enum ReadOrientation { COMPLETE, FORWARD, REVERSE }

    private AlignmentContextUtils() {
        // cannot be instantiated
    }

    /**
     * Returns a potentially derived subcontext containing only forward, reverse, or in fact all reads
     * in alignment context context.
     *
     * @param context
     * @param type
     * @return
     */
    public static AlignmentContext stratify(AlignmentContext context, ReadOrientation type) {
        switch(type) {
            case COMPLETE:
                return context;
            case FORWARD:
                return new AlignmentContext(context.getLocation(),context.getPileup().getPositiveStrandPileup());
            case REVERSE:
                return new AlignmentContext(context.getLocation(),context.getPileup().getNegativeStrandPileup());
            default:
                throw new ReviewedGATKException("Unable to get alignment context for type = " + type);
        }
    }

    public static Map<String, AlignmentContext> splitContextBySampleName(AlignmentContext context) {
        return splitContextBySampleName(context, null);
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample, but referencd by sample name instead
     * of sample object.
     *
     * @param context                the original pileup
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static Map<String, AlignmentContext> splitContextBySampleName(AlignmentContext context, String assumedSingleSample) {
        GenomeLoc loc = context.getLocation();
        HashMap<String, AlignmentContext> contexts = new HashMap<String, AlignmentContext>();

        for(String sample: context.getPileup().getSamples()) {
            ReadBackedPileup pileupBySample = context.getPileup().getPileupForSample(sample);

            // Don't add empty pileups to the split context.
            if(pileupBySample.getNumberOfElements() == 0)
                continue;

            if(sample != null)
                contexts.put(sample, new AlignmentContext(loc, pileupBySample));
            else {
                if(assumedSingleSample == null) {
                    throw new UserException.ReadMissingReadGroup(pileupBySample.iterator().next().getRead());
                }
                contexts.put(assumedSingleSample,new AlignmentContext(loc, pileupBySample));
            }
        }

        return contexts;
    }

    /**
     * Splits the AlignmentContext into one context per read group
     *
     * @param context the original pileup
     * @return a Map of ReadGroup to AlignmentContext, or an empty map if context has no base pileup
     *
     **/
    public static Map<SAMReadGroupRecord, AlignmentContext> splitContextByReadGroup(AlignmentContext context, Collection<SAMReadGroupRecord> readGroups) {
        HashMap<SAMReadGroupRecord, AlignmentContext> contexts = new HashMap<SAMReadGroupRecord, AlignmentContext>();

        for (SAMReadGroupRecord rg : readGroups) {
            ReadBackedPileup rgPileup = context.getBasePileup().getPileupForReadGroup(rg.getReadGroupId());
            if ( rgPileup != null ) // there we some reads for RG
                contexts.put(rg, new AlignmentContext(context.getLocation(), rgPileup));
        }

        return contexts;
    }

    public static Map<String, AlignmentContext> splitContextBySampleName(ReadBackedPileup pileup) {
        return splitContextBySampleName(new AlignmentContext(pileup.getLocation(), pileup));
    }


    public static AlignmentContext joinContexts(Collection<AlignmentContext> contexts) {
        // validation
        GenomeLoc loc = contexts.iterator().next().getLocation();
        for(AlignmentContext context: contexts) {
            if(!loc.equals(context.getLocation()))
                throw new ReviewedGATKException("Illegal attempt to join contexts from different genomic locations");
        }

        List<PileupElement> pe = new ArrayList<PileupElement>();
        for(AlignmentContext context: contexts) {
            for(PileupElement pileupElement: context.basePileup)
                pe.add(pileupElement);
        }
        return new AlignmentContext(loc, new ReadBackedPileupImpl(loc,pe));
    }
}