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

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Useful class for storing different AlignmentContexts
 * User: ebanks
 */
public class StratifiedAlignmentContext {

    public enum StratifiedContextType { OVERALL, FORWARD, REVERSE }

    private AlignmentContext overall = null;
    private AlignmentContext forward = null;
    private AlignmentContext reverse = null;
    private GenomeLoc loc;

    private ArrayList<SAMRecord> allReads = new ArrayList<SAMRecord>();
    private ArrayList<SAMRecord> forwardReads = new ArrayList<SAMRecord>();
    private ArrayList<SAMRecord> reverseReads = new ArrayList<SAMRecord>();

    private ArrayList<Integer> allOffsets = new ArrayList<Integer>();
    private ArrayList<Integer> forwardOffsets = new ArrayList<Integer>();
    private ArrayList<Integer> reverseOffsets = new ArrayList<Integer>();


    public StratifiedAlignmentContext(GenomeLoc loc) {
        this.loc = loc;
    }

    public AlignmentContext getContext(StratifiedContextType context) {
        switch ( context ) {
            case OVERALL: return getOverallContext();
            case FORWARD: return getForwardContext();
            case REVERSE: return getReverseContext();
        }
        return null;
    }

    private AlignmentContext getOverallContext() {
        if ( overall == null )
            overall = new AlignmentContext(loc, new ReadBackedPileup(loc, allReads, allOffsets));
        return overall;
    }

    private AlignmentContext getForwardContext() {
        if ( forward == null )
            forward = new AlignmentContext(loc, new ReadBackedPileup(loc, forwardReads, forwardOffsets));
        return forward;
    }

    private AlignmentContext getReverseContext() {
        if ( reverse == null )
            reverse = new AlignmentContext(loc, new ReadBackedPileup(loc, reverseReads, reverseOffsets));
        return reverse;
    }

    public void add(SAMRecord read, int offset) {
        if ( read.getReadNegativeStrandFlag() ) {
            reverseReads.add(read);
            reverseOffsets.add(offset);
        } else {
            forwardReads.add(read);
            forwardOffsets.add(offset);
        }
        allReads.add(read);
        allOffsets.add(offset);
     }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample.
     *
     * @param context               the original AlignmentContext
     * @param assumedSingleSample   if not null, any read without a readgroup will be given this sample name
     * @param collapseToThisSample  if not null, all reads will be assigned this read group regardless of their actual read group
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static Map<String, StratifiedAlignmentContext> splitContextBySample(AlignmentContext context, String assumedSingleSample, String collapseToThisSample) {

        HashMap<String, StratifiedAlignmentContext> contexts = new HashMap<String, StratifiedAlignmentContext>();

        ReadBackedPileup pileup = context.getPileup();
        for (PileupElement p : pileup ) {

            // get the read
            SAMRecord read = p.getRead();

            // find the sample
            String sample;
            if ( collapseToThisSample != null ) {
                sample = collapseToThisSample;
            } else {
                SAMReadGroupRecord readGroup = read.getReadGroup();
                if ( readGroup == null ) {
                    if ( assumedSingleSample == null )
                        throw new StingException("Missing read group for read " + read.getReadName());
                    sample = assumedSingleSample;
                } else {
                    sample = readGroup.getSample();
                }
            }

            // create a new context object if this is the first time we're seeing a read for this sample
            StratifiedAlignmentContext myContext = contexts.get(sample);
            if ( myContext == null ) {
                myContext = new StratifiedAlignmentContext(context.getLocation());
                contexts.put(sample, myContext);
            }

            // add the read to this sample's context
            // note that bad bases are added to the context (for DoC calculations later)
            myContext.add(read, p.getOffset());
        }

        return contexts;
    }
}