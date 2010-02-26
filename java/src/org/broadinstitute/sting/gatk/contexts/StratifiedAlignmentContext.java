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
 * Modified: chartl (split by read group)
 */
public class StratifiedAlignmentContext {

    // Definitions:
    //   COMPLETE = full alignment context
    //   FORWARD  = reads on forward strand
    //   REVERSE  = reads on forward strand
    //
    public enum StratifiedContextType { COMPLETE, FORWARD, REVERSE }

    private GenomeLoc loc;
    private AlignmentContext[] contexts = new AlignmentContext[StratifiedContextType.values().length];
    private ArrayList<SAMRecord>[] reads = new ArrayList[StratifiedContextType.values().length];
    private ArrayList<Integer>[] offsets = new ArrayList[StratifiedContextType.values().length];


    public StratifiedAlignmentContext(GenomeLoc loc) {
        this.loc = loc;
        for ( int i = 0; i < StratifiedContextType.values().length; i++) {
            reads[i] = new ArrayList<SAMRecord>();
            offsets[i] = new ArrayList<Integer>();
        }
    }

    public AlignmentContext getContext(StratifiedContextType context) {
        int index = context.ordinal();
        if ( contexts[index] == null )
            contexts[index] = new AlignmentContext(loc, new ReadBackedPileup(loc, reads[index], offsets[index]));
        return contexts[index];
    }

    public void add(SAMRecord read, int offset) {
        if ( read.getReadNegativeStrandFlag() ) {
            reads[StratifiedContextType.REVERSE.ordinal()].add(read);
            offsets[StratifiedContextType.REVERSE.ordinal()].add(offset);
        } else {
            reads[StratifiedContextType.FORWARD.ordinal()].add(read);
            offsets[StratifiedContextType.FORWARD.ordinal()].add(offset);
        }
        reads[StratifiedContextType.COMPLETE.ordinal()].add(read);
        offsets[StratifiedContextType.COMPLETE.ordinal()].add(offset);
     }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample.
     *
     * @param pileup                the original pileup
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static Map<String, StratifiedAlignmentContext> splitContextBySample(ReadBackedPileup pileup) {
        return splitContextBySample(pileup, null, null);
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample.
     *
     * @param pileup                the original pileup
     * @param assumedSingleSample   if not null, any read without a readgroup will be given this sample name
     * @param collapseToThisSample  if not null, all reads will be assigned this read group regardless of their actual read group
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static Map<String, StratifiedAlignmentContext> splitContextBySample(ReadBackedPileup pileup, String assumedSingleSample, String collapseToThisSample) {

        HashMap<String, StratifiedAlignmentContext> contexts = new HashMap<String, StratifiedAlignmentContext>();

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
                myContext = new StratifiedAlignmentContext(pileup.getLocation());
                contexts.put(sample, myContext);
            }

            // add the read to this sample's context
            // note that bad bases are added to the context (for DoC calculations later)
            myContext.add(read, p.getOffset());
        }

        return contexts;
    }

     /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per read group.
     *
     * @param pileup                the original pileup
     * @return a Map of sample name to StratifiedAlignmentContext
     * @todo - support for collapsing or assuming read groups if they are missing
     *
     **/
    public static Map<String,StratifiedAlignmentContext> splitContextByReadGroup(ReadBackedPileup pileup) {
        HashMap<String,StratifiedAlignmentContext> contexts = new HashMap<String,StratifiedAlignmentContext>();

        for ( PileupElement p : pileup ) {
            SAMRecord read = p.getRead();

            SAMReadGroupRecord readGroup = read.getReadGroup();
            if ( readGroup == null ) {
                throw new StingException("Missing read group for read " + read.getReadName());
            }

            String group = readGroup.getReadGroupId();

            StratifiedAlignmentContext myContext = contexts.get(group);

            if ( myContext == null ) {
                myContext = new StratifiedAlignmentContext(pileup.getLocation());
                contexts.put(group,myContext);
            }

            myContext.add(read,p.getOffset());
        }

        return contexts;
    }
}