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
import org.broadinstitute.sting.utils.pileup.*;

import java.util.*;

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
    private boolean isExtended = false; // tells whether this alignment context is an extended event context

    // todo -- why are you storing reads separately each time?  There's a ReadBackedPileup object that's supposed to handle this
//    private ArrayList<SAMRecord>[] reads = new ArrayList[StratifiedContextType.values().length];
//    private ArrayList<Integer>[] offsets = new ArrayList[StratifiedContextType.values().length];

    private ArrayList<PileupElement>[] pileupElems = new ArrayList[StratifiedContextType.values().length];
    //
    // accessors
    //
    public GenomeLoc getLocation() { return loc; }
//    public ArrayList<SAMRecord> getReads(StratifiedContextType type) { return reads[type.ordinal()]; }
//    public ArrayList<Integer> getOffsets(StratifiedContextType type) { return offsets[type.ordinal()]; }

    public ArrayList<PileupElement> getPileupElements(StratifiedContextType type) {
        return pileupElems[type.ordinal()];
    }

//    public ArrayList<ExtendedEventPileupElement> getExtendedPileupElements(StratifiedContextType type) {
//        if ( ! isExtended ) throw new StingException("Extended read backed pileups requested from StratifiedAlignmentContext that holds simple pileups");
//
//        return (ArrayList<ExtendedEventPileupElement>)(pileupElems[type.ordinal()]);
//    }

    public StratifiedAlignmentContext(GenomeLoc loc) {
        this(loc,false);
    }

    public StratifiedAlignmentContext(GenomeLoc loc, boolean isExtended) {
        this.loc = loc;
        this.isExtended = isExtended;
        for ( int i = 0; i < StratifiedContextType.values().length; i++) {
            if ( isExtended ) pileupElems[i] = new ArrayList<PileupElement>();
            else pileupElems[i] = new ArrayList<PileupElement>();
        }
    }

    public AlignmentContext getContext(StratifiedContextType type) {
        int index = type.ordinal();
        if ( contexts[index] == null ) {
            if ( isExtended ) {
                contexts[index] = new AlignmentContext(loc , new ReadBackedExtendedEventPileup(loc, (ArrayList<ExtendedEventPileupElement>)((ArrayList<? extends PileupElement>)getPileupElements(type))));
            } else {
                contexts[index] = new AlignmentContext(loc, new ReadBackedPileup(loc, getPileupElements(type)));
            }
        }
        return contexts[index];
    }

    public void add(SAMRecord read, int offset) {
        if ( isExtended ) throw new StingException("Can not add read/offset without event type specified to the context holding extended events");
        if ( read.getReadNegativeStrandFlag() ) {
            pileupElems[StratifiedContextType.REVERSE.ordinal()].add(new PileupElement(read,offset));
        } else {
            pileupElems[StratifiedContextType.FORWARD.ordinal()].add(new PileupElement(read,offset));
        }
        pileupElems[StratifiedContextType.COMPLETE.ordinal()].add(new PileupElement(read,offset));
     }

    public void add(PileupElement p) {
//        if ( isExtended ) throw new StingException("Can not add simple pileup element to the context holding extended events");
        SAMRecord read = p.getRead();
        if ( read.getReadNegativeStrandFlag() ) {
            pileupElems[StratifiedContextType.REVERSE.ordinal()].add(p);
        } else {
            pileupElems[StratifiedContextType.FORWARD.ordinal()].add(p);
        }
        pileupElems[StratifiedContextType.COMPLETE.ordinal()].add(p);
     }

    public void add(SAMRecord read, int offset, int length, byte [] bases) {
        if ( ! isExtended ) throw new StingException("Can not add read/offset with event type specified to the context holding simple events");
        if ( read.getReadNegativeStrandFlag() ) {
            pileupElems[StratifiedContextType.REVERSE.ordinal()].add(new ExtendedEventPileupElement(read,offset,length,bases));
        } else {
            pileupElems[StratifiedContextType.FORWARD.ordinal()].add(new ExtendedEventPileupElement(read,offset,length,bases));
        }
        pileupElems[StratifiedContextType.COMPLETE.ordinal()].add(new ExtendedEventPileupElement(read,offset,length,bases));
     }

//    public void add(ExtendedEventPileupElement p) {
//        if ( ! isExtended ) throw new StingException("Can not add extended pileup element to the context holding simple events");
//        SAMRecord read = p.getRead();
//        if ( read.getReadNegativeStrandFlag() ) {
//            pileupElems[StratifiedContextType.REVERSE.ordinal()].add(p);
//        } else {
//            pileupElems[StratifiedContextType.FORWARD.ordinal()].add(p);
//        }
//        pileupElems[StratifiedContextType.COMPLETE.ordinal()].add(p);
//     }

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
        GenomeLoc loc = pileup.getLocation();

        for (PileupElement p : pileup )
            addToContext(contexts, p, loc, assumedSingleSample, collapseToThisSample);

        return contexts;
    }

    /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per sample.
     *
     * @param pileup                the original pileup
     *
     * @return a Map of sample name to StratifiedAlignmentContext
     *
     **/
    public static Map<String, StratifiedAlignmentContext> splitContextBySample(ReadBackedExtendedEventPileup pileup) {
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
    public static Map<String, StratifiedAlignmentContext> splitContextBySample(ReadBackedExtendedEventPileup pileup, String assumedSingleSample, String collapseToThisSample) {

        HashMap<String, StratifiedAlignmentContext> contexts = new HashMap<String, StratifiedAlignmentContext>();
        GenomeLoc loc = pileup.getLocation();

        for (PileupElement p : pileup )
            addToContext(contexts, p, loc, assumedSingleSample, collapseToThisSample,true);

        return contexts;
    }

    private static void addToContext(HashMap<String, StratifiedAlignmentContext> contexts, PileupElement p, GenomeLoc loc, String assumedSingleSample, String collapseToThisSample) {
        addToContext(contexts, p, loc, assumedSingleSample, collapseToThisSample, false);
    }

    private static void addToContext(HashMap<String, StratifiedAlignmentContext> contexts, PileupElement p, GenomeLoc loc, String assumedSingleSample, String collapseToThisSample, boolean isExtended) {

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
            myContext = new StratifiedAlignmentContext(loc,isExtended);
            contexts.put(sample, myContext);
        }

        // add the read to this sample's context
        // note that bad bases are added to the context (for DoC calculations later)
        myContext.add(p);
    }

     /**
     * Splits the given AlignmentContext into a StratifiedAlignmentContext per read group.
     *
     * @param pileup                the original pileup
     * @return a Map of sample name to StratifiedAlignmentContext
     * TODO - support for collapsing or assuming read groups if they are missing
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

            myContext.add(p);
        }

        return contexts;
    }

    public static AlignmentContext joinContexts(Collection<StratifiedAlignmentContext> contexts, StratifiedContextType type) {
        ArrayList<PileupElement> pe = new ArrayList();

        if ( contexts.size() == 0  )
            throw new StingException("BUG: joinContexts requires at least one context to join");


        Iterator<StratifiedAlignmentContext> it = contexts.iterator();
        StratifiedAlignmentContext context = it.next();
        boolean isExtended = context.isExtended;
        GenomeLoc loc = context.getLocation();
        pe.addAll(context.getPileupElements(type));

        while ( it.hasNext()) {
            context = it.next();
            if ( ! loc.equals( context.getLocation() ) )
                    throw new StingException("Illegal attempt to join contexts from different genomic locations");
            if ( context.isExtended != isExtended )
                throw new StingException("Illegal attempt to join simple and extended contexts");
            pe.addAll(context.getPileupElements(type));
        }

        // dirty trick below. generics do not allow to cast pe (ArrayList<PileupElement>) directly to ArrayList<ExtendedEventPileupElement>,
        // so we first cast to "? extends" wildcard, then to what we actually need.
        if ( isExtended ) return new AlignmentContext(loc, new ReadBackedExtendedEventPileup(loc, (ArrayList< ExtendedEventPileupElement>)((ArrayList<? extends PileupElement>)pe)) );
        else return new AlignmentContext(loc, new ReadBackedPileup(loc,pe));
    }
}