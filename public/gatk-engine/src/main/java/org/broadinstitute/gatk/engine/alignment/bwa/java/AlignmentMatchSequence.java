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

package org.broadinstitute.gatk.engine.alignment.bwa.java;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;

/**
 * Represents a sequence of matches.
 *
 * @author mhanna
 * @version 0.1
 */
public class AlignmentMatchSequence implements Cloneable {
    /**
     * Stores the particular match entries in the order they occur.
     */
    private Deque<AlignmentMatchSequenceEntry> entries = new ArrayDeque<AlignmentMatchSequenceEntry>();

    /**
     * Clone the given match sequence.
     * @return A deep copy of the current match sequence.
     */
    public AlignmentMatchSequence clone() {
        AlignmentMatchSequence copy = null;
        try {
            copy = (AlignmentMatchSequence)super.clone(); 
        }
        catch( CloneNotSupportedException ex ) {
            throw new ReviewedGATKException("Unable to clone AlignmentMatchSequence.");
        }

        copy.entries = new ArrayDeque<AlignmentMatchSequenceEntry>();
        for( AlignmentMatchSequenceEntry entry: entries )
            copy.entries.add(entry.clone());

        return copy;
    }

    public Cigar convertToCigar(boolean negativeStrand) {
        Cigar cigar = new Cigar();
        Iterator<AlignmentMatchSequenceEntry> iterator = negativeStrand ? entries.descendingIterator() : entries.iterator();
        while( iterator.hasNext() ) {
            AlignmentMatchSequenceEntry entry = iterator.next();
            CigarOperator operator;
            switch( entry.getAlignmentState() ) {
                case MATCH_MISMATCH: operator = CigarOperator.MATCH_OR_MISMATCH; break;
                case INSERTION: operator = CigarOperator.INSERTION; break;
                case DELETION: operator = CigarOperator.DELETION; break;
                default: throw new ReviewedGATKException("convertToCigar: cannot process state: " + entry.getAlignmentState());
            }
            cigar.add( new CigarElement(entry.count,operator) );
        }
        return cigar;
    }

    /**
     * All a new alignment of the given state.
     * @param state State to add to the sequence.
     */
    public void addNext( AlignmentState state ) {
        AlignmentMatchSequenceEntry last = entries.peekLast();
        // If the last entry is the same as this one, increment it.  Otherwise, add a new entry.
        if( last != null && last.alignmentState == state )
            last.increment();
        else
            entries.add(new AlignmentMatchSequenceEntry(state));
    }

    /**
     * Gets the current state of this alignment (what's the state of the last base?)
     * @return State of the most recently aligned base.
     */
    public AlignmentState getCurrentState() {
        if( entries.size() == 0 )
            return AlignmentState.MATCH_MISMATCH;        
        return entries.peekLast().getAlignmentState();
    }

    /**
     * How many bases in the read match the given state.
     * @param state State to test.
     * @return number of bases which match that state.
     */
    public int getNumberOfBasesMatchingState(AlignmentState state) {
        int matches = 0;
        for( AlignmentMatchSequenceEntry entry: entries ) {
            if( entry.getAlignmentState() == state )
                matches += entry.count;
        }
        return matches;
    }

    /**
     * Stores an individual match sequence entry.
     */
    private class AlignmentMatchSequenceEntry implements Cloneable {
        /**
         * The state of the alignment throughout a given point in the sequence.
         */
        private final AlignmentState alignmentState;

        /**
         * The number of bases having this particular state.
         */
        private int count;

        /**
         * Create a new sequence entry with the given state.
         * @param alignmentState The state that this sequence should contain.
         */
        AlignmentMatchSequenceEntry( AlignmentState alignmentState ) {
            this.alignmentState = alignmentState;
            this.count = 1;
        }

        /**
         * Clone the given match sequence entry.
         * @return A deep copy of the current match sequence entry.
         */
        public AlignmentMatchSequenceEntry clone() {
            try {
                return (AlignmentMatchSequenceEntry)super.clone(); 
            }
            catch( CloneNotSupportedException ex ) {
                throw new ReviewedGATKException("Unable to clone AlignmentMatchSequenceEntry.");
            }
        }

        /**
         * Retrieves the current state of the alignment.
         * @return The state of the current sequence.
         */
        AlignmentState getAlignmentState() {
            return alignmentState;
        }

        /**
         * Increment the count of alignments having this particular state.
         */
        void increment() {
            count++;
        }
    }
}

