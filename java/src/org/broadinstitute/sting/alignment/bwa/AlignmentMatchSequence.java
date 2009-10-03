package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.utils.StingException;

import java.util.Deque;
import java.util.ArrayDeque;

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
            throw new StingException("Unable to clone AlignmentMatchSequence.");
        }

        copy.entries = new ArrayDeque<AlignmentMatchSequenceEntry>();
        for( AlignmentMatchSequenceEntry entry: entries )
            copy.entries.add(entry.clone());

        return copy;
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
                throw new StingException("Unable to clone AlignmentMatchSequenceEntry.");
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

