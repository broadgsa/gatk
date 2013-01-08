/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.locusiterator;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.LinkedList;
import java.util.List;

public final class AlignmentState {
    /**
     * Our read
     */
    private final SAMRecord read;

    /**
     * how far are we offset from the start of the read bases?
     */
    private final int readOffset;

    /**
     * how far are we offset from the alignment start on the genome?
     */
    private final int genomeOffset;

    /**
     * Our cigar element
     */
    private final CigarElement cigarElement;

    /**
     * how far are we into our cigarElement?
     */
    private final int cigarElementCounter;

    private LinkedList<CigarElement> betweenPrevPosition = null, betweenNextPosition = null;
    private AlignmentState prev = null, next = null;

    public static AlignmentState makeInternalNode(final SAMRecord read, int readOffset,
                                                  int genomeOffset, CigarElement cigarElement,
                                                  int cigarElementCounter, final LinkedList<CigarElement> betweenPrevAndThis) {
        final AlignmentState state = new AlignmentState(read, readOffset, genomeOffset, cigarElement, cigarElementCounter);
        state.setBetweenPrevPosition(betweenPrevAndThis);
        return state;
    }

    public static AlignmentState makeLeftEdge(final SAMRecord read) {
        return new AlignmentState(read, -1, 1, null, -1);
    }

    public static AlignmentState makeRightEdge(final SAMRecord read, final AlignmentState current, final LinkedList<CigarElement> betweenCurrentAndThis) {
        final AlignmentState state = new AlignmentState(read, -1, 1, null, -1);
        state.setPrev(current);
        state.setBetweenPrevPosition(betweenCurrentAndThis);
        return state;
    }

    protected AlignmentState(SAMRecord read, int readOffset, int genomeOffset, CigarElement cigarElement, int cigarElementCounter) {
        this.read = read;
        this.readOffset = readOffset;
        this.genomeOffset = genomeOffset;
        this.cigarElement = cigarElement;
        this.cigarElementCounter = cigarElementCounter;
    }

    /**
     * Is this an edge state?  I.e., one that is before or after the current read?
     * @return true if this state is an edge state, false otherwise
     */
    public boolean isEdge() {
        return readOffset == -1;
    }

    public SAMRecord getRead() {
        return read;
    }

    /**
     * What is our current offset in the read's bases that aligns us with the reference genome?
     *
     * @return the current read offset position
     */
    public int getReadOffset() {
        return readOffset;
    }

    /**
     * What is the current offset w.r.t. the alignment state that aligns us to the readOffset?
     *
     * @return the current offset
     */
    public int getGenomeOffset() {
        return genomeOffset;
    }

    public int getGenomePosition() {
        return read.getAlignmentStart() + getGenomeOffset();
    }

    public GenomeLoc getLocation(final GenomeLocParser genomeLocParser) {
        return genomeLocParser.createGenomeLoc(read.getReferenceName(), getGenomePosition());
    }

    public AlignmentState getPrev() {
        return prev;
    }

    public AlignmentState getNext() {
        return next;
    }

    public boolean hasPrev() { return prev != null; }
    public boolean hasNext() { return next != null; }
    public boolean prevIsEdge() { return hasPrev() && getPrev().isEdge(); }
    public boolean nextIsEdge() { return hasNext() && getNext().isEdge(); }

    public CigarElement getCigarElement() {
        return cigarElement;
    }

    /**
     *
     * @return null if this is an edge state
     */
    public CigarOperator getCigarOperator() {
        return cigarElement == null ? null : cigarElement.getOperator();
    }

    public String toString() {
        return String.format("%s ro=%d go=%d cec=%d %s", read.getReadName(), readOffset, genomeOffset, cigarElementCounter, cigarElement);
    }

    public int getCigarElementCounter() {
        return cigarElementCounter;
    }

    // -----------------------------------------------------------------------------------------------
    // Code for setting up prev / next states
    //
    // TODO -- should these functions all be protected?
    //
    // -----------------------------------------------------------------------------------------------

    public void setBetweenPrevPosition(LinkedList<CigarElement> betweenPrevPosition) {
        this.betweenPrevPosition = betweenPrevPosition;
    }

    public void setBetweenNextPosition(LinkedList<CigarElement> betweenNextPosition) {
        this.betweenNextPosition = betweenNextPosition;
    }

    public LinkedList<CigarElement> getBetweenPrevPosition() {
        return betweenPrevPosition;
    }

    public LinkedList<CigarElement> getBetweenNextPosition() {
        return betweenNextPosition;
    }

    public void setPrev(AlignmentState prev) {
        this.prev = prev;
    }

    public void setNext(AlignmentState next) {
        this.next = next;
    }

    // -----------------------------------------------------------------------------------------------
    // Code for computing presence / absence of states in the prev / current / next
    // -----------------------------------------------------------------------------------------------

    public boolean isAfterDeletion() { return testOperator(getPrev(), CigarOperator.D); }
    public boolean isBeforeDeletion() { return testOperator(getNext(), CigarOperator.D); }
    public boolean isAfterInsertion() { return isAfter(getBetweenPrevPosition(), CigarOperator.I); }
    public boolean isBeforeInsertion() { return isBefore(getBetweenNextPosition(), CigarOperator.I); }

    public boolean isAfterSoftClip() { return isAfter(getBetweenPrevPosition(), CigarOperator.S); }
    public boolean isBeforeSoftClip() { return isBefore(getBetweenNextPosition(), CigarOperator.S); }
    public boolean isNextToSoftClip() { return isAfterSoftClip() || isBeforeSoftClip(); }

    private boolean testOperator(final AlignmentState state, final CigarOperator op) {
        return state != null && state.getCigarOperator() == op;
    }

    private boolean isAfter(final LinkedList<CigarElement> elements, final CigarOperator op) {
        return ! elements.isEmpty() && elements.peekLast().getOperator() == op;
    }

    private boolean isBefore(final List<CigarElement> elements, final CigarOperator op) {
        return ! elements.isEmpty() && elements.get(0).getOperator() == op;
    }
}
