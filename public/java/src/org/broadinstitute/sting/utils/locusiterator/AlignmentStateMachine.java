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

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.LinkedList;
import java.util.List;

/**
 * Steps a single read along its alignment to the genome
 *
 * The logical model for generating extended events is as follows: the "record state"
 * implements the traversal along the reference; thus stepForwardOnGenome() returns
 * on every and only on actual reference bases. This can be a (mis)match or a deletion
 * (in the latter case, we still return on every individual reference base the deletion spans).
 * In the extended events mode, the record state also remembers if there was an insertion, or
 * if the deletion just started *right before* the current reference base the record state is
 * pointing to upon the return from stepForwardOnGenome(). The next call to stepForwardOnGenome()
 * will clear that memory (as we remember only extended events immediately preceding
 * the current reference base).
 *
 * User: depristo
 * Date: 1/5/13
 * Time: 1:08 PM
 */
class AlignmentStateMachine {
    // TODO -- optimizations
    // TODO -- only keep 3 States, and recycle the prev state to become the next state

    /**
     * Our read
     */
    private final SAMRecord read;
    private final Cigar cigar;
    private final int nCigarElements;
    int cigarOffset = -1;

    AlignmentState prev = null, current = null, next = null;

    @Requires("read != null")
    // TODO -- should enforce contracts like the read is aligned, etc
    public AlignmentStateMachine(final SAMRecord read) {
        this.read = read;
        this.cigar = read.getCigar();
        this.nCigarElements = cigar.numCigarElements();
        this.prev = AlignmentState.makeLeftEdge(read);
    }

    public SAMRecord getRead() {
        return read;
    }

    public AlignmentState getPrev() {
        return prev;
    }

    public AlignmentState getCurrent() {
        return current;
    }

    public AlignmentState getNext() {
        return next;
    }

    @Deprecated
    public CigarElement peekForwardOnGenome() {
        return null;
    }

    @Deprecated
    public CigarElement peekBackwardOnGenome() {
        return null;
    }

    public CigarOperator stepForwardOnGenome() {
        if ( current == null ) {
            // start processing from the edge by updating current to be prev
            current = this.prev;
            current = nextAlignmentState();
        } else {
            // otherwise prev is current, and current is next
            prev = current;
            current = next;
        }

        // if the current pointer isn't the edge, update next
        if ( ! current.isEdge() )
            next = nextAlignmentState();
        else
            next = null;

        finalizeStates();

        // todo -- cleanup historical interface
        return current.isEdge() ? null : current.getCigarOperator();
    }

    private void finalizeStates() {
        // note the order of updates on the betweens.  Next has info, and then current does, so
        // the update order is next updates current, and current update prev

        if ( next != null ) {
            // next can be null because current is the edge
            assert ! current.isEdge();

            next.setPrev(current);

            // Next holds the info about what happened between
            // current and next, so we propagate it to current
            current.setBetweenNextPosition(next.getBetweenPrevPosition());
        }

        // TODO -- prev setting to current is not necessary (except in creating the left edge)
        prev.setNext(current);
        prev.setBetweenNextPosition(current.getBetweenPrevPosition());

        // current just needs to set prev and next
        current.setPrev(prev);
        current.setNext(next);

    }

    private AlignmentState nextAlignmentState() {
        int cigarElementCounter = getCurrent().getCigarElementCounter();
        CigarElement curElement = getCurrent().getCigarElement();
        int genomeOffset = getCurrent().getGenomeOffset();
        int readOffset = getCurrent().getReadOffset();

        // todo -- optimization: could keep null and allocate lazy since most of the time the between is empty
        final LinkedList<CigarElement> betweenCurrentAndNext = new LinkedList<CigarElement>();

        boolean done = false;
        while ( ! done ) {
            // we enter this method with readOffset = index of the last processed base on the read
            // (-1 if we did not process a single base yet); this can be last matching base,
            // or last base of an insertion
            if (curElement == null || ++cigarElementCounter > curElement.getLength()) {
                cigarOffset++;
                if (cigarOffset < nCigarElements) {
                    curElement = cigar.getCigarElement(cigarOffset);
                    cigarElementCounter = 0;
                    // next line: guards against cigar elements of length 0; when new cigar element is retrieved,
                    // we reenter in order to re-check cigarElementCounter against curElement's length
                } else {
                    if (curElement != null && curElement.getOperator() == CigarOperator.D)
                        throw new UserException.MalformedBAM(read, "read ends with deletion. Cigar: " + read.getCigarString() + ". Although the SAM spec technically permits such reads, this is often indicative of malformed files. If you are sure you want to use this file, re-run your analysis with the extra option: -rf BadCigar");
                    return AlignmentState.makeRightEdge(read, getCurrent(), betweenCurrentAndNext);
                }

                // in either case we continue the loop
                continue;
            }

            switch (curElement.getOperator()) {
                case H: // ignore hard clips
                case P: // ignore pads
                    cigarElementCounter = curElement.getLength();
                    betweenCurrentAndNext.add(curElement);
                    break;
                case I: // insertion w.r.t. the reference
                case S: // soft clip
                    cigarElementCounter = curElement.getLength();
                    readOffset += curElement.getLength();
                    betweenCurrentAndNext.add(curElement);
                    break;
                case D: // deletion w.r.t. the reference
                    if (readOffset < 0)             // we don't want reads starting with deletion, this is a malformed cigar string
                        throw new UserException.MalformedBAM(read, "read starts with deletion. Cigar: " + read.getCigarString() + ". Although the SAM spec technically permits such reads, this is often indicative of malformed files. If you are sure you want to use this file, re-run your analysis with the extra option: -rf BadCigar");
                    // should be the same as N case
                    genomeOffset++;
                    done = true;
                    break;
                case N: // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                    genomeOffset++;
                    done = true;
                    break;
                case M:
                case EQ:
                case X:
                    readOffset++;
                    genomeOffset++;
                    done = true;
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with cigar op: " + curElement.getOperator());
            }
        }

        return AlignmentState.makeInternalNode(read, readOffset, genomeOffset, curElement, cigarElementCounter, betweenCurrentAndNext);
    }
}
