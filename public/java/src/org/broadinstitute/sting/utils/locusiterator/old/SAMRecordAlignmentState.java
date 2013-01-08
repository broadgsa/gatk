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

package org.broadinstitute.sting.utils.locusiterator.old;

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

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
public class SAMRecordAlignmentState {
    // TODO -- one idea to clean up this functionality:
    // TODO --
    // TODO -- split functionality here into an alignment state machine and an
    // TODO -- alignment state.  The alignment state simply carries with it the
    // TODO -- state of the alignment (the current cigar op, the genome offset,
    // TODO -- the read offset, etc.  The AlignmentStateMachine produces these
    // TODO -- states, and has operations such stepForwardOnGenome, getLastState(),
    // TODO -- getCurrentState(), getNextState();

    /**
     * Our read
     */
    private final SAMRecord read;
    private final Cigar cigar;
    private final int nCigarElements;

    /**
     * how far are we offset from the start of the read bases?
     */
    int readOffset = -1;

    /**
     * how far are we offset from the alignment start on the genome?
     */
    int genomeOffset = -1;

    int cigarOffset = -1;
    CigarElement curElement = null;

    /**
     * how far are we into a single cigarElement?
     */
    int cigarElementCounter = -1;

    @Requires("read != null")
    // TODO -- should enforce contracts like the read is aligned, etc
    public SAMRecordAlignmentState(final SAMRecord read) {
        this.read = read;
        this.cigar = read.getCigar();
        this.nCigarElements = cigar.numCigarElements();
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

    public GenomeLoc getLocation(GenomeLocParser genomeLocParser) {
        return genomeLocParser.createGenomeLoc(read.getReferenceName(), getGenomePosition());
    }

    public CigarOperator getCurrentCigarOperator() {
        return curElement.getOperator();
    }

    public String toString() {
        return String.format("%s ro=%d go=%d co=%d cec=%d %s", read.getReadName(), readOffset, genomeOffset, cigarOffset, cigarElementCounter, curElement);
    }

    public CigarElement peekForwardOnGenome() {
        return ( cigarElementCounter + 1 > curElement.getLength() && cigarOffset + 1 < nCigarElements ? cigar.getCigarElement(cigarOffset + 1) : curElement );
    }

    public CigarElement peekBackwardOnGenome() {
        return ( cigarElementCounter - 1 == 0 && cigarOffset - 1 > 0 ? cigar.getCigarElement(cigarOffset - 1) : curElement );
    }

    public CigarOperator stepForwardOnGenome() {
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
                return stepForwardOnGenome();
            } else {
                if (curElement != null && curElement.getOperator() == CigarOperator.D)
                    throw new UserException.MalformedBAM(read, "read ends with deletion. Cigar: " + read.getCigarString() + ". Although the SAM spec technically permits such reads, this is often indicative of malformed files. If you are sure you want to use this file, re-run your analysis with the extra option: -rf BadCigar");

                // Reads that contain indels model the genomeOffset as the following base in the reference.  Because
                // we fall into this else block only when indels end the read, increment genomeOffset  such that the
                // current offset of this read is the next ref base after the end of the indel.  This position will
                // model a point on the reference somewhere after the end of the read.
                genomeOffset++; // extended events need that. Logically, it's legal to advance the genomic offset here:
                // we do step forward on the ref, and by returning null we also indicate that we are past the read end.

                return null;
            }
        }

        boolean done = false;
        switch (curElement.getOperator()) {
            case H: // ignore hard clips
            case P: // ignore pads
                cigarElementCounter = curElement.getLength();
                break;
            case I: // insertion w.r.t. the reference
            case S: // soft clip
                cigarElementCounter = curElement.getLength();
                readOffset += curElement.getLength();
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

        return done ? curElement.getOperator() : stepForwardOnGenome();
    }
}
