///*
// * Copyright (c) 2012 The Broad Institute
// *
// * Permission is hereby granted, free of charge, to any person
// * obtaining a copy of this software and associated documentation
// * files (the "Software"), to deal in the Software without
// * restriction, including without limitation the rights to use,
// * copy, modify, merge, publish, distribute, sublicense, and/or sell
// * copies of the Software, and to permit persons to whom the
// * Software is furnished to do so, subject to the following
// * conditions:
// *
// * The above copyright notice and this permission notice shall be
// * included in all copies or substantial portions of the Software.
// *
// * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
// * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// */
//
//package org.broadinstitute.sting.utils.locusiterator;
//
//import com.google.java.contract.Invariant;
//import net.sf.samtools.CigarElement;
//import net.sf.samtools.CigarOperator;
//import net.sf.samtools.SAMRecord;
//import org.broadinstitute.sting.utils.GenomeLoc;
//import org.broadinstitute.sting.utils.GenomeLocParser;
//
//import java.util.LinkedList;
//import java.util.List;
//
//@Invariant({
//        "read != null",
//        "readOffset >= -1",
////        "readOffset < read.getReadLength()",
//        "genomeOffset >= -1",
//        // if read offset == -1 then genome offset and cigarElementCounter must also be -1
//        //TODO "readOffset != -1 || (genomeOffset == -1 && cigarElementCounter == -1)",
//        "cigarElementCounter >= -1",
//        // either there's no cigar element of the counter < its length
//        //TODO "cigarElement == null || cigarElementCounter < cigarElement.getLength()"
//})
//public final class AlignmentState {
//    /**
//     * Our read
//     */
//    private final SAMRecord read;
//
//    private LinkedList<CigarElement> betweenPrevPosition = null, betweenNextPosition = null;
//
//    public static AlignmentState makeInternalNode(final SAMRecord read, int readOffset,
//                                                  int genomeOffset, CigarElement cigarElement,
//                                                  int cigarElementCounter, final LinkedList<CigarElement> betweenPrevAndThis) {
//        final AlignmentState state = new AlignmentState(read, readOffset, genomeOffset, cigarElement, cigarElementCounter);
//        state.setBetweenPrevPosition(betweenPrevAndThis);
//        return state;
//    }
//
//
//
//    protected void update(final int readOffset, final int genomeOffset, final CigarElement cigarElement,
//                          final int cigarElementCounter, final LinkedList<CigarElement> betweenPrevAndThis,
//                          final CigarElement prevElement, final CigarElement nextElement) {
//        this.readOffset = readOffset;
//        this.genomeOffset = genomeOffset;
//        this.currentElement = cigarElement;
//        this.cigarElementCounter = cigarElementCounter;
//        this.betweenPrevPosition = betweenPrevAndThis;
//        this.prevElement = prevElement;
//        this.nextElement = nextElement;
//    }
//
//    // -----------------------------------------------------------------------------------------------
//    // Code for computing presence / absence of states in the prev / current / next
//    // -----------------------------------------------------------------------------------------------
//
////    public boolean isAfterDeletion() { return testOperator(getPrev(), CigarOperator.D); }
////    public boolean isBeforeDeletion() { return testOperator(getNext(), CigarOperator.D); }
////    public boolean isAfterInsertion() { return isAfter(getBetweenPrevPosition(), CigarOperator.I); }
////    public boolean isBeforeInsertion() { return isBefore(getBetweenNextPosition(), CigarOperator.I); }
////
////    public boolean isAfterSoftClip() { return isAfter(getBetweenPrevPosition(), CigarOperator.S); }
////    public boolean isBeforeSoftClip() { return isBefore(getBetweenNextPosition(), CigarOperator.S); }
////    public boolean isNextToSoftClip() { return isAfterSoftClip() || isBeforeSoftClip(); }
////
////    private boolean testOperator(final AlignmentState state, final CigarOperator op) {
////        return state != null && state.getCigarOperator() == op;
////    }
////
////    private boolean isAfter(final LinkedList<CigarElement> elements, final CigarOperator op) {
////        return ! elements.isEmpty() && elements.peekLast().getOperator() == op;
////    }
////
////    private boolean isBefore(final List<CigarElement> elements, final CigarOperator op) {
////        return ! elements.isEmpty() && elements.get(0).getOperator() == op;
////    }
//}
