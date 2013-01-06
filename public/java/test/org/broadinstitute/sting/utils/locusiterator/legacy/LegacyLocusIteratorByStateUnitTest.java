package org.broadinstitute.sting.utils.locusiterator.legacy;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.locusiterator.legacy.LegacyLocusIteratorByState;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

class FakeCloseableIterator<T> implements CloseableIterator<T> {
    Iterator<T> iterator;

    public FakeCloseableIterator(Iterator<T> it) {
        iterator = it;
    }

    @Override
    public void close() {}

    @Override
    public boolean hasNext() {
        return iterator.hasNext();
    }

    @Override
    public T next() {
        return iterator.next();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Don't remove!");
    }
}


final class LIBS_position {

    SAMRecord read;

    final int numOperators;
    int currentOperatorIndex = 0;
    int currentPositionOnOperator = 0;
    int currentReadOffset = 0;

    boolean isBeforeDeletionStart = false;
    boolean isBeforeDeletedBase = false;
    boolean isAfterDeletionEnd = false;
    boolean isAfterDeletedBase = false;
    boolean isBeforeInsertion = false;
    boolean isAfterInsertion = false;
    boolean isNextToSoftClip = false;

    boolean sawMop = false;

    public LIBS_position(final SAMRecord read) {
        this.read = read;
        numOperators = read.getCigar().numCigarElements();
    }

    public int getCurrentReadOffset() {
        return Math.max(0, currentReadOffset - 1);
    }

    /**
     * Steps forward on the genome.  Returns false when done reading the read, true otherwise.
     */
    public boolean stepForwardOnGenome() {
        if ( currentOperatorIndex == numOperators )
            return false;

        CigarElement curElement = read.getCigar().getCigarElement(currentOperatorIndex);
        if ( currentPositionOnOperator >= curElement.getLength() ) {
            if ( ++currentOperatorIndex == numOperators )
                return false;

            curElement = read.getCigar().getCigarElement(currentOperatorIndex);
            currentPositionOnOperator = 0;
        }

        switch ( curElement.getOperator() ) {
            case I: // insertion w.r.t. the reference
                if ( !sawMop )
                    break;
            case S: // soft clip
                currentReadOffset += curElement.getLength();
            case H: // hard clip
            case P: // padding
                currentOperatorIndex++;
                return stepForwardOnGenome();

            case D: // deletion w.r.t. the reference
            case N: // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                currentPositionOnOperator++;
                break;

            case M:
            case EQ:
            case X:
                sawMop = true;
                currentReadOffset++;
                currentPositionOnOperator++;
                break;
            default:
                throw new IllegalStateException("No support for cigar op: " + curElement.getOperator());
        }

        final boolean isFirstOp = currentOperatorIndex == 0;
        final boolean isLastOp = currentOperatorIndex == numOperators - 1;
        final boolean isFirstBaseOfOp = currentPositionOnOperator == 1;
        final boolean isLastBaseOfOp = currentPositionOnOperator == curElement.getLength();

        isBeforeDeletionStart = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.D, isLastOp, isLastBaseOfOp);
        isBeforeDeletedBase = isBeforeDeletionStart || (!isLastBaseOfOp && curElement.getOperator() == CigarOperator.D);
        isAfterDeletionEnd = isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.D, isFirstOp, isFirstBaseOfOp);
        isAfterDeletedBase  = isAfterDeletionEnd || (!isFirstBaseOfOp && curElement.getOperator() == CigarOperator.D);
        isBeforeInsertion = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.I, isLastOp, isLastBaseOfOp)
                || (!sawMop && curElement.getOperator() == CigarOperator.I);
        isAfterInsertion = isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.I, isFirstOp, isFirstBaseOfOp);
        isNextToSoftClip = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.S, isLastOp, isLastBaseOfOp)
                || isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.S, isFirstOp, isFirstBaseOfOp);

        return true;
    }

    private static boolean isBeforeOp(final Cigar cigar,
                                      final int currentOperatorIndex,
                                      final CigarOperator op,
                                      final boolean isLastOp,
                                      final boolean isLastBaseOfOp) {
        return  !isLastOp && isLastBaseOfOp && cigar.getCigarElement(currentOperatorIndex+1).getOperator() == op;
    }

    private static boolean isAfterOp(final Cigar cigar,
                                     final int currentOperatorIndex,
                                     final CigarOperator op,
                                     final boolean isFirstOp,
                                     final boolean isFirstBaseOfOp) {
        return  !isFirstOp && isFirstBaseOfOp && cigar.getCigarElement(currentOperatorIndex-1).getOperator() == op;
    }
}
