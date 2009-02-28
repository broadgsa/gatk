package edu.mit.broad.sam;

import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.util.NonDestructiveIterator;

/**
 * Wrapper around SAMRecord iterator that skips over non-primary elements.
 */
public class NotPrimarySkippingIterator {
    private final NonDestructiveIterator<SAMRecord, CloseableIterator<SAMRecord>> it;

    public NotPrimarySkippingIterator(final CloseableIterator<SAMRecord> underlyingIt) {
        it = new NonDestructiveIterator<SAMRecord, CloseableIterator<SAMRecord>>(underlyingIt);
        skipAnyNotprimary();
    }

    public boolean hasCurrent() {
        return it.hasCurrent();
    }

    public SAMRecord getCurrent() {
        assert(hasCurrent());
        return it.getCurrent();
    }

    public boolean advance() {
        it.advance();
        skipAnyNotprimary();
        return hasCurrent();
    }

    private void skipAnyNotprimary() {
        while (it.hasCurrent() && it.getCurrent().getNotPrimaryAlignmentFlag()) {
            it.advance();
        }
    }
}
