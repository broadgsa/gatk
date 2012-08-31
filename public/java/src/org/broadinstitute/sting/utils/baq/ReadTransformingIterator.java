package org.broadinstitute.sting.utils.baq;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Iterator;

/**
 * Iterator that applies a ReadTransformer to a stream of reads
 */
public class ReadTransformingIterator implements StingSAMIterator {
    private final StingSAMIterator it;
    private final ReadTransformer transformer;

    /**
     * Creates a new ReadTransforming iterator
     */
    @Requires({"it != null", "transformer != null", "transformer.isInitialized()"})
    public ReadTransformingIterator(final StingSAMIterator it, final ReadTransformer transformer) {
        if ( ! transformer.isInitialized() )
            throw new IllegalStateException("Creating a read transformer stream for an uninitialized read transformer: " + transformer);
        if ( transformer.getApplicationTime() == ReadTransformer.ApplicationTime.FORBIDDEN )
            throw new IllegalStateException("Creating a read transformer stream for a forbidden transformer " + transformer);

        this.it = it;
        this.transformer = transformer;
    }

    @Requires("hasNext()")
    @Ensures("result != null")
    public SAMRecord next()     {
        final GATKSAMRecord read = (GATKSAMRecord)it.next();
        return transformer.apply(read);
    }

    public boolean hasNext()    { return this.it.hasNext(); }
    public void remove()        { throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!"); }
    public void close()         { it.close(); }
    public Iterator<SAMRecord> iterator() { return this; }
}
