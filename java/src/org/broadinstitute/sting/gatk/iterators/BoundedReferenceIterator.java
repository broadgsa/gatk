package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.Iterator;

import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

/**
 *
 * User: aaron
 * Date: Apr 2, 2009
 * Time: 2:12:12 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 2, 2009
 * <p/>
 * Class BoundedReferenceIterator
 * <p/>
 * This class is a decorator class from Reference Iterator (though it is constrained
 * by the fact that referenceIterator.seekForwardOffset explicitly returns a referenceIterator
 * for now).
 * <p/>
 * TODO: Fix the underlying iterator and this class to model a real decorator pattern
 */
public class BoundedReferenceIterator implements Iterator<Character> {
    private ReferenceSequence referenceSequence = null;

    // The start and end point of the read relative to the chromosome.
    private long start;
    private long end;

    // 0-based index into the base array.
    private int position;
    // the location to screen over

    /**
     * Default constructor
     *
     * @param referenceSequenceFile sequence file over which to iterate
     * @param loc subsequence of the reference over which to iterate.
     */
    public BoundedReferenceIterator(IndexedFastaSequenceFile referenceSequenceFile, GenomeLoc loc) {
        start = loc.getStart();
        end = loc.getStop();
        position = 0;

        referenceSequence = referenceSequenceFile.getSubsequenceAt( loc.getContig(), position, end );
    }


    // our adapted next function
    public boolean hasNext() {
        return (start + position) > end;
    }

    public Character next() {
        return StringUtil.bytesToString( referenceSequence.getBases(), position++, 1 ).charAt(0);
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a reference iterator");
    }


}
