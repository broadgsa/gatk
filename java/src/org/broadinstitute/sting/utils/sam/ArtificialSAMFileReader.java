package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.io.InputStream;
import java.io.ByteArrayInputStream;
import java.io.UnsupportedEncodingException;
import java.util.List;
import java.util.Arrays;
import java.util.Iterator;
import java.util.ArrayList;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
/**
 * User: hanna
 * Date: Jun 11, 2009
 * Time: 9:35:31 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Pass specified reads into the given walker.
 */

public class ArtificialSAMFileReader extends SAMFileReader {
    /**
     * Backing data store of reads.
     */
    private List<SAMRecord> reads = null;

    /**
     * Construct an artificial SAM file reader.
     */
    public ArtificialSAMFileReader(SAMRecord... reads) {
        super( createEmptyInputStream(),true );
        this.reads = Arrays.asList(reads);
    }

    /**
     * @{inheritDoc}
     */
    @Override
    public CloseableIterator<SAMRecord> query(final String sequence, final int start, final int end, final boolean contained) {
        GenomeLoc region = GenomeLocParser.createGenomeLoc(sequence, start, end);
        List<SAMRecord> coveredSubset = new ArrayList<SAMRecord>();

        for( SAMRecord read: reads ) {
            GenomeLoc readPosition = GenomeLocParser.createGenomeLoc(read);
            if( contained && region.containsP(readPosition) ) coveredSubset.add(read);
            else if( !contained && readPosition.overlapsP(region) ) coveredSubset.add(read);
        }

        final Iterator<SAMRecord> iterator = coveredSubset.iterator();
        return new CloseableIterator<SAMRecord>() {
            public boolean hasNext() { return iterator.hasNext(); }
            public SAMRecord next() { return iterator.next(); }
            public void close() {}
            public void remove() { iterator.remove(); } 
        };
    }

    /**
     * Builds an empty input stream for faking out the sam file reader.
     * Derive it from a string so that, in the future, it might be possible
     * to fake the text of a sam file from samtools output, et.c
     * @return Stream that returns no characters.
     */
    private static InputStream createEmptyInputStream() {
        try {
            byte[] byteArray = "".getBytes("ISO-8859-1");
            return new ByteArrayInputStream(byteArray);
        }
        catch( UnsupportedEncodingException ex ) {
            throw new StingException("Unable to build empty input stream",ex);
        }
    }
}
