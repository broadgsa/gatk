package edu.mit.broad.sting.utils;

import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.picard.util.TabbedTextFileParser;

import java.io.File;
import java.io.InputStream;
import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.util.Iterator;
import java.util.HashMap;

/**
 * Class for representing arbitrary reference ordered data sets
 *
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceOrderedData implements Iterable<ReferenceOrderedDatum> {
    private File file = null;

    public ReferenceOrderedData(File file) {
        this.file = file;
    }

    // ----------------------------------------------------------------------
    //
    // Iteration
    //
    // ----------------------------------------------------------------------
    private class RODIterator implements Iterator<ReferenceOrderedDatum> {
        TabbedTextFileParser parser = null;
        public RODIterator() {
            parser = new TabbedTextFileParser(true, file);
        }

        public boolean hasNext() {
            return parser.hasNext();
        }

        public ReferenceOrderedDatum next() {
            String parts[] = parser.next();
            return parseGFFLine(parts);
        }

        public void remove () {
            throw new UnsupportedOperationException();
        }
    }

    public RODIterator iterator() {
        return new RODIterator();
    }

    // ----------------------------------------------------------------------
    //
    // Testing
    //
    // ----------------------------------------------------------------------
    public void testMe() {
        for ( ReferenceOrderedDatum rec : this ) {
            System.out.println(rec.toString());
        }
    }


    // ----------------------------------------------------------------------
    //
    // Parsing
    //
    // ----------------------------------------------------------------------
    ReferenceOrderedDatum parseGFFLine(final String[] parts) {
        //System.out.printf("Parsing GFFLine %s%n", Utils.join(" ", parts));

        final String contig = parts[0];
        final String source = parts[1];
        final String feature = parts[2];
        final long start = Long.parseLong(parts[3]);
        final long stop = Long.parseLong(parts[4]);

        double score = Double.NaN;
        if ( ! parts[5].equals(".") )
            score = Double.parseDouble(parts[5]);

        final String strand = parts[6];
        final String frame = parts[7];
        HashMap<String, String> attributes = null;
        return new ReferenceOrderedDatum(contig, source, feature, start, stop, score, strand, frame, attributes);
    }
}
