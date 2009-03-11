package org.broadinstitute.sting.utils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Collections;

import edu.mit.broad.picard.util.TabbedTextFileParser;

/**
 * Class for representing arbitrary reference ordered data sets
 *
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceOrderedData<ROD extends ReferenceOrderedDatum> implements Iterable<ROD> {
    private File file = null;
    private Class<ROD> type = null; // runtime type information for object construction

    public ReferenceOrderedData(File file, Class<ROD> type ) {
        this.file = file;
        this.type = type;
    }

    public RODIterator iterator() {
        return new RODIterator(new SimpleRODIterator());
    }

    // ----------------------------------------------------------------------
    //
    // Testing
    //
    // ----------------------------------------------------------------------
    public void testMe() {
        ReferenceOrderedDatum last = null;
        for ( ReferenceOrderedDatum rec : this ) {
            if ( last == null || ! last.getLocation().onSameContig(rec.getLocation()) ) {
                System.out.println(rec.toString());
            }
            last = rec;
        }
        System.exit(1);
    }

    // ----------------------------------------------------------------------
    //
    // Manipulations of all of the data
    //
    // ----------------------------------------------------------------------
    public ArrayList<ReferenceOrderedDatum> readAll() {
        ArrayList<ReferenceOrderedDatum> elts = new ArrayList<ReferenceOrderedDatum>();
        for ( ReferenceOrderedDatum rec : this ) {
            elts.add(rec);
        }
        elts.trimToSize();
        return elts;
    }

    public static void sortRODDataInMemory(ArrayList<ReferenceOrderedDatum> data) {
        Collections.sort(data);
    }

    public static void write(ArrayList<ReferenceOrderedDatum> data, File output) throws IOException {
        final FileWriter out = new FileWriter(output);

        for ( ReferenceOrderedDatum rec : data ) {
            out.write(rec.repl() + "\n");
        }

        out.close();
    }

    public boolean validateFile() throws Exception {
        ReferenceOrderedDatum last = null;
        for ( ReferenceOrderedDatum rec : this ) {
             if ( last != null && last.compareTo(rec) == 1 ) {
                 // It's out of order
                 throw new Exception("Out of order elements at \n" + last.toString() + "\n" + rec.toString());
             }
             last = rec;
        }
        return true;
    }

    public void indexFile() {
        // Fixme -- get access to the linear index system from Jim
    }

    // ----------------------------------------------------------------------
    //
    // Iteration
    //
    // ----------------------------------------------------------------------
    private class SimpleRODIterator implements Iterator<ROD> {
        //private WhitespaceTextFileParser parser = null;
        private TabbedTextFileParser parser = null;

        public SimpleRODIterator() {
            parser = new TabbedTextFileParser(true, file);
        }

        public boolean hasNext() {
            return parser.hasNext();
        }

        public ROD next() {
            String parts[] = parser.next();
            return parseLine(parts);
        }
 
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    public class RODIterator implements Iterator<ROD> {
        private PushbackIterator<ROD> it;
        private ROD prev = null;
        
        public RODIterator(SimpleRODIterator it) {
            this.it = new PushbackIterator<ROD>(it);
        }

        public boolean hasNext() { return it.hasNext(); }
        public ROD next() {
            prev = it.next();
            return prev; 
        }

        //
        // Seeks forward in the file until we reach (or cross) a record at contig / pos
        // If we don't find anything and cross beyond contig / pos, we return null
        // Otherwise we return the first object who's start is at pos
        //
        public ROD seekForward(final GenomeLoc loc) {
            return seekForward(loc.getContig(), loc.getStart());
        }

        protected ROD seekForward(final String contigName, final long pos) {
            final boolean DEBUG = false; 
        
            ROD result = null;
            
            if ( DEBUG ) System.out.printf("  *** starting seek to %s %d %s%n", contigName, pos, prev);
            while ( hasNext() ) {
                ROD current = next();
                //System.out.printf("    -> Seeking to %s %d AT %s %d%n", contigName, pos, current.getContig(), current.getStart());
                int strCmp = GenomeLoc.compareContigs( contigName, prev.getContig() );// contigName.compareTo( prev.getContig() );
                if ( strCmp == 0 ) {
                    // The contigs are equal
                    if ( current.getStart() > pos ) {
                        // There was nothing to find, push back next and return null
                        it.pushback(current);
                        break;
                    }
                    else if ( pos == current.getStart() ) {
                        // We found a record at contig / pos, return it
                        result = current;
                        break;
                    }
                }
                else if ( strCmp < 0 ) {
                    if ( DEBUG ) System.out.printf("    -> Jumped to contig %s%n", contigName);
                    // We've gone past the desired contig, break
                    it.pushback(current);
                    break;
                }
            }

            if ( DEBUG ) {
                if ( result == null )
                    ;
                    //System.out.printf("    --- seek result to %s %d is NULL%n", contigName, pos);
                else
                    System.out.printf("    ### Found %s %d%n", result.getContig(), result.getStart());
            }
            

             // we ran out of elements or found something
            return result;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    // ----------------------------------------------------------------------
    //
    // Parsing
    //
    // ----------------------------------------------------------------------
    ROD parseLine(final String[] parts) {
        //System.out.printf("Parsing GFFLine %s%n", Utils.join(" ", parts));
        try {
            ROD obj = type.newInstance();
            obj.parseLine(parts);
            return obj;
        } catch ( java.lang.InstantiationException e ) {
            System.out.println(e);
            return null; // wow, unsafe!
        } catch ( java.lang.IllegalAccessException e ) {
            System.out.println(e);
            return null; // wow, unsafe!
        }       
    }
}
