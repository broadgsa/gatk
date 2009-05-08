package org.broadinstitute.sting.gatk.refdata;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.util.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.xReadLines;
import org.broadinstitute.sting.utils.Utils;
import org.apache.log4j.Logger;

/**
 * Class for representing arbitrary reference ordered data sets
 *
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceOrderedData<ROD extends ReferenceOrderedDatum> implements Iterable<ROD> {
    private String name;
    private File file = null;
    private Class<ROD> type = null; // runtime type information for object construction

    // ----------------------------------------------------------------------
    //
    // Static ROD type management
    //
    // ----------------------------------------------------------------------
    public static class RODBinding {
        public final String name;
        public final Class<? extends ReferenceOrderedDatum> type;
        public RODBinding(final String name, final Class<? extends ReferenceOrderedDatum> type) {
            this.name = name;
            this.type = type;
        }
    }

    public static HashMap<String, RODBinding> Types = new HashMap<String, RODBinding>();
    public static void addModule(final String name, final Class<? extends ReferenceOrderedDatum> rodType) {
        System.out.printf("* Adding rod class %s%n", name);
        Types.put(name.toLowerCase(), new RODBinding(name, rodType));
    }

    static {
        // All known ROD types
        addModule("GFF", rodGFF.class);
        addModule("dbSNP", rodDbSNP.class);
        addModule("HapMapAlleleFrequencies", HapMapAlleleFrequenciesROD.class);
        addModule("SAMPileup", rodSAMPileup.class);
    }


    /**
     * Parse the ROD bindings.  These are of the form of a single list of strings, each triplet of the
     * form <name>,<type>,<file>.  After this function, the List of RODs contains new RODs bound to each of
     * name, of type, ready to read from the file.  This function does check for the strings to be well formed
     * and such.
     *
     * @param logger
     * @param bindings
     * @param rods
     */
    public static void parseBindings(Logger logger, ArrayList<String> bindings, List<ReferenceOrderedData<? extends ReferenceOrderedDatum> > rods)
    {
        // Loop over triplets
        for( String binding: bindings ) {
            String[] bindingTokens = binding.split(",");
            logger.info("Processing ROD bindings: " + bindings.size() + " -> " + Utils.join(" : ", bindingTokens));                                
            if( bindingTokens.length != 3 )
                Utils.scareUser(String.format("Invalid ROD specification: requires triplets of <name>,<type>,<file> but got %s", Utils.join(",", bindings)));

            final String name = bindingTokens[0];
            final String typeName = bindingTokens[1];
            final String fileName = bindingTokens[2];

            ReferenceOrderedData<?> rod = parse1Binding(logger, name, typeName, fileName);

            // check that we're not generating duplicate bindings
            for ( ReferenceOrderedData rod2 : rods )
                if ( rod2.getName().equals(rod.getName()) )
                    Utils.scareUser(String.format("Found duplicate rod bindings", rod.getName()));

            rods.add(rod);            
        }
    }

    /**
     * Helpful function that parses a single triplet of <name> <type> <file> and returns the corresponding ROD with
     * <name>, of type <type> that reads its input from <file>.
     * 
     * @param logger
     * @param trackName
     * @param typeName
     * @param fileName
     * @return
     */
    private static ReferenceOrderedData<?> parse1Binding( Logger logger, final String trackName, final String typeName, final String fileName )
    {
        // Gracefully fail if we don't have the type
        if ( ReferenceOrderedData.Types.get(typeName.toLowerCase()) == null )
            Utils.scareUser(String.format("Unknown ROD type: %s", typeName));

        // Lookup the type
        Class rodClass = ReferenceOrderedData.Types.get(typeName.toLowerCase()).type;

        // Create the ROD
        ReferenceOrderedData<?> rod = new ReferenceOrderedData<ReferenceOrderedDatum>(trackName.toLowerCase(), new File(fileName), rodClass );
        logger.info(String.format("Created binding from %s to %s of type %s", trackName.toLowerCase(), fileName, rodClass));
        return rod;
    }

    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public ReferenceOrderedData(final String name, File file, Class<ROD> type ) {
        this.file = file;
        this.type = type;
        this.name = name;
    }

    public String getName() { return name; }

    public RODIterator iterator() {
        return new RODIterator(new SimpleRODIterator());
    }

    // ----------------------------------------------------------------------
    //
    // Testing
    //
    // ----------------------------------------------------------------------
    public void testMe() {
        for ( ReferenceOrderedDatum rec : this ) {
            System.out.println(rec.toString());

            rodGFF gff = (rodGFF)rec;
            String[] keys = {"LENGTH", "ALT", "FOBARBAR"};
            for ( String key : keys) {
                System.out.printf("  -> %s is (%s)%n", key, gff.containsAttribute(key) ? gff.getAttribute(key) : "none");
            }
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
        private xReadLines parser = null;

        public SimpleRODIterator() {
            try {
                parser = new xReadLines(file);
            } catch ( FileNotFoundException e ) {
                Utils.scareUser("Couldn't open file: " + file);
            }
        }

        public boolean hasNext() {
            //System.out.printf("Parser has next: %b%n", parser.hasNext());
            return parser.hasNext();
        }

        public ROD next() {
            final String line = parser.next();
            //System.out.printf("Line is %s%n", line);
            String parts[] = line.split("\t");
            return parseLine(parts);
        }
 
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

//    private class SimpleRODIterator implements Iterator<ROD> {
//        //private WhitespaceTextFileParser parser = null;
//        private TabbedTextFileParser parser = null;
//
//        public SimpleRODIterator() {
//            parser = new TabbedTextFileParser(true, file);
//        }
//
//        public boolean hasNext() {
//            return parser.hasNext();
//        }
//
//        public ROD next() {
//            String parts[] = parser.next();
//            return parseLine(parts);
//        }
//
//        public void remove() {
//            throw new UnsupportedOperationException();
//        }
//    }

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
            //ROD obj = type.newInstance();
            Constructor<ROD> c = type.getConstructor(String.class);
            ROD obj = (ROD)c.newInstance(name);
            obj.parseLine(parts);
            return obj;
        } catch ( java.lang.InstantiationException e ) {
            System.out.println(e);
            return null; // wow, unsafe!
        } catch ( java.lang.IllegalAccessException e ) {
            System.out.println(e);
            return null; // wow, unsafe!
        } catch ( java.lang.NoSuchMethodException e ) {
            System.out.println(e);
            return null; // wow, unsafe!
        } catch ( InvocationTargetException e ) {
            System.out.println(e);
            return null; // wow, unsafe!
        }
    }
}
