package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.*;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Class for representing arbitrary reference ordered data sets
 * <p/>
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceOrderedData<ROD extends ReferenceOrderedDatum> implements Iterable<ReferenceOrderedDatum> {
    private String name;
    private File file = null;
//    private String fieldDelimiter;

    /** Header object returned from the datum */
//    private Object header = null;

    private Class<ROD> type = null; // runtime type information for object construction

    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(ReferenceOrderedData.class);

    /**
     * given an existing file, open it and append all the valid triplet lines to an existing list
     *
     * @param rodTripletList the list of existing triplets
     * @param filename       the file to attempt to extract ROD triplets from
     */
    protected static void extractRodsFromFile(List<String> rodTripletList, String filename) {
        BufferedReader str;
        try {
            str = new BufferedReader(new FileReader(new File(filename)));
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(new File(filename), "Unable to load the ROD input file", e);
        }
        String line = "NO LINES READ IN";
        try {
            while ((line = str.readLine()) != null) {
                if (line.matches(".+,.+,.+")) rodTripletList.add(line.trim());
                else logger.warn("the following file line didn't parsing into a triplet -> " + line);
            }
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(new File(filename), "Failed reading the input rod file; last line read was " + line, e);
        }
    }


    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public ReferenceOrderedData(final String name, File file, Class<ROD> type ) {
        this.name = name;
        this.file = file;
        this.type = type;
//        this.header = initializeROD(name, file, type);
//        this.fieldDelimiter = newROD(name, type).delimiterRegex();
    }

    public String getName() { return name; }

    public File getFile() { return file; }

    public Class<ROD> getType() { return type; }

    /**
     * Special equals override to see if this ROD is compatible with the given
     * name and type.  'Compatible' means that this ROD has the name that's passed
     * in and its data can fit into the container specified by type.
     *
     * @param name Name to check.
     * @param type Type to check.
     *
     * @return True if these parameters imply this rod.  False otherwise.
     */
    public boolean matches(String name, Class<? extends ReferenceOrderedDatum> type) {
        return this.name.equals(name) && type.isAssignableFrom(this.type);
    }

    public Iterator<ReferenceOrderedDatum> iterator() {
        Iterator<ReferenceOrderedDatum> it;
        try {
            Method m = type.getDeclaredMethod("createIterator", String.class, java.io.File.class);
            it = (Iterator<ReferenceOrderedDatum>) m.invoke(null, name, file);
        } catch (java.lang.NoSuchMethodException e) {
            it = new RODRecordIterator(file,name,type);
        } catch (java.lang.NullPointerException e) {
            throw new RuntimeException(e);
        } catch (java.lang.SecurityException e) {
            throw new RuntimeException(e);
        } catch (java.lang.IllegalAccessException e) {
            throw new RuntimeException(e);
        } catch (java.lang.IllegalArgumentException e) {
            throw new RuntimeException(e);
        } catch (java.lang.reflect.InvocationTargetException e) {
            throw new RuntimeException(e);
        }
  //      return new RODIterator<ROD>(it);
        return it;
    }

    // ----------------------------------------------------------------------
    //
    // Manipulations of all of the data
    //
    // ----------------------------------------------------------------------

    public static void write(ArrayList<ReferenceOrderedDatum> data, File output) throws IOException {
        final FileWriter out = new FileWriter(output);

        for (ReferenceOrderedDatum rec : data) {
            out.write(rec.repl() + "\n");
        }

        out.close();
    }


}
