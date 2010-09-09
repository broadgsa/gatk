/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.exceptions.UserError;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;

import java.util.Iterator;
import java.util.regex.Pattern;
import java.io.FileNotFoundException;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Constructor;

/**
 * This is a low-level iterator designed to provide system-wide generic support for reading record-oriented data
 * files. The only assumption made is that every line in the file provides a complete and separate data record. The records
 * can be associated with coordinates or coordinate intervals, there can be one or more records associated with a given
 * position/interval, or intervals can overlap. The records must be comprised of delimited fields, but the format is
 * otherwise free. For any specific line-based data format, an appropriate implementation of ReferenceOrderedDatum must be
 * provided that is capable of parsing itself from a single line of data. This implementation will be used,
 * through reflection mechanism, as a callback to do all the work.
 *
 * The model is, hence, as follows:
 *
 *  String dataRecord <---> RodImplementation ( ::parseLine(dataRecord.split(delimiter)) is aware of the format and fills
 *                                              an instance of RodImplementation with data values from dataRecord line).
 *
 *
 *  instantiation of RODRecordIterator(dataFile, trackName, RodImplementation.class) will immediately provide an iterator
 * that walks along the dataFile line by line, and on each call to next() returns a new RodImplementation object
 * representing a single line (record) of data. The returned object will be initialized with "track name" trackName -
 * track names (as returned by ROD.getName()) are often used in other parts of the code to distinguish between
 * multiple streams of (possibly heterogeneous) annotation data bound to an application.
 *
 * This generic iterator skips and ignores a) empty lines, b) lines starting with '#' (comments): they are never sent back
 * to the ROD implementation class for processing.
 *
 * This iterator does not actually check if the ROD records (lines) in the file are indeed ordedered by coordinate,
 * and it does not depend on such an order as it still implements a low-level line-based traversal of the data. Higher-level
 * iterators/wrappers will perform all the necessary checks.
 *
 * Note: some data formats/ROD implementations may require a header line in the file. In this case the current (ugly)
 * mechanism is as follows:
 * 1) rod implementation's ::initialize(file) method should be able to open the file, find and read the header line
 *    and return the header object (to be kept by the iterator)
 * 2) rod implementation's ::parseLine(header,line) method should be capable of making use of that saved header object now served to it
 * and
 * 3) ::parseLine(header,line) should be able to recognize the original header line in the file and skip it (after ROD's initialize()
 *    method is called, the iterator will re-open the file and start reading it from the very beginning; there is no
 *    other way, except for "smart" ::parseLine(), to avoid reading in the header line as "data"). 
 *
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 10, 2009
 * Time: 1:22:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class RODRecordIterator<ROD extends ReferenceOrderedDatum> implements Iterator<ROD> {

    private PushbackIterator<String> reader;

    // stores name of the track this iterator reads (will be also returned by getName() of ROD objects
    // generated by this iterator)
    private String name;

    // we keep the file object, only to use file name in error reports
    private File file;

    // rod type; this is what we will instantiate for RODs at runtime
    private Class<ROD> type;

    private Object header = null; // Some RODs may use header

    // field delimiter in the file. Should it be the job of the iterator to split the lines though? RODs can do that!
    private String fieldDelimiter;

    // constructor for the ROD objects we are going to return. Constructor that takes the track name as its single arg is required.
    private Constructor<ROD> named_constructor;

    // keep track of the lines we are reading. used for error messages only.
    private long linenum = 0;

    private boolean allow_empty = true;
    private boolean allow_comments = true;
    public static Pattern EMPTYLINE_PATTERN = Pattern.compile("^\\s*$");

    public RODRecordIterator(File file, String name, Class<ROD> type) {
        try {
            reader = new PushbackIterator<String>(new XReadLines(file));
        } catch (FileNotFoundException e) {
            throw new UserError.CouldNotReadInputFile(file, e);
        }
        this.file = file;
        this.name = name;
        this.type = type;
        try {
            named_constructor = type.getConstructor(String.class);
        }
        catch (java.lang.NoSuchMethodException e) {
            throw new GATKException("ROD class "+type.getName()+" does not have constructor that accepts a single String argument (track name)");
        }
        ROD rod = instantiateROD(name);
        fieldDelimiter = rod.delimiterRegex(); // get delimiter from the ROD itself
        try {
            header = rod.initialize(file);
        } catch (FileNotFoundException e) {
            throw new UserError.CouldNotReadInputFile(file, "ROD "+type.getName() + " failed to initialize properly from file "+file);
        }

    }


    /**
     * Returns <tt>true</tt> if the iteration has more elements. (In other
     * words, returns <tt>true</tt> if <tt>next</tt> would return an element
     * rather than throwing an exception.)
     *
     * @return <tt>true</tt> if the iterator has more elements.
     */
    public boolean hasNext() {
        if ( allow_empty || allow_comments ) {
            while ( reader.hasNext() ) {
                String line = reader.next();
                if ( allow_empty && EMPTYLINE_PATTERN.matcher(line).matches() ) continue; // skip empty line
                if ( allow_comments && line.charAt(0) == '#' ) continue; // skip comment lines
                // the line is not empty and not a comment line, so we have next after all
                reader.pushback(line);
                return true;
            }
            return false; // oops, we end up here if there's nothing left
        } else {
            return reader.hasNext();
        }
    }

    /**
     * Returns the next valid ROD record in the file, skipping empty and comment lines.
     *
     * @return the next element in the iteration.
     * @throws java.util.NoSuchElementException
     *          iteration has no more elements.
     */
    public ROD next() {
        ROD n = null;
        boolean parsed_ok = false;
        String line ;

        while ( ! parsed_ok && reader.hasNext() ) {
            line = reader.next();
            linenum++;
            while (   allow_empty && EMPTYLINE_PATTERN.matcher(line).matches() ||
                      allow_comments && line.charAt(0) == '#' ) {
                if ( reader.hasNext() ) {
                    line = reader.next();
                    linenum++;
                } else {
                    line = null;
                    break;
                }
            }

            if ( line == null ) break; // if we ran out of lines while skipping empty lines/comments, then we are done
            
            String parts[] = line.split(fieldDelimiter);

            try {
                 n = instantiateROD(name);
                 parsed_ok =  n.parseLine(header,parts) ;
             }
             catch ( Exception e ) {
                 throw new UserError.MalformedFile(file, "Failed to parse ROD data ("+type.getName()+") from file "+ file + " at line #"+linenum+
                         "\nOffending line: "+line+
                         "\nReason ("+e.getClass().getName()+")", e);
             }
        }


        return n;
    }

    /**
     * Removes from the underlying collection the last element returned by the
     * iterator (optional operation).  This method can be called only once per
     * call to <tt>next</tt>.  The behavior of an iterator is unspecified if
     * the underlying collection is modified while the iteration is in
     * progress in any way other than by calling this method.
     *
     * @throws UnsupportedOperationException if the <tt>remove</tt>
     *                                       operation is not supported by this Iterator.
     * @throws IllegalStateException         if the <tt>next</tt> method has not
     *                                       yet been called, or the <tt>remove</tt> method has already
     *                                       been called after the last call to the <tt>next</tt>
     *                                       method.
     */
    public void remove() {
        throw new UnsupportedOperationException("remove() operation is not supported by RODRecordIterator");
    }

    /** Instantiates appropriate implementation of the ROD used by this iteratot. The 'name' argument is the name
     * of the ROD track.
     * @param name
     * @return
     */
    private ROD instantiateROD(final String name) {
        try {
            return (ROD) named_constructor.newInstance(name);
        } catch (java.lang.InstantiationException e) {
            throw new GATKException("Failed to instantiate ROD object of class "+type.getName());
        } catch (java.lang.IllegalAccessException e) {
            throw new GATKException("Access violation attempt while instantiating ROD object of class "+type.getName());
        } catch (InvocationTargetException e) {
            throw new GATKException("InvocationTargetException: Failed to instantiate ROD object of class "+type.getName());
        }
    }

}
