/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.util;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.sam.util.AsciiLineReader;
import edu.mit.broad.sam.util.RuntimeIOException;

import java.io.*;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

/**
 * TextFileParser which reads a single text file.
 *
 * @author Kathleen Tibbetts
 */
public class BasicTextFileParser extends AbstractTextFileParser
{
    // TODO: Replace this with AsciiStreamReader when Alec creates it.
    private AsciiLineReader reader;
    private final ArrayList<File> files = new ArrayList<File>();
    String currentFileName = null;

    /**
     * Constructor.  Opens up a buffered reader and reads the first line.
     *
     * @param files  the file(s) to parse, in order
     */
    public BasicTextFileParser(boolean treatGroupedDelimitersAsOne, File... files) {
        if (files.length == 0) {
            throw new IllegalArgumentException("At least one file must be specified.");
        }
        this.files.addAll(Arrays.asList(files));
        File f = this.files.remove(0);
        currentFileName = f.getAbsolutePath();
        reader = new AsciiLineReader(IoUtil.openFileForReading(f));
        this.setTreatGroupedDelimitersAsOne(treatGroupedDelimitersAsOne);
    }

    /**
     * Constructor.  In addition to opening and priming the files, it sets the number of
     * whitespace-separated "words" per line.
     *
     * @param files      the file(s) to parse
     * @param wordCount number of whitespace-separated "words" per line
     */
    public BasicTextFileParser(boolean treatGroupedDelimitersAsOne, int wordCount, File... files) {
        this(treatGroupedDelimitersAsOne, files);
        setWordCount(wordCount);
    }
    /**
     * Workhorse method that reads the next line from the underlying reader
     *
     * @return  String or null if there is no next line
     */
    protected byte[] readNextLine()
    {
        try {
            String line = reader.readLine();
            if (line != null) {
                return line.getBytes();
            }
            if (files.size() > 0) {
                currentFileName = files.get(0).getAbsolutePath();
                reader = new AsciiLineReader(IoUtil.openFileForReading(files.remove(0)));
                return readNextLine();
            }
            return null;
        }
        catch(RuntimeIOException ioe) {
            throw new PicardException("Error reading from file " + currentFileName, ioe);
        }
    }

    /**
     * Closes the underlying stream
     */
    public void close() {
        if (reader != null)  {
            reader.close();
        }
    }

    /**
     * Gets the name of the file being parsed
     *
     * @return  the name of the file being parsed
     */
    protected String getFileName() {
        return this.currentFileName;
    }
}
