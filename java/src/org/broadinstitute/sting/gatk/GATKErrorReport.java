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

package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.utils.StingException;
import org.simpleframework.xml.Element;
import org.simpleframework. xml.ElementList;
import org.simpleframework.xml.Serializer;
import org.simpleframework.xml.core.Persister;
import org.simpleframework.xml.stream.Format;
import org.simpleframework.xml.stream.HyphenStyle;

import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class GATKErrorReport
 *         <p/>
 *         A basic description of what went wrong during a run of the GATK,
 *         what the parameters were, etc.
 */
public class GATKErrorReport {
    // the listing of the fields is somewhat important; this is the order that the simple XML will output them    
    @ElementList(required = true, name = "gatk_header_Information")
    private static List<String> mGATKHeader;

    @Element(required = false, name = "exception")
    private final ExceptionToXML mException;

    @Element(required = true, name = "date_time")
    private final String mDateTime;

    @Element(required = true, name = "argument_collection")
    private final GATKArgumentCollection mCollection;

    @Element(required = true, name = "working_directory")
    private static String currentPath;

    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH.mm.ss");

    static {
        GATKErrorReport.mGATKHeader = CommandLineGATK.createApplicationHeader();
        currentPath = System.getProperty("user.dir");
    }

    public GATKErrorReport(Exception e, GATKArgumentCollection collection) {
        this.mCollection = collection;
        this.mException = new ExceptionToXML(e);
        java.util.Date date = new java.util.Date();
        mDateTime = dateFormat.format(date);
    }

    public void reportToStream(PrintStream stream) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        try {
            serializer.write(this, stream);
        } catch (Exception e) {
            throw new StingException("Failed to marshal the data to the file " + stream, e);
        }
    }

    class ExceptionToXML {
        @ElementList(required = false)
        final List<String> exceptionDetails = new ArrayList<String>();

        public ExceptionToXML(Exception e) {
            for (StackTraceElement element : e.getStackTrace()) {
                exceptionDetails.add(element.toString());
            }
        }
    }
}
