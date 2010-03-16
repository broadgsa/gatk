/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.utils.report.formats;

import org.broadinstitute.sting.playground.utils.report.TableContainer;
import org.broadinstitute.sting.playground.utils.report.tags.*;
import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class SimpleTextOutputer
 *
 * a very small demo of text output, just to use for testing and development
 */
public class SimpleTextOutputer implements OutputFormater {
    protected PrintWriter out;
    protected String terminator = "\n";
    protected String spacer = "\t\t\t";
    /**
     * create a simple text output format given the output file
     * @param outputFile
     */
    public SimpleTextOutputer(File outputFile) {
        try {
            FileOutputStream stream = new FileOutputStream(outputFile);
            out = new PrintWriter(stream);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to find file ",e);
        }
    }

    /**
     * start the report, given the report tag
     *
     * @param reportTag the report tag, which describes the basics of the report
     */
    @Override
    public void startReport(String reportTag) {
        out.write("----------------" + reportTag + "----------------" + terminator);
    }

    /**
     * add an analysis to the report
     *
     * @param analysisTag the analysis tag
     */
    @Override
    public void addAnalysis(Analysis analysisTag) {
        out.write("----------------" + analysisTag.name() + "----------------" + terminator);
        out.write(String.format("Name%s%s%s",spacer,analysisTag.name(),terminator));
        out.write(String.format("Description\t\t%s%s",analysisTag.description(),terminator));
        out.write(String.format("Version%s%s%s",spacer,analysisTag.version(),terminator));
    }

    /**
     * add a parameter value to the output location
     *
     * @param value
     */
    @Override
    public void addParam(Param paramTag, String value) {
        out.write(String.format("Param %s%s%s%s",paramTag.name(),spacer,value,terminator));
    }

    /**
     * add a parameter value to the output location
     *
     * @param value
     */
    @Override
    public void addTable(Table tableTag, TableContainer value) {
        out.write(tableTag.name());
        boolean first = false;
        for (List<Object> row : value.toRows()) {
            if (!first) first = true;
            else out.write("\t");
            out.write("\t\t");
            for (Object obj: row) out.write(obj+",");
            out.write(terminator);
        }

    }

    /**
     * add a datum object to the output location
     *
     * @param value
     */
    @Override
    public void addDatum(Datum datumTag, String value) {
        out.write("Datum "+datumTag.name()+spacer+value+terminator);
    }

    /**
     * end the report
     *     
     */
    @Override
    public void endReport() {
        out.write("----------------" + "----------------" + terminator);
        out.close();
    }
}
