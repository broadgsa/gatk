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


/**
 * 
 * @author aaron 
 * 
 * Class OutputFormater
 *
 * The formatter defines a perticular output format style.  The formatter methods
 * are required to be called in the following order:
 *
 * 1. startReport()
 * 2. for each analysis module:
 *      addAnalysis()
 *          with many addParam(), addTable(), and addDatum() in no particular order promised
 * 3. endReport()
 */
public interface OutputFormater {
    /**
     * start the report, given the report tag
     * @param reportName the report name.
     */
    public void startReport(String reportName);

    /**
     * add an analysis to the report
     * @param analysisTag the analysis tag
     */
    public void addAnalysis(Analysis analysisTag);

    /**
     * add a parameter value to the output location
     * @param value
     */
    public void addParam(Param paramTag, String value);

    /**
     * add a parameter value to the output location
     * @param value
     */
    public void addTable(Table tableTag, TableContainer value);

    /**
     * add a datum object to the output location
     * @param value
     */
    public void addDatum(Datum datumTag, String value);

    /**
     * end the report 
     */
    public void endReport();
}
