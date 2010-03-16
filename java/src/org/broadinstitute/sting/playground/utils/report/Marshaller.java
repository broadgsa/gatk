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

package org.broadinstitute.sting.playground.utils.report;

import org.broadinstitute.sting.playground.utils.report.formats.OutputFormater;
import org.broadinstitute.sting.playground.utils.report.tags.*;
import org.broadinstitute.sting.utils.StingException;

import java.util.Collection;

/**
 * @author aaron
 *         <p/>
 *         Class Marshaller
 *         <p/>
 *         Given an analysis scan, an object, and a formatter, this class emits
 *         a report to disk (or other source).
 */
public class Marshaller {
    private String report;
    private OutputFormater formatter;

    /**
     * create a report
     * @param reportName the report name
     * @param formater the formater to use
     */
    public void createReport(String reportName, OutputFormater formater) {
        this.formatter = formater;
        report = reportName;
        formatter.startReport(report);
    }

    /**
     * add an analysis module to the output source
     * @param toMarshall the object to marshall
     */
    public void write(Object toMarshall) {
        AnalysisModuleScanner moduleScanner = new AnalysisModuleScanner(toMarshall);
        formatter.addAnalysis(moduleScanner.getAnalysis());
        for (Param p : moduleScanner.getParameters().keySet())
            try {
                formatter.addParam(p, moduleScanner.getParameters().get(p).get(toMarshall).toString());
            } catch (IllegalAccessException e) {
                throw new StingException("Unable to access variable " + moduleScanner.getParameters().get(p), e);
            }
        for (Datum d : moduleScanner.getData().keySet())
            try {
                formatter.addDatum(d, moduleScanner.getData().get(d).get(toMarshall).toString());
            } catch (IllegalAccessException e) {
                throw new StingException("Unable to access variable " + moduleScanner.getParameters().get(d), e);
            }
        for (Table t : moduleScanner.getTables().keySet())
            try {
                formatter.addTable(t, new TableContainer((Collection<Object>) moduleScanner.getTables().get(t).get(toMarshall),t.columns()));
            } catch (IllegalAccessException e) {
                throw new StingException("Unable to access variable " + moduleScanner.getParameters().get(t), e);
            }

    }

    public void endReport() {
        formatter.endReport();    
    }


}
