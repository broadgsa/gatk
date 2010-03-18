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

package org.broadinstitute.sting.playground.utils.report.chunks;

import org.broadinstitute.sting.playground.utils.report.AnalysisModuleScanner;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;


/**
 * @author aaron
 *         <p/>
 *         Class AnalysisChunk
 *         <p/>
 *         The basic data we've decorated an analysis with
 */
public class AnalysisChunk implements Chunk {
    public String name;
    public String description;
    public String version;
    public String className;
    /**
     * extract the analysis data from an analysis annotation using the module scanner.
     * @param moduleScanner the module scanner
     */
    public AnalysisChunk(AnalysisModuleScanner moduleScanner) {
            Analysis a = moduleScanner.getAnalysis();
            name = (a.name().equals("") ? moduleScanner.getModuleClass().getSimpleName() : a.name());
            description = a.description();
            version = a.version().equals("") ? "unknown" : a.version();
            className = moduleScanner.getModuleClass().getName();

        }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public String getVersion() {
        return version;
    }

    public String getClassName() {
        return className;
    }

    /**
     * is the chunk we've created valid?  Invalid chunk would contain null data, or various other
     * factors that eliminate a chunk from being outputted.
     *
     * @return true if it's valid, false if not.
     */
    @Override
    public boolean isValid() {
        if (name == null || name.equals("")) return false;
        if (description == null || description.equals("")) return false;
        if (version == null) return false;
        if (className == null || className.equals("")) return false;
        return true;
    }
}
