package org.broadinstitute.sting.playground.somaticcoverage;

import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *
 *
 * a executable command line for the Somatic Coverage Walker.
 */
public class SomaticCoverageTool extends CommandLineExecutable {    
    // our genome analysis engine
    private GenomeAnalysisEngine GATKEngine = null;

    // the two sam/bam files, one for cancer, one for normal
    @Argument(fullName = "bam_file", shortName = "I", doc = "The bam files, one for the tumor one for the normal", required = true)
    public List<File> samFiles = new ArrayList<File>();

    /** Required main method implementation. */
    public static void main( String[] argv ) {
        try {
            SomaticCoverageTool instance = new SomaticCoverageTool();
            start(instance, argv);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }

    @Override
    protected GenomeAnalysisEngine getGATKEngine() {
        if( GATKEngine == null )
            GATKEngine = new GenomeAnalysisEngine( null );
        return GATKEngine;
    }


    /**
     * a required method, returns the analysis name.  This is usually the walker
     * name with 'Walker' stripped off.
     *
     * @return the name of the analysis we're running
     */
    protected String getAnalysisName() {
        return "SomaticCoverage";
    }

    /** override any arguments we see fit. */
    protected GATKArgumentCollection getArgumentCollection() {
        GATKArgumentCollection argCollection = GATKArgumentCollection.unmarshal(getClass().getClassLoader().getResourceAsStream("SomaticCoverage.xml"));
        argCollection.samFiles = samFiles;
        return argCollection;
    }
}
