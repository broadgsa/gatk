/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.reads.performance;

import net.sf.picard.filter.SamRecordFilter;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.io.File;
import java.util.Collections;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 25, 2011
 * Time: 10:16:54 AM
 * To change this template use File | Settings | File Templates.
 */
class GATKWalkerInvoker extends ReadProcessor {
    /**
     * Walker to run over the existing dataset.
     */
    private Walker<?,?> walker;

    public GATKWalkerInvoker(BAMProcessingPerformanceMeter performanceMeter) {
        super(performanceMeter);
    }

    @Override
    public String getTestName() { return "GATK-CountReads"; }

    public void setWalker(Walker<?,?> walker) {
        this.walker = walker;
    }

    @Override
    public void execute(File samFile, File fastaFile) {
        GenomeAnalysisEngine engine = new GenomeAnalysisEngine();

        // Establish the argument collection
        GATKArgumentCollection argCollection = new GATKArgumentCollection();
        argCollection.referenceFile = fastaFile;
        argCollection.samFiles = Collections.singletonList(samFile.getAbsolutePath());

        engine.setArguments(argCollection);
        // Bugs in the engine mean that this has to be set twice.
        engine.setSAMFileIDs(Collections.singletonList(new SAMReaderID(samFile,new Tags())));
        engine.setFilters(Collections.<SamRecordFilter>emptyList());
        engine.setReferenceMetaDataFiles(Collections.<RMDTriplet>emptyList());

        // Create the walker
        engine.setWalker(walker);

        startTest();
        engine.execute();
        stopTest();
    }

}
