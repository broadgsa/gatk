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

import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

import net.sf.picard.filter.SamRecordFilter;

/**
 * @author aaron
 */
public abstract class CommandLineExecutable extends CommandLineProgram {
    public Object GATKResult = null;

    /**
     * The actual engine which performs the analysis.
     */
    protected GenomeAnalysisEngine GATKEngine = new GenomeAnalysisEngine();

    // get the analysis name
    protected abstract String getAnalysisName();

    /**
     * Gets the GATK argument bundle.
     * @return A structure consisting of whatever arguments should be used to initialize the GATK engine.
     */
    protected abstract GATKArgumentCollection getArgumentCollection();

    // override select arguments
    protected void overrideArguments() { }

    /**
     * this is the function that the inheriting class can expect to have called
     * when the command line system has initialized.
     *
     * @return the return code to exit the program with
     */
    protected int execute() throws Exception {
        Walker<?,?> mWalker = GATKEngine.getWalkerByName(getAnalysisName());
        Collection<SamRecordFilter> filters = GATKEngine.createFiltersForWalker(getArgumentCollection(),mWalker);

        // load the arguments into the walker / filters.
        loadArgumentsIntoObject(mWalker);
        for (SamRecordFilter filter: filters)
            loadArgumentsIntoObject(filter);

        // set the analysis name in the argument collection
        try {
            GATKResult = GATKEngine.execute(getArgumentCollection(), mWalker, filters);
            generateGATKRunReport(mWalker);
        } catch ( Exception e ) {
            generateGATKRunReport(mWalker, e);
            throw e;
        }

        // always return 0
        return 0;
    }

    /**
     * generate an error log
     * @param e the exception, can be null if no exception occurred
     */
    private void generateGATKRunReport(Walker<?,?> mWalker, Exception e) {
        if ( getArgumentCollection().phoneHomeType != GATKRunReport.PhoneHomeOption.NO_ET ) {
            GATKRunReport report = new GATKRunReport(mWalker, e, GATKEngine);
            if ( getArgumentCollection().phoneHomeType == GATKRunReport.PhoneHomeOption.STDOUT )
                report.postReport(System.out);
            else
                report.postReport();
        }
    }

    private void generateGATKRunReport(Walker<?,?> mWalker) {
        generateGATKRunReport(mWalker, null);
    }

    /**
     * Subclasses of CommandLinePrograms can provide their own types of command-line arguments.
     * @return A collection of type descriptors generating implementation-dependent placeholders.
     */
    protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
        return GATKEngine.getArgumentTypeDescriptors();
    }

    /**
     * GATK can add arguments dynamically based on analysis type.
     *
     * @return true
     */
    @Override
    protected boolean canAddArgumentsDynamically() {
        return true;
    }

    /**
     * GATK provides the walker as an argument source.
     * @return List of walkers to load dynamically.
     */
    @Override
    protected Class[] getArgumentSources() {
        // No walker info?  No plugins.
        if (getAnalysisName() == null) return new Class[] {};

        Collection<Class> argumentSources = new ArrayList<Class>();

        Walker walker = GATKEngine.getWalkerByName(getAnalysisName());
        argumentSources.add(walker.getClass());

        Collection<SamRecordFilter> filters = GATKEngine.createFiltersForWalker(getArgumentCollection(),walker);
        for(SamRecordFilter filter: filters)
            argumentSources.add(filter.getClass());

        Class[] argumentSourcesAsArray = new Class[argumentSources.size()];
        return argumentSources.toArray(argumentSourcesAsArray);
    }

    @Override
    protected String getArgumentSourceName( Class argumentSource ) {
        return GATKEngine.getWalkerName((Class<Walker>)argumentSource);
    }

    /**
     * Supply command-line argument tags to the GATK engine.
     * @param key Key to use, created by the command-line argument system.
     * @param tags List of freeform tags.
     */
    @Override
    protected void addTags(Object key, List<String> tags) {
        GATKEngine.addTags(key,tags);
    }
}
