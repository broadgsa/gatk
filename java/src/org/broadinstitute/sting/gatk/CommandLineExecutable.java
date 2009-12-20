package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.ArgumentTypeDescriptor;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.xReadLines;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileReaderArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.GenotypeWriterArgumentTypeDescriptor;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import net.sf.picard.filter.SamRecordFilter;

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
    protected int execute() {
        GATKResult = executeGATK();
        return 0;
    }

    protected Object executeGATK() {
        Walker<?,?> mWalker = GATKEngine.getWalkerByName(getAnalysisName());
        Collection<SamRecordFilter> filters = GATKEngine.createFiltersForWalker(getArgumentCollection(),mWalker);

        // load the arguments into the walker / filters.
        loadArgumentsIntoObject(mWalker);
        for(SamRecordFilter filter: filters)
            loadArgumentsIntoObject(filter);

        // process any arguments that need a second pass
        GATKArgumentCollection arguments = getArgumentCollection();
        processArguments(arguments);

        // set the analysis name in the argument collection
        return GATKEngine.execute(arguments, mWalker, filters);
    }

    /**
     * Subclasses of CommandLinePrograms can provide their own types of command-line arguments.
     * @return A collection of type descriptors generating implementation-dependent placeholders.
     */
    protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
        return Arrays.asList( new GenotypeWriterArgumentTypeDescriptor(GATKEngine),
                              new SAMFileReaderArgumentTypeDescriptor(GATKEngine),
                              new SAMFileWriterArgumentTypeDescriptor(GATKEngine),
                              new OutputStreamArgumentTypeDescriptor(GATKEngine) );
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
     * Preprocess the arguments before submitting them to the GATK engine.
     *
     * @param argCollection Collection of arguments to preprocess.
     */
    private void processArguments( GATKArgumentCollection argCollection ) {
        argCollection.samFiles = unpackList( argCollection.samFiles );
    }

    /**
     * Unpack the files to be processed, given a list of files.  That list of files can
     * itself contain lists of other files to be read.
     *
     * @param inputFiles
     *
     * @return
     */
    public static List<File> unpackList( List<File> inputFiles ) {
        List<File> unpackedReads = new ArrayList<File>();
        for( File inputFile: inputFiles ) {
            if (inputFile.getName().endsWith(".list") ) {
                try {
                    for( String fileName : new xReadLines(inputFile) )
                        unpackedReads.add( new File(fileName) );
                }
                catch( FileNotFoundException ex ) {
                    throw new StingException("Unable to find file while unpacking reads", ex);
                }
            }
            else
                unpackedReads.add( inputFile );
        }
        return unpackedReads;
    }
}
