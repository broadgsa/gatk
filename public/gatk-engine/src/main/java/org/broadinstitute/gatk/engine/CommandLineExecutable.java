/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.gatk.engine.io.stubs.VCFWriterArgumentTypeDescriptor;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.text.ListFileUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

/**
 * @author aaron
 */
public abstract class CommandLineExecutable extends CommandLineProgram {
    /**
     * The actual engine which performs the analysis.
     */
    protected GenomeAnalysisEngine engine = new GenomeAnalysisEngine();

    // get the analysis name
    public abstract String getAnalysisName();

    /**
     * Gets the GATK argument bundle.
     * @return A structure consisting of whatever arguments should be used to initialize the GATK engine.
     */
    protected abstract GATKArgumentCollection getArgumentCollection();

    /**
     * A list of all the arguments initially used as sources.
     */
    private final Collection<Object> argumentSources = new ArrayList<Object>();

    protected static Logger logger = Logger.getLogger(CommandLineExecutable.class);

    /**
     * this is the function that the inheriting class can expect to have called
     * when the command line system has initialized.
     *
     * @return the return code to exit the program with
     */
    protected int execute() throws Exception {
        engine.setParser(parser);
        argumentSources.add(this);

        Walker<?,?> walker = engine.getWalkerByName(getAnalysisName());

        try {
            engine.setArguments(getArgumentCollection());

            // File lists can require a bit of additional expansion.  Set these explicitly by the engine. 
            final Collection<SAMReaderID> bamFileList=ListFileUtils.unpackBAMFileList(getArgumentCollection().samFiles,parser);
            engine.setSAMFileIDs(bamFileList);
            if(getArgumentCollection().showFullBamList){
                logger.info(String.format("Adding the following input SAM Files: %s",bamFileList.toString()));
            }

            engine.setWalker(walker);
            walker.setToolkit(engine);

            Collection<ReadFilter> filters = engine.createFilters();
            engine.setFilters(filters);

            // load the arguments into the walker / filters.
            // TODO: The fact that this extra load call exists here when all the parsing happens at the engine
            // TODO: level indicates that we're doing something wrong.  Turn this around so that the GATK can drive
            // TODO: argument processing.
            loadArgumentsIntoObject(walker);
            argumentSources.add(walker);

            Collection<RMDTriplet> rodBindings = ListFileUtils.unpackRODBindings(parser.getRodBindings(), parser);
            engine.setReferenceMetaDataFiles(rodBindings);

            for (ReadFilter filter: filters) {
                loadArgumentsIntoObject(filter);
                argumentSources.add(filter);
            }

            engine.execute();
        } catch ( Exception e ) {
            throw e;
        }

        // always return 0
        return 0;
    }

    /**
     * Subclasses of CommandLinePrograms can provide their own types of command-line arguments.
     * @return A collection of type descriptors generating implementation-dependent placeholders.
     */
    protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
        return Arrays.asList( new VCFWriterArgumentTypeDescriptor(engine,System.out,argumentSources),
                              new SAMFileWriterArgumentTypeDescriptor(engine,System.out),
                              new OutputStreamArgumentTypeDescriptor(engine,System.out) );
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

        Walker walker = engine.getWalkerByName(getAnalysisName());
        engine.setArguments(getArgumentCollection());
        engine.setWalker(walker);
        walker.setToolkit(engine);
        argumentSources.add(walker.getClass());

        Collection<ReadFilter> filters = engine.createFilters();
        for(ReadFilter filter: filters)
            argumentSources.add(filter.getClass());

        Class[] argumentSourcesAsArray = new Class[argumentSources.size()];
        return argumentSources.toArray(argumentSourcesAsArray);
    }

    @Override
    protected String getArgumentSourceName( Class argumentSource ) {
        return engine.getWalkerName((Class<Walker>)argumentSource);
    }

}
