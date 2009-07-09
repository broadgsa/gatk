package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.ArgumentException;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.xReadLines;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.io.File;
import java.io.FileNotFoundException;
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
 *         <p/>
 *         Generate a executable class for the SomaticCoverageWalker
 */
public abstract class CommandLineExecutable extends CommandLineProgram {

    // our argument collection, the collection of command line args we accept
    protected GATKArgumentCollection argCollection = new GATKArgumentCollection();

    /** the type of analysis to run - switch to the type of analysis */
    private String analysisName = "SomaticCoverageWalker";

    // our genome analysis engine
    private GenomeAnalysisEngine GATKEngine = new GenomeAnalysisEngine();

    // get the analysis name
    protected abstract String getAnalysisName();

    // override select arguments
    protected abstract void overrideArguments();

    /**
     * this is the function that the inheriting class can expect to have called
     * when the command line system has initialized.
     *
     * @return the return code to exit the program with
     */
    protected int execute() {
        Walker<?,?> mWalker = GATKEngine.getWalkerByName(analysisName);

        // load the arguments into the walkers
        loadArgumentsIntoObject(argCollection);
        loadArgumentsIntoObject(mWalker);

        // process any arguments that need a second pass
        processArguments(argCollection);

        // set the analysis name in the argument collection
        this.argCollection.analysisName = this.analysisName;
        try {
            GATKEngine.execute(argCollection, mWalker);
        }
        catch (ArgumentException ex) {
            // Rethrow argument exceptions.  Let the command-line argument do what it's designed to do.
            throw ex;
        }
        catch (StingException exp) {
            System.err.println("Caught StingException. It's message is " + exp.getMessage());
            exp.printStackTrace();
            return -1;
        }
        return 0;
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
     * GATK provides the walker as an argument source.  As a side-effect, initializes the walker variable.
     *
     * @return List of walkers to load dynamically.
     */
    @Override
    protected Class[] getArgumentSources() {
        // No walker info?  No plugins.
        if (analysisName == null) return new Class[] {};
        return new Class[] { GATKEngine.getWalkerByName(analysisName).getClass() };
    }


    /**
     * Preprocess the arguments before submitting them to the GATK engine.
     *
     * @param argCollection Collection of arguments to preprocess.
     */
    private void processArguments( GATKArgumentCollection argCollection ) {
        argCollection.samFiles = unpackReads(argCollection.samFiles);
    }

    /**
     * Unpack the files to be processed, given a list of files.  That list of files can
     * itself contain lists of other files to be read.
     *
     * @param inputFiles
     *
     * @return
     */
    private List<File> unpackReads( List<File> inputFiles ) {
        List<File> unpackedReads = new ArrayList<File>();
        for (File inputFile : inputFiles) {
            if (inputFile.getName().endsWith(".list")) {
                try {
                    for (String fileName : new xReadLines(inputFile))
                        unpackedReads.add(new File(fileName));
                }
                catch (FileNotFoundException ex) {
                    throw new StingException("Unable to find file while unpacking reads", ex);
                }
            } else
                unpackedReads.add(inputFile);
        }
        return unpackedReads;
    }

    @Override
    protected String getArgumentSourceName( Class argumentSource ) {
        return WalkerManager.getWalkerName((Class<Walker>) argumentSource);
    }

    public GATKArgumentCollection getArgCollection() {
        return argCollection;
    }

    public void setArgCollection( GATKArgumentCollection argCollection ) {
        this.argCollection = argCollection;
    }


}
