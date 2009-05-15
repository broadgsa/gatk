package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.xReadLines;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.ArgumentCollection;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.FileNotFoundException;
import java.io.File;
import java.util.List;
import java.util.ArrayList;

/**
 *
 * User: aaron
 * Date: May 8, 2009
 * Time: 10:50:58 AM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date May 8, 2009
 * <p/>
 * Class CommandLineGATK
 * <p/>
 * We run command line GATK programs using this class.  It gets the command line args, parses them, and hands the
 * gatk all the parsed out information.  Pretty much anything dealing with the underlying system should go here,
 * the gatk engine should  deal with any data related information.
 */
public class CommandLineGATK extends CommandLineProgram {

    @Argument(fullName = "analysis_type", shortName = "T", doc = "Type of analysis to run")
    public String analysisName = null;

    @ArgumentCollection // our argument collection, the collection of command line args we accept
    public GATKArgumentCollection argCollection = new GATKArgumentCollection();

    public String pluginPathName = null;

    // our genome analysis engine
    GenomeAnalysisEngine GATKEngine = null;

    // our walker manager
    private WalkerManager walkerManager = null;


    /** Required main method implementation. */
    public static void main(String[] argv) {
        try {
            CommandLineGATK instance = new CommandLineGATK();
            start(instance, argv);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }


    /**
     * this is the function that the inheriting class can expect to have called
     * when the command line system has initialized.
     *
     * @return the return code to exit the program with
     */
    protected int execute() {
        Walker<?, ?> mWalker = null;
        try {
            mWalker = walkerManager.createWalkerByName(analysisName);
        } catch (InstantiationException ex) {
            throw new RuntimeException("Unable to instantiate walker.", ex);
        }
        catch (IllegalAccessException ex) {
            throw new RuntimeException("Unable to access walker", ex);
        }
        loadArgumentsIntoObject(argCollection);
        loadArgumentsIntoObject(mWalker);

        processArguments(argCollection);

        this.argCollection.analysisName = this.analysisName;
        try {
            GATKEngine = new GenomeAnalysisEngine(argCollection, mWalker);
        } catch (StingException exp) {
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
        if (analysisName == null)
            throw new IllegalArgumentException("Must provide analysis name");

        walkerManager = new WalkerManager(pluginPathName);

        if (!walkerManager.doesWalkerExist(analysisName))
            throw new IllegalArgumentException("Invalid analysis name");

        return new Class[]{walkerManager.getWalkerClassByName(analysisName)};
    }

    @Override
    protected String getArgumentSourceName(Class argumentSource) {
        return WalkerManager.getWalkerName((Class<Walker>) argumentSource);
    }

    public GATKArgumentCollection getArgCollection() {
        return argCollection;
    }

    public void setArgCollection(GATKArgumentCollection argCollection) {
        this.argCollection = argCollection;
    }


    /**
     * Preprocess the arguments before submitting them to the GATK engine.
     * @param argCollection Collection of arguments to preprocess.
     */
    private void processArguments( GATKArgumentCollection argCollection ) {
        argCollection.samFiles = unpackReads( argCollection.samFiles );    
    }

    /**
     * Unpack the files to be processed, given a list of files.  That list of files can
     * itself contain lists of other files to be read.
     * @param inputFiles
     * @return
     */
    private List<File> unpackReads( List<File> inputFiles ) {
        List<File> unpackedReads = new ArrayList<File>();
        for( File inputFile: inputFiles ) {
            if( inputFile.getName().endsWith(".list") ) {
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
