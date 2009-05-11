package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.ArgumentCollection;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.gatk.GATKArgumentCollection;

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
            mWalker = walkerManager.createWalkerByName(argCollection.analysisName);
        } catch (InstantiationException ex) {
            throw new RuntimeException("Unable to instantiate walker.", ex);
        }
        catch (IllegalAccessException ex) {
            throw new RuntimeException("Unable to access walker", ex);
        }
        loadArgumentsIntoObject(mWalker);
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
        loadArgumentsIntoObject(this.argCollection);
        if (argCollection.analysisName == null)
            throw new IllegalArgumentException("Must provide analysis name");

        walkerManager = new WalkerManager(pluginPathName);

        if (!walkerManager.doesWalkerExist(argCollection.analysisName))
            throw new IllegalArgumentException("Invalid analysis name");

        return new Class[]{walkerManager.getWalkerClassByName(argCollection.analysisName)};
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

}
