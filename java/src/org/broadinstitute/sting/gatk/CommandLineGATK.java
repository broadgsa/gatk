package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.cmdLine.*;

import java.util.List;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;

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
public class CommandLineGATK extends CommandLineExecutable {
    // our genome analysis engine
    private GenomeAnalysisEngine GATKEngine = null;

    @Argument(fullName = "analysis_type", shortName = "T", doc = "Type of analysis to run")
    private String analysisName = null;

    // our argument collection, the collection of command line args we accept
    @ArgumentCollection
    private GATKArgumentCollection argCollection = new GATKArgumentCollection();    

    /**
     * Get pleasing info about the GATK.
     * @return A list of Strings that contain pleasant info about the GATK.
     */
    @Override
    protected List<String> getApplicationHeader() {
        List<String> header = new ArrayList<String>();
        header.add("The Genome Analysis Toolkit (GATK)");
        header.add("Copyright (c) 2009 The Broad Institute");
        header.add("Please view our documentation at http://www.broadinstitute.org/gsa/wiki");
        header.add("For support, email gsadevelopers@broadinstitute.org");
        header.add("");
        return header;
    }

    /**
     * Lazy load the GATK engine.  This current CANNOT happen until after the command-line arguments are populated.
     * TODO: Make this chain of events more explicit.  Perhaps an abstract initialize method after clp arguments are parsed?
     * @return The GATK engine that will power the requested traversal.
     */
    @Override
    protected GenomeAnalysisEngine getGATKEngine() {
        if( GATKEngine == null )
            GATKEngine = new GenomeAnalysisEngine();
        return GATKEngine;
    }


    @Override
    protected String getAnalysisName() {
        return analysisName;
    }

    @Override
    protected GATKArgumentCollection getArgumentCollection() {
        return argCollection;
    }

    /** Required main method implementation. */
    public static void main(String[] argv) {
        try {
            CommandLineGATK instance = new CommandLineGATK();
            start(instance, argv);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }
}
