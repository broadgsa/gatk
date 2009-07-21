package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.cmdLine.*;

import java.util.*;

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
    protected ApplicationDetails getApplicationDetails() {
        return new ApplicationDetails( createApplicationHeader(),
                                       ApplicationDetails.createDefaultRunningInstructions(getClass()),
                                       getAdditionalHelp() );
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

    /**
     * Creates the a short blurb about the GATK, copyright info, and where to get documentation.
     * @return The application header.
     */
    private List<String> createApplicationHeader() {
        List<String> header = new ArrayList<String>();
        header.add("The Genome Analysis Toolkit (GATK)");
        header.add("Copyright (c) 2009 The Broad Institute");
        header.add("Please view our documentation at http://www.broadinstitute.org/gsa/wiki");
        header.add("For support, email gsadevelopers@broadinstitute.org");
        header.add("");
        return header;
    }

    /**
     * Retrieves additional information about GATK walkers.
     * TODO: This functionality is very similar to that employed by the HelpFormatter.  Generalize
     *       the code in HelpFormatter and supply it as a helper to this method.
     * @return A string summarizing the walkers available in this distribution.
     */
    private String getAdditionalHelp() {
        // Get the list of walker names from the walker manager.
        Set<String> walkerNames = GATKEngine.getWalkerNames();

        // Sort the list of walker names.
        walkerNames = new TreeSet<String>( walkerNames );

        // Construct a help string to output available walkers.         
        StringBuilder additionalHelp = new StringBuilder();
        Formatter formatter = new Formatter( additionalHelp );

        formatter.format( "Available analyses:%n" );

        // Compute the max size of any walker name
        int maxNameLength = 0;
        for( String walkerName: walkerNames ) {
            if( maxNameLength < walkerName.length() )
                maxNameLength = walkerName.length();
        }
        
        final int fieldWidth = maxNameLength + HelpFormatter.FIELD_SEPARATION_WIDTH;
        final int walkersPerLine = Math.min(HelpFormatter.LINE_WIDTH / fieldWidth, 4 );
        final int columnSpacing = (HelpFormatter.LINE_WIDTH - (fieldWidth * walkersPerLine)) / walkersPerLine;

        int currentWalkerName = 0;
        for( String walkerName: walkerNames ) {
            formatter.format( "%-" + HelpFormatter.FIELD_SEPARATION_WIDTH + "s" +
                              "%-" + fieldWidth + "s" +
                              "%-" + columnSpacing + "s", "", walkerName, "" );
            if( ++currentWalkerName % walkersPerLine == 0 )
                formatter.format("%n");
        }

        return additionalHelp.toString();                
    }
}
