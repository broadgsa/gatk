package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.GATKErrorReport;
import org.broadinstitute.sting.utils.TextFormattingUtils;
import org.broadinstitute.sting.utils.cmdLine.*;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.io.PrintStream;
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
     *
     * @return A list of Strings that contain pleasant info about the GATK.
     */
    @Override
    protected ApplicationDetails getApplicationDetails() {
        return new ApplicationDetails(createApplicationHeader(),
                ApplicationDetails.createDefaultRunningInstructions(getClass()),
                getAdditionalHelp());
    }

    /**
     * generate an error log, given the stream to write to and the exception that was thrown
     *
     * @param stream the output stream
     * @param e      the exception
     */
    @Override
    public void generateErrorLog(PrintStream stream, Exception e) {
        GATKErrorReport report = new GATKErrorReport(e, this.argCollection);
        report.reportToStream(stream);
    }

    @Override
    protected String getAnalysisName() {
        return analysisName;
    }

    @Override
    protected GATKArgumentCollection getArgumentCollection() {
        return argCollection;
    }

    /**
     * Required main method implementation.
     */
    public static void main(String[] argv) {
        try {
            CommandLineGATK instance = new CommandLineGATK();
            start(instance, argv);
            System.exit(CommandLineProgram.result); // todo -- this is a painful hack
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }

    /**
     * Creates the a short blurb about the GATK, copyright info, and where to get documentation.
     *
     * @return The application header.
     */
    public static List<String> createApplicationHeader() {
        List<String> header = new ArrayList<String>();
        header.add("The Genome Analysis Toolkit (GATK)");
        header.add("Copyright (c) 2009 The Broad Institute");
        header.add("Please view our documentation at http://www.broadinstitute.org/gsa/wiki");
        header.add("For support, email gsahelp@broadinstitute.org");
        header.add("");
        return header;
    }

    /**
     * Retrieves additional information about GATK walkers.
     * TODO: This functionality is very similar to that employed by the HelpFormatter.  Generalize
     * the code in HelpFormatter and supply it as a helper to this method.
     *
     * @return A string summarizing the walkers available in this distribution.
     */
    private String getAdditionalHelp() {

        //HelpFormatter.LINE_WIDTH;

        final int PACKAGE_INDENT = 1;
        final int WALKER_INDENT = 3;
        final String FIELD_SEPARATOR = " | ";

        // Construct a help string to output available walkers.
        StringBuilder additionalHelp = new StringBuilder();
        Formatter formatter = new Formatter(additionalHelp);

        formatter.format("Available analyses:%n%n");

        // Get the list of walker names from the walker manager.
        WalkerManager walkerManager = GATKEngine.getWalkerManager();
        Map<String,Collection<Class<? extends Walker>>> walkers = walkerManager.getWalkerNamesByPackage();

        int longestPackageName = 0;
        int longestWalkerName = 0;
        for(Map.Entry<String,Collection<Class<? extends Walker>>> walkersByPackage: walkers.entrySet()) {
            longestPackageName = Math.max(longestPackageName,walkerManager.getPackageDisplayName(walkersByPackage.getKey()).length());
            for(Class<? extends Walker> walkerType: walkersByPackage.getValue())
                longestWalkerName = Math.max(longestWalkerName,walkerManager.getName(walkerType).length());
        }

        final int headerWidth = Math.max(longestPackageName+PACKAGE_INDENT,longestWalkerName+WALKER_INDENT);

        // Sort the list of walker names.
        walkers = new TreeMap<String,Collection<Class<? extends Walker>>>(walkers);

        for(String packageName: walkers.keySet()) {
            String packageDisplayName = walkerManager.getPackageDisplayName(packageName);
            String packageHelpText = walkerManager.getPackageHelpText(packageName);
            printDescriptorLine(formatter,PACKAGE_INDENT,packageDisplayName,headerWidth,FIELD_SEPARATOR,packageHelpText, TextFormattingUtils.DEFAULT_LINE_WIDTH);
            
            for(Class<? extends Walker> walkerType: walkers.get(packageName)) {
                String walkerName = walkerManager.getName(walkerType);
                String walkerHelpText = walkerManager.getWalkerHelpText(walkerType);
                printDescriptorLine(formatter,WALKER_INDENT,walkerName,headerWidth,FIELD_SEPARATOR,walkerHelpText, TextFormattingUtils.DEFAULT_LINE_WIDTH);
            }

            // Print a blank line between sets of walkers.
            printDescriptorLine(formatter,0,"",headerWidth,FIELD_SEPARATOR,"", TextFormattingUtils.DEFAULT_LINE_WIDTH);
        }

        return additionalHelp.toString();
    }

    private void printDescriptorLine(Formatter formatter,
                                     int headerIndentWidth,
                                     String header,
                                     int headerWidth,
                                     String fieldSeparator,
                                     String description,
                                     int lineWidth) {
        final int headerPaddingWidth = headerWidth - header.length() - headerIndentWidth;
        final int descriptionWidth = lineWidth - fieldSeparator.length() - headerWidth;
        List<String> wordWrappedText = TextFormattingUtils.wordWrap(description,descriptionWidth);

        String headerIndentFormatString  = headerIndentWidth  > 0 ? "%" + headerIndentWidth  + "s" : "%s";
        String headerPaddingFormatString = headerPaddingWidth > 0 ? "%" + headerPaddingWidth + "s" : "%s";
        String headerWidthFormatString   = headerWidth        > 0 ? "%" + headerWidth        + "s" : "%s";

        // Output description line.
        formatter.format(headerIndentFormatString + "%s" + headerPaddingFormatString + "%s%s%n",
                "", header, "", fieldSeparator, wordWrappedText.size()>0?wordWrappedText.get(0):"");
        for(int i = 1; i < wordWrappedText.size(); i++)
            formatter.format(headerWidthFormatString + "%s%s%n", "", fieldSeparator, wordWrappedText.get(i));
    }
}
