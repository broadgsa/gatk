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
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.help.ApplicationDetails;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.util.*;

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

    @Override
    public String getAnalysisName() {
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
        } catch (UserException e) {
            exitSystemWithUserError(e);
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
        header.add(String.format("The Genome Analysis Toolkit (GATK) v%s, Compiled %s",getVersionNumber(), getBuildTime()));
        header.add("Copyright (c) 2009 The Broad Institute");
        header.add("Please view our documentation at http://www.broadinstitute.org/gsa/wiki");
        header.add("For support, please view our support site at http://getsatisfaction.com/gsa");
        return header;
    }

    public static String getVersionNumber() {
        ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        return headerInfo.containsKey("org.broadinstitute.sting.gatk.version") ? headerInfo.getString("org.broadinstitute.sting.gatk.version") : "<unknown>";
    }

    public static String getBuildTime() {
        ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        return headerInfo.containsKey("build.timestamp") ? headerInfo.getString("build.timestamp") : "<unknown>";
    }

    /**
     * Retrieves additional information about GATK walkers.
     * the code in HelpFormatter and supply it as a helper to this method.
     *
     * @return A string summarizing the walkers available in this distribution.
     */
    private String getAdditionalHelp() {
        String additionalHelp;

        // If no analysis name is present, fill in extra help on the walkers.
        WalkerManager walkerManager = engine.getWalkerManager();
        String analysisName = getAnalysisName();
        if(analysisName != null && walkerManager.exists(getAnalysisName()))
            additionalHelp = getWalkerHelp(walkerManager.getWalkerClassByName(getAnalysisName()));
        else
            additionalHelp = getAllWalkerHelp();

        return additionalHelp;
    }

    private static final int PACKAGE_INDENT = 1;
    private static final int WALKER_INDENT = 3;
    private static final String FIELD_SEPARATOR = "  ";

    private String getWalkerHelp(Class<? extends Walker> walkerType) {
        // Construct a help string to output details on this walker.
        StringBuilder additionalHelp = new StringBuilder();
        Formatter formatter = new Formatter(additionalHelp);

        formatter.format("Description:%n");

        WalkerManager walkerManager = engine.getWalkerManager();
        String walkerHelpText = walkerManager.getWalkerDescriptionText(walkerType);

        printDescriptorLine(formatter,WALKER_INDENT,"",WALKER_INDENT,FIELD_SEPARATOR,walkerHelpText,TextFormattingUtils.DEFAULT_LINE_WIDTH);

        return additionalHelp.toString();
    }

    /**
     * Load in additional help information about all available walkers.
     * @return A string representation of the additional help.
     */
    private String getAllWalkerHelp() {
        // Construct a help string to output available walkers.
        StringBuilder additionalHelp = new StringBuilder();
        Formatter formatter = new Formatter(additionalHelp);

        formatter.format("Available analyses:%n");

        // Get the list of walker names from the walker manager.
        WalkerManager walkerManager = engine.getWalkerManager();

        // Build a list sorted by walker display name.  As this information is collected, keep track of the longest
        // package / walker name for later formatting.
        SortedSet<HelpEntry> helpText = new TreeSet<HelpEntry>(new HelpEntryComparator());
        
        int longestPackageName = 0;
        int longestWalkerName = 0;
        for(Map.Entry<String,Collection<Class<? extends Walker>>> walkersByPackage: walkerManager.getWalkerNamesByPackage(true).entrySet()) {
            // Get the display name.
            String packageName = walkersByPackage.getKey();
            String packageDisplayName = walkerManager.getPackageDisplayName(walkersByPackage.getKey());
            String packageHelpText = walkerManager.getPackageSummaryText(packageName);

            // Compute statistics about which names is longest.
            longestPackageName = Math.max(longestPackageName,packageDisplayName.length());

            SortedSet<HelpEntry> walkersInPackage = new TreeSet<HelpEntry>(new HelpEntryComparator());
            for(Class<? extends Walker> walkerType: walkersByPackage.getValue()) {
                String walkerName = walkerType.getName();
                String walkerDisplayName = walkerManager.getName(walkerType);
                String walkerHelpText = walkerManager.getWalkerSummaryText(walkerType);                

                longestWalkerName = Math.max(longestWalkerName,walkerManager.getName(walkerType).length());

                walkersInPackage.add(new HelpEntry(walkerName,walkerDisplayName,walkerHelpText));
            }

            // Dump the walkers into the sorted set.
            helpText.add(new HelpEntry(packageName,packageDisplayName,packageHelpText,Collections.unmodifiableSortedSet(walkersInPackage)));
        }

        final int headerWidth = Math.max(longestPackageName+PACKAGE_INDENT,longestWalkerName+WALKER_INDENT);


        for(HelpEntry packageHelp: helpText) {
            printDescriptorLine(formatter,PACKAGE_INDENT,packageHelp.displayName,headerWidth,FIELD_SEPARATOR,packageHelp.summary,TextFormattingUtils.DEFAULT_LINE_WIDTH);
            
            for(HelpEntry walkerHelp: packageHelp.children)
                printDescriptorLine(formatter,WALKER_INDENT,walkerHelp.displayName,headerWidth,FIELD_SEPARATOR,walkerHelp.summary,TextFormattingUtils.DEFAULT_LINE_WIDTH);

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

/**
 * Represents a given help entry; contains a display name, a summary and optionally some children.
 */
class HelpEntry {
    public final String uid;
    public final String displayName;
    public final String summary;
    public final SortedSet<HelpEntry> children;

    /**
     * Create a new help entry with the given display name, summary and children.
     * @param uid a unique identifier.  Usually, the java package.
     * @param displayName display name for this help entry.
     * @param summary summary for this help entry.
     * @param children children for this help entry.
     */
    public HelpEntry(String uid, String displayName, String summary, SortedSet<HelpEntry> children)  {
        this.uid = uid;
        this.displayName = displayName;
        this.summary = summary;
        this.children = children;
    }

    /**
     * Create a new help entry with the given display name, summary and children.
     * @param uid a unique identifier.  Usually, the java package.
     * @param displayName display name for this help entry.
     * @param summary summary for this help entry.
     */
    public HelpEntry(String uid, String displayName, String summary) {
        this(uid,displayName,summary,null);
    }

}

/**
 * Compare two help entries by display name.
 */
class HelpEntryComparator implements Comparator<HelpEntry> {
    private static TextFormattingUtils.CaseInsensitiveComparator textComparator = new TextFormattingUtils.CaseInsensitiveComparator();

    /**
     * Compares the order of lhs to rhs, not taking case into account.
     * @param lhs First object to compare.
     * @param rhs Second object to compare.
     * @return 0 if objects are identical; -1 if lhs is before rhs, 1 if rhs is before lhs.  Nulls are treated as after everything else.
     */
    public int compare(HelpEntry lhs, HelpEntry rhs) {
        if(lhs == null && rhs == null) return 0;
        if(lhs == null || lhs.displayName.equals("")) return 1;
        if(rhs == null || rhs.displayName.equals("")) return -1;
        return lhs.displayName.equals(rhs.displayName) ? textComparator.compare(lhs.uid,rhs.uid) : textComparator.compare(lhs.displayName,rhs.displayName);
    }


}