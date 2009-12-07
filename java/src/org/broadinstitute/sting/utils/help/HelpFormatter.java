package org.broadinstitute.sting.utils.help;

import org.broadinstitute.sting.utils.cmdLine.ArgumentDefinition;
import org.broadinstitute.sting.utils.cmdLine.ArgumentDefinitionGroup;
import org.broadinstitute.sting.utils.cmdLine.ArgumentDefinitions;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.TextFormattingUtils;
import org.apache.log4j.Logger;

import java.util.Formatter;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Comparator;
import java.util.Collection;
import java.util.Collections;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
/**
 * User: hanna
 * Date: May 6, 2009
 * Time: 10:16:43 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Print out help for Sting command-line applications.
 */

public class HelpFormatter {
    /** our log, which we want to capture anything from org.broadinstitute.sting */
    private static Logger logger = Logger.getRootLogger();    

    public static final int FIELD_SEPARATION_WIDTH = 3;

    /**
     * Prints the help, given a collection of argument definitions.
     * @param argumentDefinitions Argument definitions for which help should be printed.
     */
    public void printHelp( ApplicationDetails applicationDetails, ArgumentDefinitions argumentDefinitions ) {
        List<ArgumentDefinitionGroup> argumentGroups = prepareArgumentGroups( argumentDefinitions );

        String synopsis = getSynopsis(applicationDetails.runningInstructions,argumentGroups);
        String additionalDetails = applicationDetails.additionalHelp != null ? applicationDetails.additionalHelp : "";
        String detailedDescription = getDetailed(argumentGroups);

        System.out.printf("%s%n%s%n%s%n",synopsis,detailedDescription,additionalDetails );
    }

    /**
     * Gets the synopsis: the actual command to run.
     * @param runningInstructions Instructions on how to run hte application.
     * @param argumentGroups Program arguments sorted in order of definition group displays.
     * @return A synopsis line.
     */
    private String getSynopsis( String runningInstructions,
                                List<ArgumentDefinitionGroup> argumentGroups ) {
        // Build out the synopsis all as one long line.        
        StringBuilder lineBuilder = new StringBuilder();
        Formatter lineFormatter = new Formatter( lineBuilder );

        lineFormatter.format("java %s", runningInstructions);

        for( ArgumentDefinitionGroup argumentGroup: argumentGroups ) {
            for( ArgumentDefinition argumentDefinition: argumentGroup.argumentDefinitions ) {
                lineFormatter.format(" ");
                if( !argumentDefinition.required ) lineFormatter.format("[");
                if( argumentDefinition.shortName != null )
                    lineFormatter.format("-%s", argumentDefinition.shortName);
                else
                    lineFormatter.format("--%s", argumentDefinition.fullName);
                if( !argumentDefinition.source.isFlag() )
                    lineFormatter.format(" <%s>", argumentDefinition.fullName);
                if( !argumentDefinition.required ) lineFormatter.format("]");
            }
        }

        // Word wrap the synopsis.
        List<String> wrappedSynopsis = TextFormattingUtils.wordWrap( lineBuilder.toString(), TextFormattingUtils.DEFAULT_LINE_WIDTH );

        String header = "usage: ";
        int headerLength = header.length();

        StringBuilder synopsisBuilder = new StringBuilder();
        Formatter synopsisFormatter = new Formatter(synopsisBuilder);
        for( String synopsisLine: wrappedSynopsis ) {
            synopsisFormatter.format("%" + headerLength + "s%s%n", header, synopsisLine);
            header = "";
        }

        return synopsisBuilder.toString();
    }

    /**
     * Gets detailed output about each argument type.
     * @param argumentGroups Collection of program arguments sorted according to how they should be shown. 
     * @return Detailed text about all arguments.
     */
    private String getDetailed( List<ArgumentDefinitionGroup> argumentGroups ) {
        StringBuilder builder = new StringBuilder();

        for( ArgumentDefinitionGroup argumentGroup: argumentGroups ) {
            if( argumentGroup.groupName != null && argumentGroup.argumentDefinitions.size() != 0 )
                builder.append( String.format("%nArguments for %s:%n", argumentGroup.groupName ) );
            builder.append( getDetailForGroup( argumentGroup.argumentDefinitions ) );
        }

        return builder.toString();
    }

    /**
     * Gets a detailed description for a given argument group.
     * @param argumentDefinitions The argument definitions contained withina group.
     * @return A string giving detailed info about the contents of this group.
     */
    private String getDetailForGroup( List<ArgumentDefinition> argumentDefinitions ) {
        StringBuilder builder = new StringBuilder();
        Formatter formatter = new Formatter( builder );

        // Try to fit the entire argument definition across the screen, but impose an arbitrary cap of 3/4 *
        // LINE_WIDTH in case the length of the arguments gets out of control.
        int argWidth = Math.min( findLongestArgumentCallingInfo(argumentDefinitions), (TextFormattingUtils.DEFAULT_LINE_WIDTH*3)/4 - FIELD_SEPARATION_WIDTH );
        int docWidth = TextFormattingUtils.DEFAULT_LINE_WIDTH - argWidth - FIELD_SEPARATION_WIDTH;

        for( ArgumentDefinition argumentDefinition: argumentDefinitions ) {
            Iterator<String> wordWrappedArgs = TextFormattingUtils.wordWrap( getArgumentCallingInfo(argumentDefinition), argWidth ).iterator();
            Iterator<String> wordWrappedDoc  = TextFormattingUtils.wordWrap( argumentDefinition.doc, docWidth ).iterator();

            while( wordWrappedArgs.hasNext() || wordWrappedDoc.hasNext() ) {
                String arg = wordWrappedArgs.hasNext() ? wordWrappedArgs.next() : "";
                String doc = wordWrappedDoc.hasNext() ? wordWrappedDoc.next() : "";

                String formatString = "%-" + argWidth + "s%" + FIELD_SEPARATION_WIDTH + "s%s%n";
                formatter.format( formatString, arg, "", doc );
            }
        }

        return builder.toString();
    }

    /**
     * Gets a string indicating how this argument should be passed to the application.
     * @param argumentDefinition Argument definition for which help should be printed.
     * @return Calling information for this argument.
     */
    private String getArgumentCallingInfo( ArgumentDefinition argumentDefinition ) {
        StringBuilder builder = new StringBuilder();
        Formatter formatter = new Formatter( builder );

        formatter.format(" ");
        if( argumentDefinition.shortName != null )
            formatter.format("-%s,", argumentDefinition.shortName);
        formatter.format("--%s", argumentDefinition.fullName);
        if( !argumentDefinition.source.isFlag() )
            formatter.format(" <%s>", argumentDefinition.fullName);

        return builder.toString();
    }

    /**
     * Crude implementation which finds the longest argument portion
     * given a set of arguments.
     * @param argumentDefinitions argument definitions to inspect.
     * @return longest argument length.
     */
    private int findLongestArgumentCallingInfo( Collection<ArgumentDefinition> argumentDefinitions ) {
        int longest = 0;
        for( ArgumentDefinition argumentDefinition: argumentDefinitions ) {
            String argumentText = getArgumentCallingInfo( argumentDefinition );
            if( longest < argumentText.length() )
                longest = argumentText.length();
        }
        return longest;
    }

    /**
     * Extract the argument definition groups from the argument definitions and arrange them appropriately.
     * For help, we want the arguments sorted as they are declared in the class.  However, required arguments
     * should appear before optional arguments.
     * @param argumentDefinitions Argument definitions from which to extract argument groups.
     * @return A list of argument groups sorted in display order.
     */
    private List<ArgumentDefinitionGroup> prepareArgumentGroups( ArgumentDefinitions argumentDefinitions ) {
        // Sort the list of argument definitions according to how they should be shown.
        // Put the sorted results into a new cloned data structure.
        Comparator definitionComparator = new Comparator<ArgumentDefinition>() {
            public int compare( ArgumentDefinition lhs, ArgumentDefinition rhs ) {
                if( lhs.required && rhs.required ) return 0;
                if( lhs.required ) return -1;
                if( rhs.required ) return 1;
                return 0;
            }
        };

        List<ArgumentDefinitionGroup> argumentGroups = new ArrayList();        
        for( ArgumentDefinitionGroup argumentGroup: argumentDefinitions.getArgumentDefinitionGroups() ) {
            List<ArgumentDefinition> sortedDefinitions = new ArrayList( argumentGroup.argumentDefinitions );
            Collections.sort( sortedDefinitions, definitionComparator );
            argumentGroups.add( new ArgumentDefinitionGroup(argumentGroup.groupName,sortedDefinitions) );
        }

        // Sort the argument groups themselves with main arguments first, followed by plugins sorted in name order.
        Comparator groupComparator = new Comparator<ArgumentDefinitionGroup>() {
            public int compare( ArgumentDefinitionGroup lhs, ArgumentDefinitionGroup rhs ) {
                if( lhs.groupName == null && rhs.groupName == null ) return 0;
                if( lhs.groupName == null ) return -1;
                if( rhs.groupName == null ) return 1;
                return lhs.groupName.compareTo(rhs.groupName);
            }
        };
        Collections.sort( argumentGroups, groupComparator );


        return argumentGroups;
    }

    /**
     * generateHeaderInformation
     * <p/>
     * <p/>
     * Generate a standard header for the logger
     *
     * @param applicationDetails details of the application to run.
     * @param args the command line arguments passed in
     */
    public static void generateHeaderInformation(ApplicationDetails applicationDetails, String[] args) {

        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        java.util.Date date = new java.util.Date();

        logger.info("-------------------------------------------------------");
        for (String headerLine : applicationDetails.applicationHeader)
            logger.info(headerLine);
        String output = "";
        for (String str : args) {
            output = output + str + " ";
        }
        logger.info("Program Args: " + output);
        logger.info("Date/Time: " + dateFormat.format(date));
        logger.info("-------------------------------------------------------");
    }
}
