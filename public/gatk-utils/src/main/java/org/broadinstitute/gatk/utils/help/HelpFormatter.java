/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.help;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.text.TextFormattingUtils;

import java.net.InetAddress;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
/**
 * Print out help for GATK command-line applications.
 */

public class HelpFormatter {
    /** our log, which we want to capture anything from org.broadinstitute.gatk */
    private static Logger logger = Logger.getLogger(HelpFormatter.class);

    public static final int FIELD_SEPARATION_WIDTH = 3;

    /**
     * Prints the help, given a collection of argument definitions.
     * @param applicationDetails Application details
     * @param argumentDefinitions Argument definitions for which help should be printed.
     */
    public void printHelp( ApplicationDetails applicationDetails, ArgumentDefinitions argumentDefinitions ) {
        List<ArgumentDefinitionGroup> argumentGroups = prepareArgumentGroups( argumentDefinitions );

        List<String> header = applicationDetails.applicationHeader;
        String barrier = createBarrier(header);

        System.out.printf("%s%n",barrier);
        for(String headerLine: header)
            System.out.printf("%s%n",headerLine);
        System.out.printf("%s%n",barrier);
        for(String attributionLine: applicationDetails.attribution)
            System.out.printf("%s%n",attributionLine);
        System.out.printf("%s%n",barrier);

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
                if(argumentDefinition.isHidden)
                    continue;
                lineFormatter.format(" ");
                if( !argumentDefinition.required ) lineFormatter.format("[");
                if( argumentDefinition.shortName != null )
                    lineFormatter.format("-%s", argumentDefinition.shortName);
                else
                    lineFormatter.format("--%s", argumentDefinition.fullName);
                if( !argumentDefinition.isFlag )
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

        for( ArgumentDefinitionGroup argumentGroup: argumentGroups )
            builder.append( getDetailForGroup( argumentGroup ) );

        return builder.toString();
    }

    /**
     * Gets a detailed description for a given argument group.
     * @param argumentDefinitionGroup The group of argument definitions to render.
     * @return A string giving detailed info about the contents of this group.
     */
    private String getDetailForGroup( ArgumentDefinitionGroup argumentDefinitionGroup ) {
        if(argumentDefinitionGroup.allHidden())
            return "";

        StringBuilder builder = new StringBuilder();
        Formatter formatter = new Formatter( builder );

        if( argumentDefinitionGroup.groupName != null && argumentDefinitionGroup.argumentDefinitions.size() != 0 )
            builder.append( String.format("%nArguments for %s:%n", argumentDefinitionGroup.groupName ) );

        List<ArgumentDefinition> argumentDefinitions = new ArrayList<ArgumentDefinition>();
        for(ArgumentDefinition argumentDefinition: argumentDefinitionGroup.argumentDefinitions) {
            if(!argumentDefinition.isHidden)
                argumentDefinitions.add(argumentDefinition);
        }

        // Try to fit the entire argument definition across the screen, but impose an arbitrary cap of 3/4 *
        // LINE_WIDTH in case the length of the arguments gets out of control.
        int argWidth = Math.min( findLongestArgumentCallingInfo(argumentDefinitions), (TextFormattingUtils.DEFAULT_LINE_WIDTH*3)/4 - FIELD_SEPARATION_WIDTH );
        int docWidth = TextFormattingUtils.DEFAULT_LINE_WIDTH - argWidth - FIELD_SEPARATION_WIDTH;

        for( ArgumentDefinition argumentDefinition: argumentDefinitions ) {
            Iterator<String> wordWrappedArgs = TextFormattingUtils.wordWrap( getArgumentCallingInfo(argumentDefinition), argWidth ).iterator();
            Iterator<String> wordWrappedDoc  = TextFormattingUtils.wordWrap( getArgumentDoc(argumentDefinition), docWidth ).iterator();

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
        if( !argumentDefinition.isFlag )
            formatter.format(" <%s>", argumentDefinition.fullName);

        return builder.toString();
    }

    /**
     * Gets a string of argument documentation.
     * @param argumentDefinition Argument definition for which help should be printed.
     * @return Brief description for this argument.
     */
    private String getArgumentDoc( ArgumentDefinition argumentDefinition ) {
        StringBuilder builder = new StringBuilder();
        builder.append(argumentDefinition.doc);
        if( argumentDefinition.validOptions != null ) {
            builder.append(" (");
            builder.append(Utils.join("|",argumentDefinition.validOptions));
            builder.append(")");
        }
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
        Comparator<ArgumentDefinition> definitionComparator = new Comparator<ArgumentDefinition>() {
            public int compare( ArgumentDefinition lhs, ArgumentDefinition rhs ) {
                if( lhs.required && rhs.required ) return 0;
                if( lhs.required ) return -1;
                if( rhs.required ) return 1;
                return 0;
            }
        };

        List<ArgumentDefinitionGroup> argumentGroups = new ArrayList<ArgumentDefinitionGroup>();
        for( ArgumentDefinitionGroup argumentGroup: argumentDefinitions.getArgumentDefinitionGroups() ) {
            List<ArgumentDefinition> sortedDefinitions = new ArrayList<ArgumentDefinition>( argumentGroup.argumentDefinitions );
            Collections.sort( sortedDefinitions, definitionComparator );
            argumentGroups.add( new ArgumentDefinitionGroup(argumentGroup.groupName,sortedDefinitions) );
        }

        // Sort the argument groups themselves with main arguments first, followed by plugins sorted in name order.
        Comparator<ArgumentDefinitionGroup> groupComparator = new Comparator<ArgumentDefinitionGroup>() {
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
     * @param parsedArgs the arguments passed in
     */
    public static void generateHeaderInformation(ApplicationDetails applicationDetails, Map<ArgumentMatchSource, ParsedArgs> parsedArgs) {

        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        java.util.Date date = new java.util.Date();

        String barrier = createBarrier(applicationDetails.applicationHeader);

        logger.info(barrier);
        for (String headerLine : applicationDetails.applicationHeader)
            logger.info(headerLine);
        logger.debug("Current directory: " + System.getProperty("user.dir"));
        for (Map.Entry<ArgumentMatchSource, ParsedArgs> entry: parsedArgs.entrySet()) {
            ArgumentMatchSource matchSource = entry.getKey();
            final String sourceName;
            switch (matchSource.getType()) {
                case CommandLine: sourceName = "Program"; break;
                case Provider: sourceName = matchSource.getDescription(); break;
                default: throw new RuntimeException("Unexpected argument match source type: " + matchSource.getType());
            }

            String output = sourceName + " Args: " + entry.getValue().getDescription();
            logger.info(output);
        }
        logger.info(generateUserHelpData());
        logger.info("Date/Time: " + dateFormat.format(date));
        logger.info(barrier);

        for(String attribution: applicationDetails.attribution)
            logger.info(attribution);
        logger.info(barrier);
    }

    /**
     * Create the user-related help information.
     * @return a non-null, non-empty String with the relevant information.
     */
    private static String generateUserHelpData() {
	try {
	    return "Executing as " +
		System.getProperty("user.name") + "@" + InetAddress.getLocalHost().getHostName() +
		" on " + System.getProperty("os.name") + " " + System.getProperty("os.version") +
		" " + System.getProperty("os.arch") + "; " + System.getProperty("java.vm.name") +
		" " + System.getProperty("java.runtime.version") + ".";
	} catch (Exception e) {
	    // don't fail
	    return "";
	}
    }

    /**
     * Create a barrier to use to distinguish the header from the rest of the output.
     * @param text A collection of lines to output as part of a header.
     * @return A barrier consisting of the '-' character.
     */
    private static String createBarrier(List<String> text) {
        int barrierWidth = 0;
        for(String headerLine: text)
            barrierWidth = Math.max(headerLine.length(),barrierWidth);
        return String.format("%0" + barrierWidth + "d",0).replace('0','-');
    }
}
