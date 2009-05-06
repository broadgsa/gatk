package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;

import java.util.Formatter;
import java.util.Locale;
import java.util.Formattable;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.HashMap;
import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;
import java.util.SortedSet;
import java.util.Collection;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
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
    /**
     * Target this line width.
     */
    private static final int LINE_WIDTH = 100;
    private static final int ARG_DOC_SEPARATION_WIDTH = 3;

    /**
     * Prints the help, given a collection of argument definitions.
     * @param argumentDefinitions Argument definitions for which help should be printed.
     */
    public void printHelp( ArgumentDefinitions argumentDefinitions ) {
        SortedSet<ArgumentDefinition> mainArguments = getMainArguments( argumentDefinitions );
        Map<String,SortedSet<ArgumentDefinition>> pluginsByGroup = getPluginArguments( argumentDefinitions );

        System.out.printf("%s%s%n",
                getSynopsis(mainArguments,pluginsByGroup),
                getDetailed(mainArguments,pluginsByGroup) );
    }

    /**
     * Retrieve a sorted set of the command-line arguments for the main application from the ArgumentDefinitions object.
     * @param argumentDefinitions Predefined argument definitions.
     * @return A list of argument definitions.
     */
    private SortedSet<ArgumentDefinition> getMainArguments( ArgumentDefinitions argumentDefinitions ) {
        for( ArgumentDefinitionGroup argumentGroup: argumentDefinitions.getArgumentDefinitionGroups() ) {
            // A null groupName is the signal that these are the application's main arguments.
            if( argumentGroup.groupName == null ){
                SortedSet<ArgumentDefinition> sortedDefinitions = new TreeSet<ArgumentDefinition>( new ArgumentDefinitionComparator() );
                sortedDefinitions.addAll( argumentGroup.argumentDefinitions );
                return sortedDefinitions;
            }
        }
        throw new StingException( "Unable to retrieve main command-line arguments." );
    }

    /**
     * Retrieve a sorted set of the command-line arguments for application plugin objects.
     * @param argumentDefinitions Predefined argument definitions.
     * @return A map of argument group name -> a list of argument definitions.
     */
    private Map<String,SortedSet<ArgumentDefinition>> getPluginArguments( ArgumentDefinitions argumentDefinitions ) {
        Map<String,SortedSet<ArgumentDefinition>> pluginsByGroup = new HashMap<String,SortedSet<ArgumentDefinition>>();
        for( ArgumentDefinitionGroup argumentGroup: argumentDefinitions.getArgumentDefinitionGroups() ) {
            // A null groupName is the signal that these are the application's main arguments.
            if( argumentGroup.groupName == null )
                continue;
            SortedSet<ArgumentDefinition> sortedDefinitions = new TreeSet<ArgumentDefinition>( new ArgumentDefinitionComparator() );
            sortedDefinitions.addAll( argumentGroup.argumentDefinitions );
            pluginsByGroup.put( argumentGroup.groupName, sortedDefinitions );
        }
        return pluginsByGroup;
    }

    /**
     * Gets the synopsis: the actual command to run.
     * @param mainArguments Main program arguments.
     * @praam pluginArgumentGroups Groups of plugin arguments
     * @return A synopsis line.
     */
    private String getSynopsis( SortedSet<ArgumentDefinition> mainArguments, Map<String,SortedSet<ArgumentDefinition>> pluginArgumentGroups ) {
        // Build out the synopsis all as one long line.        
        StringBuilder lineBuilder = new StringBuilder();
        Formatter lineFormatter = new Formatter( lineBuilder );

        lineFormatter.format("java -jar dist/GenomeAnalysisTK.jar");

        List<ArgumentDefinition> argumentDefinitions = new ArrayList<ArgumentDefinition>();
        argumentDefinitions.addAll( mainArguments );
        for( SortedSet pluginArguments: pluginArgumentGroups.values() )
            argumentDefinitions.addAll( pluginArguments );

        for( ArgumentDefinition argumentDefinition: argumentDefinitions ) {
            lineFormatter.format(" ");
            if( !argumentDefinition.required ) lineFormatter.format("[");
            if( argumentDefinition.shortName != null )
                lineFormatter.format("-%s", argumentDefinition.shortName);
            else
                lineFormatter.format("--%s", argumentDefinition.fullName);
            if( !argumentDefinition.isFlag() )
                lineFormatter.format(" <%s>", argumentDefinition.fullName);
            if( !argumentDefinition.required ) lineFormatter.format("]");
        }

        // Word wrap the synopsis.
        List<String> wrappedSynopsis = wordWrap( lineBuilder.toString(), LINE_WIDTH );

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
     * @param mainArguments Main program arguments.
     * @param pluginArgumentGroups Groups of plugin arguments
     * @return Detailed text about all arguments.
     */
    private String getDetailed( SortedSet<ArgumentDefinition> mainArguments, Map<String,SortedSet<ArgumentDefinition>> pluginArgumentGroups ) {
        StringBuilder builder = new StringBuilder();

        builder.append( getDetailForGroup( mainArguments ) );

        for( String groupName: pluginArgumentGroups.keySet() ) {
            builder.append( String.format("Arguments for %s:%n", groupName ) );
            builder.append( getDetailForGroup( pluginArgumentGroups.get(groupName) ) );
        }
        return builder.toString();
    }

    private String getDetailForGroup( SortedSet<ArgumentDefinition> argumentDefinitions ) {
        StringBuilder builder = new StringBuilder();
        Formatter formatter = new Formatter( builder );

        // Try to fit the entire argument definition across the screen, but impose an arbitrary cap of 3/4 *
        // LINE_WIDTH in case the length of the arguments gets out of control.
        int argWidth = Math.min( findLongestArgumentCallingInfo(argumentDefinitions), (LINE_WIDTH*3)/4 - ARG_DOC_SEPARATION_WIDTH );
        int docWidth = LINE_WIDTH - argWidth - ARG_DOC_SEPARATION_WIDTH;

        for( ArgumentDefinition argumentDefinition: argumentDefinitions ) {
            Iterator<String> wordWrappedArgs = wordWrap( getArgumentCallingInfo(argumentDefinition), argWidth ).iterator();
            Iterator<String> wordWrappedDoc  = wordWrap( argumentDefinition.doc, docWidth ).iterator();

            while( wordWrappedArgs.hasNext() || wordWrappedDoc.hasNext() ) {
                String arg = wordWrappedArgs.hasNext() ? wordWrappedArgs.next() : "";
                String doc = wordWrappedDoc.hasNext() ? wordWrappedDoc.next() : "";

                String formatString = "%-" + argWidth + "s%" + ARG_DOC_SEPARATION_WIDTH + "s%s%n";
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
        if( !argumentDefinition.isFlag() )
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
     * Simple implementation of word-wrap for a line of text.  Idea and
     * regexp shamelessly stolen from http://joust.kano.net/weblog/archives/000060.html.
     * Regexp can probably be simplified for our application.
     * @param text Text to wrap.
     * @param width Maximum line width.
     * @return A list of word-wrapped lines.
     */
    private List<String> wordWrap( String text, int width ) {
        Pattern wrapper = Pattern.compile( String.format(".{0,%d}(?:\\S(?: |$)|$)", width-1) );
        Matcher matcher = wrapper.matcher( text );

        List<String> wrapped = new ArrayList<String>();
        while( matcher.find() ) {
            // Regular expression is supersensitive to whitespace.
            // Assert that content is present before adding the line.
            String line = matcher.group().trim();
            if( line.length() > 0 )
                wrapped.add( matcher.group() );
        }
        return wrapped;
    }

    /**
     * A comparator to reorder arguments in alphabetical order.
     */
    private class ArgumentDefinitionComparator implements Comparator<ArgumentDefinition> {
        public int compare( ArgumentDefinition lhs, ArgumentDefinition rhs ) {
            if( lhs.shortName == null && rhs.shortName == null )
                return lhs.fullName.compareTo( rhs.fullName );
            // short names always come before long names.
            else if ( lhs.shortName == null )
                return -1;
            else if ( rhs.shortName == null )
                return 1;
            else if ( lhs.shortName.equals(rhs.shortName) )
                return lhs.fullName.compareTo( rhs.fullName );
            else
                return lhs.shortName.compareTo( rhs.shortName );
        }
    }
}
