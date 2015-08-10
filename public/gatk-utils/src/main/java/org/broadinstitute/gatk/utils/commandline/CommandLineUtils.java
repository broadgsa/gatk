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

package org.broadinstitute.gatk.utils.commandline;

import org.apache.log4j.Appender;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.lang.annotation.Annotation;
import java.util.Collections;
import java.util.Enumeration;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Static utility methods for working with command-line arguments.
 *
 * @author mhanna
 * @version 0.1
 */
public class CommandLineUtils {

    /**
     * Returns a key-value mapping of the command-line arguments passed into the GATK.
     * Will be approximate; this class doesn't have all the required data to completely
     * reconstruct the list of command-line arguments from the given objects.
     *
     * @param parsingEngine      The parsing engine
     * @param argumentProviders  The providers of command-line arguments.
     * @return A key-value mapping of argument full names to argument values.  Produces best string representation
     *         possible given the information available.
     */
    public static Map<String,String> getApproximateCommandLineArguments(ParsingEngine parsingEngine, Object... argumentProviders) {
        return getApproximateCommandLineArguments(parsingEngine, false, argumentProviders);
    }

    /**
     * Returns a key-value mapping of the command-line arguments passed into the GATK.
     * Will be approximate; this class doesn't have all the required data to completely
     * reconstruct the list of command-line arguments from the given objects.
     * 
     * @param parsingEngine      The parsing engine
     * @param skipObjectPointers Should we skip arguments whose values are pointers (and don't print nicely)?
     * @param argumentProviders  The providers of command-line arguments.
     * @return A key-value mapping of argument full names to argument values.  Produces best string representation
     *         possible given the information available.
     */
    public static Map<String,String> getApproximateCommandLineArguments(ParsingEngine parsingEngine, boolean skipObjectPointers, Object... argumentProviders) {
        Map<String,String> commandLineArguments = new LinkedHashMap<String,String>();

        for(Object argumentProvider: argumentProviders) {
            Map<ArgumentSource, Object> argBindings = parsingEngine.extractArgumentBindings(argumentProvider);
            for(Map.Entry<ArgumentSource, Object> elt: argBindings.entrySet()) {
                Object argumentValue = elt.getValue();

                String argumentValueString = argumentValue != null ? argumentValue.toString() : null;
                if ( skipObjectPointers && isObjectPointer(argumentValueString) )
                    continue;

                for(ArgumentDefinition definition: elt.getKey().createArgumentDefinitions()) {
                    String argumentName = definition.fullName;
                    commandLineArguments.put(argumentName,argumentValueString);
                }
            }
        }

        return commandLineArguments;
    }

    /**
     * Create an approximate list of command-line arguments based on the given argument providers.
     * @param parsingEngine      The parsing engine
     * @param argumentProviders  Argument providers to inspect.
     * @return A string representing the given command-line arguments.
     */
    public static String createApproximateCommandLineArgumentString(ParsingEngine parsingEngine, Object... argumentProviders) {
        return createApproximateCommandLineArgumentString(parsingEngine, true, argumentProviders);
    }

    /**
     * Create an approximate list of command-line arguments based on the given argument providers.
     * @param parsingEngine      The parsing engine
     * @param skipObjectPointers Should we skip arguments whose values are pointers (and don't print nicely)?
     * @param argumentProviders  Argument providers to inspect.
     * @return A string representing the given command-line arguments.
     */
    public static String createApproximateCommandLineArgumentString(ParsingEngine parsingEngine, boolean skipObjectPointers, Object... argumentProviders) {
        Map<String,String> commandLineArgs = getApproximateCommandLineArguments(parsingEngine, skipObjectPointers, argumentProviders);
        StringBuffer sb = new StringBuffer();

        boolean first = true;
        for ( Map.Entry<String, String> commandLineArg : commandLineArgs.entrySet() ) {
            if ( !first )
                sb.append(" ");
            sb.append(commandLineArg.getKey());
            sb.append("=");
            sb.append(commandLineArg.getValue());
            first = false;
        }

        return sb.toString();
    }

    /**
     * A hack to get around the fact that Java doesn't like inheritance in Annotations.
     * @param annotation to run the method on
     * @param method the method to invoke
     * @return the return value of the method
     */
    public static Object getValue(Annotation annotation, String method) {
        try {
            return annotation.getClass().getMethod(method).invoke(annotation);
        } catch (Exception e) {
            throw new ReviewedGATKException("Unable to access method " + method + " on annotation " + annotation.getClass(), e);
        }
    }

    // The problem here is that some of the fields being output are Objects - and those
    //  Objects don't overload toString() so that the output is just the memory pointer
    //  to the Object.  Because those values are non-deterministic, they don't merge well
    //  into BAM/VCF headers (plus, it's just damn ugly).  Perhaps there's a better way to
    //  do this, but at least this one works for the moment.
    private static final String pointerRegexp = ".+@[0-9a-fA-F]+$";
    private static boolean isObjectPointer(String s) {
        return s != null && s.matches(pointerRegexp);
    }

    /**
     * Returns the root logger for all GATK code.
     * @return the root logger for all GATK  code.
     */
    public static Logger getStingLogger() {
        return Logger.getLogger("org.broadinstitute.gatk");
    }

    /**
     * Enables console logging.
     */
    @SuppressWarnings("unchecked")
    public static void configureConsoleLogging() {
        // Check to see if a console logger has already been enabled.
        for (Logger logger = getStingLogger(); logger != null; logger = (Logger)logger.getParent()) {
            Enumeration<Appender> e = (Enumeration<Appender>) logger.getAllAppenders();
            for (Appender appender: Collections.list(e)) {
                if (appender instanceof ConsoleAppender)
                    return;
            }
        }
        // Extracted from BasicConfigurator.configure(), but only applied to the GATK logger.
        Logger.getRootLogger().addAppender(new ConsoleAppender(
                    new PatternLayout(PatternLayout.TTCC_CONVERSION_PATTERN), ConsoleAppender.SYSTEM_ERR));
    }

    /**
     * Sets the layout of the logger.
     * @param logger The logger.
     * @param layout The layout.
     */
    @SuppressWarnings("unchecked")
    public static void setLayout(Logger logger, PatternLayout layout) {
        for (; logger != null; logger = (Logger)logger.getParent()) {
            Enumeration<Appender> e = (Enumeration<Appender>) logger.getAllAppenders();
            for (Appender appender: Collections.list(e))
                appender.setLayout(layout);
        }
    }
}
