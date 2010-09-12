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

package org.broadinstitute.sting.commandline;

import org.broadinstitute.sting.utils.exceptions.GATKException;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.util.*;
import java.lang.annotation.Annotation;

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
     * @param argumentProviders The providers of command-line arguments.
     * @return A key-value mapping of argument full names to argument values.  Produces best string representation
     *         possible given the information available.
     */
    public static Map<String,String> getApproximateCommandLineArguments(Collection<Object> argumentProviders) {
        Map<String,String> commandLineArguments = new LinkedHashMap<String,String>();

        for(Object argumentProvider: argumentProviders) {
            Map<ArgumentSource, Object> argBindings = ParsingEngine.extractArgumentBindings(argumentProvider);
            for(Map.Entry<ArgumentSource, Object> elt: argBindings.entrySet()) {
                Object argumentValue = elt.getValue();
                String argumentValueString = argumentValue != null ? argumentValue.toString() : null;

                for(ArgumentDefinition definition: elt.getKey().createArgumentDefinitions()) {
                    String argumentName = definition.fullName;
                    commandLineArguments.put(argumentName,argumentValueString);
                }
            }
        }

        return commandLineArguments;
    }

//    public static Map<String,String> getApproximateCommandLineArguments(Collection<Object> argumentProviders) {
//        Map<String,String> commandLineArguments = new LinkedHashMap<String,String>();
//
//        for(Object argumentProvider: argumentProviders) {
//            Map<ArgumentSource, Object> argBings = ParsingEngine.extractArgumentBindings(argumentProvider);
//            List<ArgumentSource> argumentSources = ParsingEngine.extractArgumentSources(argumentProvider.getClass());
//            for(ArgumentSource argumentSource: argumentSources) {
//                Object argumentValue = JVMUtils.getFieldValue(argumentSource.field,argumentProvider);
//                String argumentValueString = argumentValue != null ? argumentValue.toString() : null;
//
//                for(ArgumentDefinition definition: argumentSource.createArgumentDefinitions()) {
//                    String argumentName = definition.fullName;
//                    commandLineArguments.put(argumentName,argumentValueString);
//                }
//            }
//        }
//
//        return commandLineArguments;
//    }

    public static String createApproximateCommandLineArgumentString(GenomeAnalysisEngine toolkit, Walker walker) {
        return createApproximateCommandLineArgumentString(toolkit, null, walker);
    }

    public static String createApproximateCommandLineArgumentString(GenomeAnalysisEngine toolkit, Collection<Object> otherArgumentProviders, Walker walker) {
        StringBuffer sb = new StringBuffer();
        sb.append("analysis_type=");
        sb.append(toolkit.getWalkerName(walker.getClass()));

        ArrayList<Object> allArgumentProviders = new ArrayList<Object>();
        allArgumentProviders.add(toolkit.getArguments());
        allArgumentProviders.add(walker);
        if (otherArgumentProviders != null) allArgumentProviders.addAll(otherArgumentProviders);

        Map<String,String> commandLineArgs = getApproximateCommandLineArguments(allArgumentProviders);

        for ( Map.Entry<String, String> commandLineArg : commandLineArgs.entrySet() ) {
            sb.append(" ");
            sb.append(commandLineArg.getKey());
            sb.append("=");
            sb.append(commandLineArg.getValue());
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
            throw new GATKException("Unable to access method " + method + " on annotation " + annotation.getClass(), e);
        }
    }
}
