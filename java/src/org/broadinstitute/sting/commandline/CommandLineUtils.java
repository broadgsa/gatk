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

import org.broadinstitute.sting.utils.classloader.JVMUtils;

import java.util.Map;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Collection;

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
            List<ArgumentSource> argumentSources = ParsingEngine.extractArgumentSources(argumentProvider.getClass());
            for(ArgumentSource argumentSource: argumentSources) {
                Object argumentValue = JVMUtils.getFieldValue(argumentSource.field,argumentProvider);
                String argumentValueString = argumentValue != null ? argumentValue.toString() : null;

                for(ArgumentDefinition definition: argumentSource.createArgumentDefinitions()) {
                    String argumentName = definition.fullName;
                    commandLineArguments.put(argumentName,argumentValueString);
                }
            }
        }

        return commandLineArguments;
    }
}
