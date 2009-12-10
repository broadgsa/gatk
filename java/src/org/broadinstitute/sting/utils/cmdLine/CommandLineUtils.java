package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.JVMUtils;

import java.util.Map;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Collection;
import java.lang.reflect.Field;

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
