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

import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.classloader.JVMUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * Contains details additional details that the program can
 * supply about itself.
 *
 * @author hanna
 * @version 0.1
 */

public class ApplicationDetails {
    /**
     * Retrieve key information about the application (name, who to contact for support, etc.).
     */
    final List<String> applicationHeader;

    /**
     * Stores additional attribution for a given walker.
     */
    final List<String> attribution;

    /**
     * Extract details covering exactly how to run this executable.
     */
    final String runningInstructions;

    /**
     * Additional help particular to this command-line application.
     */
    final String additionalHelp;

    public ApplicationDetails( List<String> applicationHeader, List<String> attribution, String runningInstructions, String additionalHelp ) {
        this.applicationHeader = applicationHeader;
        this.attribution = attribution;
        this.runningInstructions = runningInstructions;
        this.additionalHelp = additionalHelp;
    }

    public static List<String> createDefaultHeader(Class<? extends CommandLineProgram> application) {
        return Collections.singletonList("Program Name: " + application.getName());
    }

    public static String createDefaultRunningInstructions(Class<? extends CommandLineProgram> application) {
        // Default implementation to find a command line that makes sense.
        // If the user is running from a jar, return '-jar <jarname>'; otherwise
        // return the full class name.
        String runningInstructions = null;
        try {
            runningInstructions = JVMUtils.getLocationFor( application ).getName();
        }
        catch( IOException ex ) {
            throw new ReviewedGATKException("Unable to determine running instructions", ex);
        }

        if( runningInstructions.endsWith(".jar") )
            runningInstructions = String.format("-jar %s", runningInstructions);
        else
            runningInstructions = application.getName();

        return runningInstructions;
    }
}
