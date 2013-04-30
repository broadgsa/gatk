/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.tools;

import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.help.HelpUtils;

/**
 * Utility program to print a list of available annotations
 *
 * <p>This is a very simple utility tool that retrieves available annotations for use with tools such as
 * UnifiedGenotyper, HaplotypeCaller and VariantAnnotator.</p>
 *
 * <h3>Important note</h3>
 * <p>This is a command-line utility that bypasses the GATK engine. As a result, the command-line you must use to
 * invoke it is a little different from other GATK tools (see usage below), and it does not accept any of the
 * classic "CommandLineGATK" arguments.</p>
 *
 * <h3>Usage</h3>
 * <pre>java -cp GenomeAnalysisTK.jar org.broadinstitute.sting.tools.ListAnnotations</pre>
 *
 * @author vdauwera
 * @since 3/14/13
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_HELPUTILS )
public class ListAnnotations extends CommandLineProgram {

    /*
     * Print usage information
     *
     * TODO: would be more convenient if we could just call the program by name instead of the full classpath
     */
    private static void printUsage() {
        System.err.println("Usage: java -cp dist/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.ListAnnotations");
        System.err.println("    Prints a list of available annotations and exits.");
    }

    // TODO: override CommandLineProgram bit that offers version, logging etc arguments. We don't need that stuff here and it makes the doc confusing.

    @Override
    protected int execute() throws Exception {

        HelpUtils.listAnnotations();
        return 0;
    }

    public static void main(String[] args){
        try {
            ListAnnotations instance = new ListAnnotations();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch ( UserException e ) {
            printUsage();
            exitSystemWithUserError(e);
        } catch ( Exception e ) {
            exitSystemWithError(e);
        }
    }
}
