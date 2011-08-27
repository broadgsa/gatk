/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.R;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.utils.PathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Generic service for executing RScripts in the GATK directory
 *
 * @author Your Name
 * @since Date created
 */
public class RScriptExecutor {
    /**
     * our log
     */
    protected static Logger logger = Logger.getLogger(RScriptExecutor.class);

    public static class RScriptArgumentCollection {
        @Advanced
        @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is maybe /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
        public String PATH_TO_RSCRIPT = "Rscript";

        @Advanced
        @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
        public List<String> PATH_TO_RESOURCES = Arrays.asList("public/R/", "private/R/");
    }

    final RScriptArgumentCollection myArgs;
    final boolean exceptOnError;

    public RScriptExecutor(final RScriptArgumentCollection myArgs, final boolean exceptOnError) {
        this.myArgs = myArgs;
        this.exceptOnError = exceptOnError;
    }

    public void callRScripts(String scriptName, String... scriptArgs) {
        callRScripts(scriptName, Arrays.asList(scriptArgs));
    }

    public void callRScripts(String scriptName, List<String> scriptArgs) {
        try {
            final File pathToScript = findScript(scriptName);
            if ( pathToScript == null ) return; // we failed but shouldn't exception out
            final String argString = Utils.join(" ", scriptArgs);
            final String cmdLine = Utils.join(" ", Arrays.asList(myArgs.PATH_TO_RSCRIPT, pathToScript, argString));
            logger.info("Executing RScript: " + cmdLine);
            Runtime.getRuntime().exec(cmdLine).waitFor();
        } catch (InterruptedException e) {
            generateException(e);
        } catch (IOException e) {
            generateException("Fatal Exception: Perhaps RScript jobs are being spawned too quickly?", e);
        }
    }

    public File findScript(final String scriptName) {
        for ( String pathToResource : myArgs.PATH_TO_RESOURCES ) {
            final File f = new File(pathToResource + "/" + scriptName);
            if ( f.exists() ) {
                if ( f.canRead() )
                    return f;
                else
                    generateException("Script exists but couldn't be read: " + scriptName);
            }
        }

        generateException("Couldn't find script: " + scriptName + " in " + myArgs.PATH_TO_RESOURCES);
        return null;
    }

    private void generateException(String msg) {
        generateException(msg, null);
    }

    private void generateException(Throwable e) {
        generateException("", e);
    }

    private void generateException(String msg, Throwable e) {
        if ( exceptOnError )
            throw new UserException(msg, e);
        else
            logger.warn(msg + (e == null ? "" : ":" + e.getMessage()));
    }
}
