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
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.io.IOUtils;
import org.broadinstitute.sting.utils.io.Resource;
import org.broadinstitute.sting.utils.runtime.ProcessController;
import org.broadinstitute.sting.utils.runtime.ProcessSettings;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Generic service for executing RScripts
 */
public class RScriptExecutor {
    /**
     * our log
     */
    private static Logger logger = Logger.getLogger(RScriptExecutor.class);

    public static class RScriptArgumentCollection {
        @Advanced
        @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. Defaults Rscript meaning to use the first available on the environment PATH. For Broad users should 'use R-2.12' or later.", required = false)
        public String PATH_TO_RSCRIPT = "Rscript";

        @Advanced
        @Argument(fullName = "path_to_Rresources", shortName = "Rresources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
        public List<String> PATH_TO_RESOURCES = Arrays.asList("public/R/", "private/R/");

        public RScriptArgumentCollection() {}

        /* For testing and convenience */
        public RScriptArgumentCollection(final String PATH_TO_RSCRIPT, final List<String> PATH_TO_RESOURCES) {
            this.PATH_TO_RSCRIPT = PATH_TO_RSCRIPT;
            this.PATH_TO_RESOURCES = PATH_TO_RESOURCES;
        }
    }

    private final RScriptArgumentCollection myArgs;
    private final boolean exceptOnError;
    private final List<RScriptLibrary> libraries = new ArrayList<RScriptLibrary>();
    private final List<Resource> scriptResources = new ArrayList<Resource>();
    private final List<File> scriptFiles = new ArrayList<File>();
    private final List<String> args = new ArrayList<String>();

    public RScriptExecutor(final RScriptArgumentCollection myArgs, final boolean exceptOnError) {
        this.myArgs = myArgs;
        this.exceptOnError = exceptOnError;
    }

    public void addLibrary(RScriptLibrary library) {
        this.libraries.add(library);
    }

    public void addScript(Resource script) {
        this.scriptResources.add(script);
    }

    public void addScript(File script) {
        this.scriptFiles.add(script);
    }

    /**
     * Adds args to the end of the Rscript command line.
     * @param args the args.
     * @throws NullPointerException if any of the args are null.
     */
    public void addArgs(Object... args) {
        for (Object arg: args)
            this.args.add(arg.toString());
    }

    public void exec() {
        List<File> tempFiles = new ArrayList<File>();
        try {
            File tempLibDir = IOUtils.tempDir("R.", ".lib");
            tempFiles.add(tempLibDir);

            StringBuilder expression = new StringBuilder("tempLibDir = '").append(tempLibDir).append("';");

            if (this.libraries.size() > 0) {
                List<String> tempLibraryPaths = new ArrayList<String>();
                for (RScriptLibrary library: this.libraries) {
                    File tempLibrary = library.writeTemp();
                    tempFiles.add(tempLibrary);
                    tempLibraryPaths.add(tempLibrary.getAbsolutePath());
                }

                expression.append("install.packages(");
                expression.append("pkgs=c('").append(StringUtils.join(tempLibraryPaths, "', '")).append("'), lib=tempLibDir, repos=NULL, type='source', ");
                // Install faster by eliminating cruft.
                expression.append("INSTALL_opts=c('--no-libs', '--no-data', '--no-help', '--no-demo', '--no-exec')");
                expression.append(");");

                for (RScriptLibrary library: this.libraries) {
                    expression.append("require('").append(library.getLibraryName()).append("', lib.loc=tempLibDir);");
                }
            }

            for (Resource script: this.scriptResources) {
                File tempScript = IOUtils.writeTempResource(script);
                tempFiles.add(tempScript);
                expression.append("source('").append(tempScript.getAbsolutePath()).append("');");
            }

            for (File script: this.scriptFiles) {
                expression.append("source('").append(script.getAbsolutePath()).append("');");
            }

            String[] cmd = new String[this.args.size() + 3];
            int i = 0;
            cmd[i++] = myArgs.PATH_TO_RSCRIPT;
            cmd[i++] = "-e";
            cmd[i++] = expression.toString();
            for (String arg: this.args)
                cmd[i++] = arg;

            ProcessSettings processSettings = new ProcessSettings(cmd);
            if (logger.isDebugEnabled()) {
                processSettings.getStdoutSettings().printStandard(true);
                processSettings.getStderrSettings().printStandard(true);
            }

            ProcessController controller = ProcessController.getThreadLocal();

            logger.debug("Executing: " + Utils.join(" ", cmd));
            logger.debug("Result: " + controller.exec(processSettings).getExitValue());

        } catch (StingException e) {
            generateException(e);
        } finally {
            for (File temp: tempFiles)
                FileUtils.deleteQuietly(temp);
        }
    }

    public void callRScripts(String scriptName, Object... scriptArgs) {
        final File pathToScript = findScript(scriptName);
        if (pathToScript == null) return; // we failed but shouldn't exception out
        addScript(pathToScript);
        addArgs(scriptArgs);
        exec();
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
