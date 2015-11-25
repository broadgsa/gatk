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

import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.help.ApplicationDetails;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.help.HelpFormatter;
import org.broadinstitute.gatk.utils.text.TextFormattingUtils;

import java.io.IOException;
import java.util.*;

public abstract class CommandLineProgram {

    /** The command-line program and the arguments it returned. */
    public ParsingEngine parser = null;

    /**
     * Setting INFO gets you INFO up to FATAL, setting ERROR gets you ERROR and FATAL level logging, and so on.
     */
    @Argument(fullName = "logging_level", shortName = "l", doc = "Set the minimum level of logging", required = false)
    protected String logging_level = "INFO";

    /**
     * File to save the logging output.
     */
    @Output(fullName = "log_to_file", shortName = "log", doc = "Set the logging location", required = false)
    protected String toFile = null;

    /**
     * This will produce a help message in the terminal with general usage information, listing available arguments
     * as well as tool-specific information if applicable.
     */
    @Argument(fullName = "help", shortName = "h", doc = "Generate the help message", required = false)
    public Boolean help = false;

    /**
     * Use this to check the version number of the GATK executable you are invoking. Note that the version number is
     * always included in the output at the start of every run as well as any error message.
     */
    @Argument(fullName = "version", shortName = "version", doc ="Output version information", required = false)
    public Boolean version = false;


    /** our logging output patterns */
    private static final String patternString = "%-5p %d{HH:mm:ss,SSS} %C{1} - %m %n";

    static {
        /**
         * The very first thing that any GATK application does is forces the JVM locale into US English, so that we don't have
         * to think about number formatting issues.
         */
        forceJVMLocaleToUSEnglish();
        // setup a basic log configuration
        CommandLineUtils.configureConsoleLogging();
    }


    /**
     * Allows a given application to return a brief description of itself.
     *
     * @return An ApplicationDetails object describing the current application.  Should not be null.
     */
    protected ApplicationDetails getApplicationDetails() {
        return new ApplicationDetails(ApplicationDetails.createDefaultHeader(getClass()),
                                      Collections.<String>emptyList(),
                                      ApplicationDetails.createDefaultRunningInstructions(getClass()),
                                      null);
    }

    /**
     * Subclasses of CommandLinePrograms can provide their own types of command-line arguments.
     * @return A collection of type descriptors generating implementation-dependent placeholders.
     */
    protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
        return Collections.emptyList();
    }

    /**
     * Will this application want to vary its argument list dynamically?
     * If so, parse the command-line options and then prompt the subclass to return
     * a list of argument providers.
     *
     * @return Whether the application should vary command-line arguments dynamically.
     */
    protected boolean canAddArgumentsDynamically() { return false; }

    /**
     * Provide a list of object to inspect, looking for additional command-line arguments.
     *
     * @return A list of objects to inspect.
     */
    protected Class[] getArgumentSources() {
        return new Class[]{};
    }

    /**
     * Name this argument source.  Provides the (full) class name as a default.
     *
     * @param source The argument source.
     *
     * @return a name for the argument source.
     */
    protected String getArgumentSourceName( Class source ) { return source.toString(); }

    /**
     * Sets the command-line parsing engine. Necessary for unit testing purposes.
     * @param parser the new command-line parsing engine
     */
    public void setParser( ParsingEngine parser ) {
        this.parser = parser;
    }

    /**
     * this is the function that the inheriting class can expect to have called
     * when all the argument processing is done
     *
     * @return the return code to exit the program with
     * @throws Exception when an exception occurs
     */
    protected abstract int execute() throws Exception;

    public static int result = -1;

    @SuppressWarnings("unchecked")
    public static void start(CommandLineProgram clp, String[] args) throws Exception {
        start(clp, args, false);
    }

    /**
     * This function is called to start processing the command line, and kick
     * off the execute message of the program.
     *
     * @param clp  the command line program to execute
     * @param args the command line arguments passed in
     * @param dryRun dry run
     * @throws Exception when an exception occurs
     */
    @SuppressWarnings("unchecked")
    public static void start(CommandLineProgram clp, String[] args, boolean dryRun) throws Exception {

        try {
            // setup our log layout
            PatternLayout layout = new PatternLayout();

            Logger logger = CommandLineUtils.getStingLogger();

            // now set the layout of all the loggers to our layout
            CommandLineUtils.setLayout(logger, layout);

            // Initialize the logger using the defaults.
            clp.setupLoggerLevel(layout);

            // setup the parser
            ParsingEngine parser = clp.parser = new ParsingEngine(clp);
            parser.addArgumentSource(clp.getClass());

            Map<ArgumentMatchSource, ParsedArgs> parsedArgs;

            // process the args
            if (clp.canAddArgumentsDynamically()) {
                // if the command-line program can toss in extra args, fetch them and reparse the arguments.
                parser.parse(args);

                // Allow invalid and missing required arguments to pass this validation step.
                //   - InvalidArgument in case these arguments are specified by plugins.
                //   - MissingRequiredArgument in case the user requested help.  Handle that later, once we've
                //                             determined the full complement of arguments.
                if ( ! dryRun )
                    parser.validate(EnumSet.of(ParsingEngine.ValidationType.MissingRequiredArgument,
                            ParsingEngine.ValidationType.InvalidArgument));
                parser.loadArgumentsIntoObject(clp);

                // Initialize the logger using the loaded command line.
                clp.setupLoggerLevel(layout);

                Class[] argumentSources = clp.getArgumentSources();
                    for (Class argumentSource : argumentSources)
                    parser.addArgumentSource(clp.getArgumentSourceName(argumentSource), argumentSource);
                parsedArgs = parser.parse(args);

                if (isVersionPresent(parser))
                    printVersionAndExit();

                if (isHelpPresent(parser))
                    printHelpAndExit(clp, parser);

                if ( ! dryRun ) parser.validate();
            } else {
                parsedArgs = parser.parse(args);

                if ( ! dryRun ) {
                    if (isHelpPresent(parser))
                        printHelpAndExit(clp, parser);

                    parser.validate();
                }
                parser.loadArgumentsIntoObject(clp);

                // Initialize the logger using the loaded command line.
                clp.setupLoggerLevel(layout);
            }

            if ( ! dryRun ) {
                // if they specify a log location, output our data there
                if (clp.toFile != null) {
                    FileAppender appender;
                    try {
                        appender = new FileAppender(layout, clp.toFile, false);
                        logger.addAppender(appender);
                    } catch (IOException e) {
                        throw new RuntimeException("Unable to re-route log output to " + clp.toFile + " make sure the destination exists");
                    }
                }

                // regardless of what happens next, generate the header information
                HelpFormatter.generateHeaderInformation(clp.getApplicationDetails(), parsedArgs);

                // call the execute
                CommandLineProgram.result = clp.execute();
            }
        }
        catch (ArgumentException e) {
            //clp.parser.printHelp(clp.getApplicationDetails());
            // Rethrow the exception to exit with an error.
            throw e;
        }
    }

    /**
     * Find fields in the object obj that look like command-line arguments, and put command-line
     * arguments into them.
     *
     * @param obj Object to inspect for command line arguments.
     */
    public void loadArgumentsIntoObject(Object obj) {
        parser.loadArgumentsIntoObject(obj);
    }

    /**
     * this function checks the logger level passed in on the command line, taking the lowest
     * level that was provided.
     * @param layout Pattern layout to format based on the logger level.
     */
    private void setupLoggerLevel(PatternLayout layout) {
        layout.setConversionPattern(patternString);

        // set the default logger level
        Level par;
        if (logging_level.toUpperCase().equals("DEBUG")) {
            par = Level.DEBUG;
        } else if (logging_level.toUpperCase().equals("INFO")) {
            par = Level.INFO;
        } else if (logging_level.toUpperCase().equals("WARN")) {
            par = Level.WARN;
        } else if (logging_level.toUpperCase().equals("ERROR")) {
            par = Level.ERROR;
        } else if (logging_level.toUpperCase().equals("FATAL")) {
            par = Level.FATAL;
        } else if (logging_level.toUpperCase().equals("OFF")) {
            par = Level.OFF;
        } else {
            // we don't understand the logging level, let's get out of here
            throw new ArgumentException("Unable to match: " + logging_level + " to a logging level, make sure it's a valid level (DEBUG, INFO, WARN, ERROR, FATAL, OFF)");
        }

        Logger.getRootLogger().setLevel(par);
    }

    public static String getVersionNumber() {
        ResourceBundle headerInfo = TextFormattingUtils.GATK_RESOURCE_BUNDLE;
        return headerInfo.containsKey("org.broadinstitute.gatk.utils.version") ? headerInfo.getString("org.broadinstitute.gatk.utils.version") : "<unknown>";
    }

    public static String getBuildTime() {
        ResourceBundle headerInfo = TextFormattingUtils.GATK_RESOURCE_BUNDLE;
        return headerInfo.containsKey("build.timestamp") ? headerInfo.getString("build.timestamp") : "<unknown>";
    }

    /**
     * a function used to indicate an error occurred in the command line tool
     */
    private static void printDocumentationReference() {
        errorPrintf("Visit our website and forum for extensive documentation and answers to %n");
        errorPrintf("commonly asked questions " + HelpConstants.BASE_GATK_URL + "%n");
    }


    /**
     * Do a cursory search for the given argument.
     *
     * @param parser Parser
     *
     * @return True if help is present; false otherwise.
     */
    private static boolean isHelpPresent(ParsingEngine parser) {
        return parser.isArgumentPresent("help");
    }

    /**
     * Print help and exit.
     *
     * @param clp    Instance of the command-line program.
     * @param parser True if help is present; false otherwise.
     */
    private static void printHelpAndExit(CommandLineProgram clp, ParsingEngine parser) {
        parser.printHelp(clp.getApplicationDetails());
        System.exit(0);
    }

    /**
     * Do a cursory search for the argument "version".
     *
     * @param parser Parser
     *
     * @return True if version is present; false otherwise.
     */
    private static boolean isVersionPresent(ParsingEngine parser) {
        return parser.isArgumentPresent("version");
    }

    /**
     * Print help and exit.
     */
    private static void printVersionAndExit() {
        System.out.println(getVersionNumber().toString());
        System.exit(0);
    }


    private static void errorPrintf(String format, Object... s) {
        String formatted = String.format(format, s);

        if ( formatted.trim().equals("") )
            System.err.println("##### ERROR");
        else {
            for ( String part : formatted.split("\n") ) {
                System.err.println("##### ERROR " + part);
            }
        }
    }


    /**
     * used to indicate an error occured
     *
     * @param msg the message
     * @param t   the error
     */
    public static void exitSystemWithError(String msg, final Throwable t) {
        errorPrintf("------------------------------------------------------------------------------------------%n");
        errorPrintf("stack trace %n");
        t.printStackTrace();

        errorPrintf("------------------------------------------------------------------------------------------%n");
        errorPrintf("A GATK RUNTIME ERROR has occurred (version %s):%n", getVersionNumber());
        errorPrintf("%n");
        errorPrintf("This might be a bug. Please check the documentation guide to see if this is a known problem.%n");
        errorPrintf("If not, please post the error message, with stack trace, to the GATK forum.%n");
        printDocumentationReference();
        if ( msg == null ) // some exceptions don't have detailed messages
            msg = "Code exception (see stack trace for error itself)";
        errorPrintf("%n");
        errorPrintf("MESSAGE: %s%n", msg.trim());
        errorPrintf("------------------------------------------------------------------------------------------%n");
        System.exit(1);
    }

    public static void exitSystemWithUserError(final Exception e) {
        if ( e.getMessage() == null )
            throw new ReviewedGATKException("UserException found with no message!", e);

        errorPrintf("------------------------------------------------------------------------------------------%n");
        errorPrintf("A USER ERROR has occurred (version %s): %n", getVersionNumber());
        errorPrintf("%n");
        errorPrintf("This means that one or more arguments or inputs in your command are incorrect.%n");
        errorPrintf("The error message below tells you what is the problem.%n");
        errorPrintf("%n");
        errorPrintf("If the problem is an invalid argument, please check the online documentation guide%n");
        errorPrintf("(or rerun your command with --help) to view allowable command-line arguments for this tool.%n");
        errorPrintf("%n");
        printDocumentationReference();
        errorPrintf("%n");
        errorPrintf("Please do NOT post this error to the GATK forum unless you have really tried to fix it yourself.%n");
        errorPrintf("%n");
        errorPrintf("MESSAGE: %s%n", e.getMessage().trim());
        errorPrintf("------------------------------------------------------------------------------------------%n");
        System.exit(1);
    }

    public static void exitSystemWithSamError(final Throwable t) {
        if ( t.getMessage() == null )
            throw new ReviewedGATKException("SamException found with no message!", t);

        errorPrintf("------------------------------------------------------------------------------------------%n");
        errorPrintf("A BAM/CRAM ERROR has occurred (version %s): %n", getVersionNumber());
        errorPrintf("%n");
        errorPrintf("This means that there is something wrong with the BAM/CRAM file(s) you provided.%n");
        errorPrintf("The error message below tells you what is the problem.%n");
        errorPrintf("%n");
        printDocumentationReference();
        errorPrintf("%n");
        errorPrintf("Please do NOT post this error to the GATK forum until you have followed these instructions:%n");
        errorPrintf("- Make sure that your BAM file is well-formed by running Picard's validator on it%n");
        errorPrintf("(see http://picard.sourceforge.net/command-line-overview.shtml#ValidateSamFile for details)%n");
        errorPrintf("- Ensure that your BAM index is not corrupted: delete the current one and regenerate it with 'samtools index'%n");
        errorPrintf("- Ensure that your CRAM index is not corrupted: delete the current one and regenerate it with%n");
        errorPrintf("'java -jar cramtools-3.0.jar index --bam-style-index --input-file <input cram file> --reference-fasta-file <reference fasta file>'%n");
        errorPrintf("(see https://github.com/enasequence/cramtools/tree/v3.0 for details)%n");
        errorPrintf("%n");
        errorPrintf("MESSAGE: %s%n", t.getMessage().trim());
        errorPrintf("------------------------------------------------------------------------------------------%n");
        System.exit(1);
    }


    /**
     * used to indicate an error occured
     *
     * @param t the exception that occurred
     */
    public static void exitSystemWithError(Throwable t) {
        exitSystemWithError(t.getMessage(), t);
    }

    /**
     * A hack to ensure that numbers are always formatted in the US style.
     */
    protected static void forceJVMLocaleToUSEnglish() {
        Locale.setDefault(Locale.US);
    }
}
