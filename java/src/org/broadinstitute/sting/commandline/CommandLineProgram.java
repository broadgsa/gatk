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

import org.apache.log4j.*;
import org.broadinstitute.sting.utils.help.ApplicationDetails;
import org.broadinstitute.sting.utils.help.HelpFormatter;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public abstract class CommandLineProgram {

    /** The command-line program and the arguments it returned. */
    private ParsingEngine parser = null;

    /** our log, which we want to capture anything from org.broadinstitute.sting */
    private static Logger logger = Logger.getRootLogger();

    /** the default log level */
    @Argument(fullName = "logging_level",
              shortName = "l",
              doc = "Set the minimum level of logging, i.e. setting INFO get's you INFO up to FATAL, setting ERROR gets you ERROR and FATAL level logging.",
              required = false)
    protected String logging_level = "INFO";


    /** where to send the output of our logger */
    @Output(fullName = "log_to_file",
              shortName = "log",
              doc = "Set the logging location",
              required = false)
    protected String toFile = null;

    /** do we want to silence the command line output */
    @Argument(fullName = "quiet_output_mode",
              shortName = "quiet",
              doc = "Set the logging to quiet mode, no output to stdout",
              required = false)
    protected Boolean quietMode = false;

    /** do we want to generate debugging information with the logs */
    @Argument(fullName = "debug_mode",
              shortName = "debug",
              doc = "Set the logging file string to include a lot of debugging information (SLOW!)",
              required = false)
    protected Boolean debugMode = false;

    /** this is used to indicate if they've asked for help */
    @Argument(fullName = "help", shortName = "h", doc = "Generate this help message", required = false)
    public Boolean help = false;

    /** our logging output patterns */
    private static String patternString = "%-5p %d{HH:mm:ss,SSS} %C{1} - %m %n";
    private static String debugPatternString = "%n[level] %p%n[date]\t\t %d{dd MMM yyyy HH:mm:ss,SSS} %n[class]\t\t %C %n[location]\t %l %n[line number]\t %L %n[message]\t %m %n";

    /**
     * Allows a given application to return a brief description of itself.
     *
     * @return An ApplicationDetails object describing the current application.  Should not be null.
     */
    protected ApplicationDetails getApplicationDetails() {
        return new ApplicationDetails(ApplicationDetails.createDefaultHeader(getClass()),
                                      ApplicationDetails.createDefaultRunningInstructions(getClass()));
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
     * this is the function that the inheriting class can expect to have called
     * when all the argument processing is done
     *
     * @return the return code to exit the program with
     */
    protected abstract int execute();

    static {
        // setup a basic log configuration
        BasicConfigurator.configure();
    }

    public static int result = -1;

    /**
     * This function is called to start processing the command line, and kick
     * off the execute message of the program.
     *
     * @param clp  the command line program to execute
     * @param args the command line arguments passed in
     */
    @SuppressWarnings("unchecked")
    public static void start(CommandLineProgram clp, String[] args) {

        try {
            // setup our log layout
            PatternLayout layout = new PatternLayout();

            // setup the parser
            ParsingEngine parser = clp.parser = new ParsingEngine(clp);
            parser.addArgumentSource(clp.getClass());

            // process the args
            if (clp.canAddArgumentsDynamically()) {
                // if the command-line program can toss in extra args, fetch them and reparse the arguments.
                parser.parse(args);

                // Allow invalid and missing required arguments to pass this validation step.
                //   - InvalidArgument in case these arguments are specified by plugins.
                //   - MissingRequiredArgument in case the user requested help.  Handle that later, once we've
                //                             determined the full complement of arguments.
                parser.validate(EnumSet.of(ParsingEngine.ValidationType.MissingRequiredArgument,
                                           ParsingEngine.ValidationType.InvalidArgument));
                parser.loadArgumentsIntoObject(clp);

                Class[] argumentSources = clp.getArgumentSources();
                for (Class argumentSource : argumentSources)
                    parser.addArgumentSource(clp.getArgumentSourceName(argumentSource), argumentSource);
                parser.parse(args);

                if (isHelpPresent(parser))
                    printHelpAndExit(clp, parser);

                parser.validate();
            } else {
                parser.parse(args);

                if (isHelpPresent(parser))
                    printHelpAndExit(clp, parser);

                parser.validate();
                parser.loadArgumentsIntoObject(clp);
            }

            // if we're in debug mode, set the mode up
            if (clp.debugMode) {
                //logger.info("Setting debug");
                layout.setConversionPattern(debugPatternString);
            } else {
                //logger.info("not Setting debug");
                layout.setConversionPattern(patternString);
            }
            // now set the layout of all the loggers to our layout
            Enumeration<Appender> en = logger.getAllAppenders();
            for (; en.hasMoreElements();) {
                Appender app = en.nextElement();
                app.setLayout(layout);
            }

            // if they set the mode to quiet
            if (clp.quietMode) {

                // the only appender we should have is stdout, the following meathod is
                // deprecated, but the standard remove all appenders doesn't seem to work
                // TODO: find the right function
                //Category root = Category.getRoot();
                //root.removeAllAppenders();
                //logger.removeAllAppenders();
            }

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

            // set the default logger level
            clp.setupLoggerLevel();

            // regardless of what happens next, generate the header information
            HelpFormatter.generateHeaderInformation(clp.getApplicationDetails(), args);

            // call the execute
            CommandLineProgram.result = clp.execute();

            // return the result
            //System.exit(result);     // todo -- is this safe -- why exit here?  I want to run the GATK like normal
        }
        catch (ArgumentException e) {
            clp.parser.printHelp(clp.getApplicationDetails());
            // Rethrow the exception to exit with an error.
            throw e;
        }
        catch (Exception e) {
            // we catch all exceptions here. if it makes it to this level, we're in trouble.  Let's bail!
            // TODO: what if the logger is the exception? hmm...
            logger.fatal("\n");
            toErrorLog(clp, e);
            throw new RuntimeException(e);
        }
    }

    /**
     * generate an error log
     * @param clp the command line program
     * @param e the exception
     */
    private static void toErrorLog(CommandLineProgram clp, Exception e) {
        File logFile = new File("GATK_Error.log");
        PrintStream stream;
        try {
            stream = new PrintStream(logFile);
        } catch (Exception e1) { // catch all the exceptions here, if we can't create the file, do the alternate path
            if ( e.getCause() != null ) logger.fatal("with cause: " + e.getCause());
            throw new RuntimeException(e);
        }
        clp.generateErrorLog(stream, e);
        stream.close();
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
     */
    private void setupLoggerLevel() {
        Level par;
        if (logging_level.toUpperCase().equals("DEBUG")) {
            par = Level.DEBUG;
        } else if (logging_level.toUpperCase().equals("ERROR")) {
            par = Level.ERROR;
        } else if (logging_level.toUpperCase().equals("FATAL")) {
            par = Level.FATAL;
        } else if (logging_level.toUpperCase().equals("INFO")) {
            par = Level.INFO;
        } else if (logging_level.toUpperCase().equals("WARN")) {
            par = Level.WARN;
        } else if (logging_level.toUpperCase().equals("OFF")) {
            par = Level.OFF;
        } else {
            // we don't understand the logging level, let's get out of here
            throw new ArgumentException("Unable to match: " + logging_level + " to a logging level, make sure it's a valid level (INFO, DEBUG, ERROR, FATAL, OFF)");
        }

        logger.setLevel(par);
    }

    /**
     * a function used to indicate an error occurred in the command line tool
     *
     * @param msg message to display
     */
    private static void printExitSystemMsg(final String msg) {
        System.out.printf("The following error has occurred:%n%n");
        System.out.printf("%s:%n%n", msg);
        System.out.printf("Please check your command line arguments for any typos or inconsistencies.%n");
        System.out.printf("Also, please review our documentation at:%n");
        System.out.printf("        http://www.broadinstitute.org/gsa/wiki %n%n");
        System.out.printf("To report bugs or to get help resolving undocumented issues, please contact us via our support site at:%n");
        System.out.printf("        http://getsatisfaction.com/gsa %n%n");
        System.out.printf("Please be sure to include the stack trace below when posting a message on the support site:%n");
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
     * used to indicate an error occured
     *
     * @param msg the message to display
     */
    public static void exitSystemWithError(final String msg) {
        printExitSystemMsg(msg);
        System.exit(1);
    }

    /**
     * used to indicate an error occured
     *
     * @param msg the message
     * @param e   the error
     */
    public static void exitSystemWithError(final String msg, Exception e) {
        System.out.printf("------------------------------------------------------------------------------------------%n");
        printExitSystemMsg(msg);
        System.out.printf("------------------------------------------------------------------------------------------%n");
        e.printStackTrace();
        System.out.printf("------------------------------------------------------------------------------------------%n");
        System.exit(1);
    }

    /**
     * used to indicate an error occured
     *
     * @param e the exception occured
     */
    public static void exitSystemWithError(Exception e) {
        exitSystemWithError(e.getMessage(), e);
    }

    /**
     * generate an error log, given the stream to write to and the execption that was generated
     *
     * @param stream the output stream
     * @param e      the exception
     */
    public void generateErrorLog(PrintStream stream, Exception e) {
        stream.println(e.getStackTrace().toString());
    }

    /**
     * A hack to ensure that numbers are always formatted in the US style.
     */
    protected static void forceJVMLocaleToUSEnglish() {
        Locale.setDefault(Locale.US);
    }
}
