package org.broadinstitute.sting.utils.cmdLine;

import org.apache.log4j.*;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.EnumSet;
import java.util.Enumeration;

/**
 * User: aaron
 * Date: Mar 19, 2009
 * Time: 3:54:56 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 * <p/>
 * <p/>
 * This class is our implementation of the command line parser, similar to Pickard's.  We instead
 * support GNU style command line arguements, and use this class to setup the global parser.
 */
public abstract class CommandLineProgram {

    /**
     * The command-line program and the arguments it returned.
     */
    private ParsingEngine parser = null;

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    private static Logger logger = Logger.getRootLogger();// .getLogger(CommandLineProgram.class);

    /**
     * the default log level
     */
    @Argument(fullName="logging_level",
              shortName="l",
              doc="Set the minimum level of logging, i.e. setting INFO get's you INFO up to FATAL, setting ERROR gets you ERROR and FATAL level logging. (DEBUG, INFO, WARN, ERROR, FATAL, OFF). ",
              required=false)    
    protected String logging_level = "ERROR";


    /**
     * where to send the output of our logger
     */
    @Argument(fullName="log_to_file",
              shortName="log",
              doc="Set the logging location",
              required=false)    
    protected String toFile = null;

    /**
     * do we want to silence the command line output
     */
    @Argument(fullName="quiet_output_mode",
              shortName="quiet",
              doc="Set the logging to quiet mode, no output to stdout",
              required=false)
    protected Boolean quietMode = false;

    /**
     * do we want to generate debugging information with the logs
     */
    @Argument(fullName="debug_mode",
              shortName="debug",
              doc="Set the logging file string to include a lot of debugging information (SLOW!)",
              required=false)    
    protected Boolean debugMode = false;

    /**
     * this is used to indicate if they've asked for help
     */
    @Argument(fullName="help",shortName="h",doc="Generate this help message",required=false)
    public Boolean help = false;    

    /**
     * our logging output patterns
     */
    private static String patternString = "%-5p %d{HH:mm:ss,SSS} %C{1} - %m %n";
    private static String debugPatternString = "%n[level] %p%n[date]\t\t %d{dd MMM yyyy HH:mm:ss,SSS} %n[class]\t\t %C %n[location]\t %l %n[line number]\t %L %n[message]\t %m %n";

    /**
     * Allows a given application to return a brief description of itself.
     * @return An ApplicationDetails object describing the current application.  Should not be null. 
     */
    protected ApplicationDetails getApplicationDetails() {
        return new ApplicationDetails( ApplicationDetails.createDefaultHeader(getClass()),
                                       ApplicationDetails.createDefaultRunningInstructions(getClass()) );
    }

    /**
     * Will this application want to vary its argument list dynamically?
     * If so, parse the command-line options and then prompt the subclass to return
     * a list of argument providers.
     * @return Whether the application should vary command-line arguments dynamically.
     */
    protected boolean canAddArgumentsDynamically() { return false; }

    /**
     * Provide a list of object to inspect, looking for additional command-line arguments.
     * @return A list of objects to inspect.
     */
    protected Class[] getArgumentSources() { return new Class[] {}; }

    /**
     * Allows arguments to be hijacked by subclasses of the program before being placed
     * into plugin classes.
     * @param source Source class for the argument.
     * @param targetInstance Instance into which the value should be ultimately injected.
     * @param value Value to inject.
     * @return True if the particular field has been hijacked; false otherwise.
     */
    protected boolean intercept( ArgumentSource source, Object targetInstance, Object value ) { return false; }

    /**
     * Name this argument source.  Provides the (full) class name as a default.
     * @param source The argument source.
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

    /**
     * This function is called to start processing the command line, and kick
     * off the execute message of the program.
     *
     * @param clp the command line program to execute
     * @param args the command line arguments passed in
     */
    public static void start(CommandLineProgram clp, String[] args) {

        try {
            // setup our log layout
            PatternLayout layout = new PatternLayout();

            // setup the parser
            ParsingEngine parser = clp.parser = new ParsingEngine(clp);
            parser.addArgumentSource( clp.getClass() );

            // process the args
            if( clp.canAddArgumentsDynamically() ) {
                // if the command-line program can toss in extra args, fetch them and reparse the arguments.
                parser.parse(args);

                // Allow invalid and missing required arguments to pass this validation step.
                //   - InvalidArgument in case these arguments are specified by plugins.
                //   - MissingRequiredArgument in case the user requested help.  Handle that later, once we've
                //                             determined the full complement of arguments.
                parser.validate( EnumSet.of(ParsingEngine.ValidationType.MissingRequiredArgument,
                                            ParsingEngine.ValidationType.InvalidArgument) );
                parser.loadArgumentsIntoObject( clp );

                Class[] argumentSources = clp.getArgumentSources();
                for( Class argumentSource: argumentSources )
                    parser.addArgumentSource( clp.getArgumentSourceName(argumentSource), argumentSource );
                parser.parse(args);

                if( isHelpPresent( clp, parser ) )
                    printHelpAndExit( clp, parser );                

                parser.validate();
            }
            else {
                parser.parse(args);

                if( isHelpPresent( clp, parser ) )
                    printHelpAndExit( clp, parser );

                parser.validate();
                parser.loadArgumentsIntoObject( clp );
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
                FileAppender appender = null;
                try {
                    appender = new FileAppender(layout, clp.toFile, false);
                    logger.addAppender(appender);
                } catch (IOException e) {
                    throw new RuntimeException("Unable to re-route log output to " + clp.toFile + " make sure the destination exists");
                }
            }

            // regardless of what happens next, generate the header information
            generateHeaderInformation(clp, args);

            // set the default logger level
            clp.setupLoggerLevel();

            // call the execute
            int result = clp.execute();

            // return the result
            System.exit(result);
        }
        catch (ArgumentException e) {
            clp.parser.printHelp( clp.getApplicationDetails() );
            // Rethrow the exception to exit with an error.
            throw e;
        }
        catch (Exception e) {
            // we catch all exceptions here. if it makes it to this level, we're in trouble.  Let's bail!
            // TODO: what if the logger is the exception? hmm...
            logger.fatal("Exception caught by base Command Line Program, with message: " + e.getMessage());
            logger.fatal("with cause: " + e.getCause());
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    /**
     * Find fields in the object obj that look like command-line arguments, and put command-line
     * arguments into them.
     * @param obj Object to inspect for command line arguments.
     */
    public void loadArgumentsIntoObject( Object obj ) {
        parser.loadArgumentsIntoObject( obj );
    }

    /**
     * a manual way to load argument providing objects into the program
     * @param clp the command line program
     * @param cls the class to load the arguments off of
     */
    public void loadAdditionalSource(CommandLineProgram clp, Class cls ) {
        parser.addArgumentSource( clp.getArgumentSourceName(cls), cls );
    }

    /**
     * generateHeaderInformation
     * <p/>
     * 
     * Generate a standard header for the logger
     * @param clp the command line program to execute
     * @param args the command line arguments passed in
     *
     **/
    protected static void generateHeaderInformation(CommandLineProgram clp, String[] args) {

        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        java.util.Date date = new java.util.Date();

        logger.info("-------------------------------------------------------");
        for( String headerLine: clp.getApplicationDetails().applicationHeader )
            logger.info(headerLine);
        String output = "";
        for (String str : args) {
            output = output + str + " ";
        }
        logger.info("Program Args: " + output);
        logger.info("Time/Date: " + dateFormat.format(date));
        logger.info("-------------------------------------------------------");
    }

    /**
     * this function checks the logger level passed in on the command line, taking the lowest
     * level that was provided.
     */
    private void setupLoggerLevel()  {

        Level par = Level.ERROR;
        if (logging_level.equals("DEBUG")) {
            par = Level.DEBUG;
        }
        else if (logging_level.equals("ERROR")) {
            par = Level.ERROR;
        }
        else if (logging_level.equals("FATAL")) {
            par = Level.FATAL;
        }
        else if (logging_level.equals("INFO")) {
            par = Level.INFO;
        }
        else if (logging_level.equals("WARN")) {
            par = Level.WARN;
        }
        else if (logging_level.equals("OFF")) {
            par = Level.OFF;
        }
        else {
            // we don't understand the logging level, let's get out of here
            throw new ArgumentException("Unable to match: " + logging_level + " to a logging level, make sure it's a valid level (INFO, DEBUG, ERROR, FATAL, OFF)");
        }

        logger.setLevel(par);
    }

    /**
     * a function used to indicate an error occured in the command line tool
     *
     * @param msg
     */
    private static void printExitSystemMsg(final String msg) {
        System.out.printf("------------------------------------------------------------------------------------------%n");
        System.out.printf("An error has occurred.  Please check your command line arguments for any typos or inconsistencies.%n%n");
        System.out.printf("For assistance, please email us at gsadevelopers@broad.mit.edu, or review our documentation at http://www.broadinstitute.org/gsa/wiki.%n");
    }

    /**
     * Do a cursory search for the given argument.
     * @param clp Instance of the command-line program.
     * @param parser Parser
     * @return True if help is present; false otherwise.
     */
    private static boolean isHelpPresent( CommandLineProgram clp, ParsingEngine parser ) {
        return parser.isArgumentPresent("help");
    }

    /**
     * Print help and exit.
     * @param clp Instance of the command-line program.
     * @param parser True if help is present; false otherwise.
     */
    private static void printHelpAndExit( CommandLineProgram clp, ParsingEngine parser ) {
        parser.printHelp( clp.getApplicationDetails() );
        System.exit(0);
    }

    /**
     * used to indicate an error occured
     * @param msg the message to display
     */
    public static void exitSystemWithError(final String msg) {
        printExitSystemMsg(msg);
        System.exit(1);
    }

    /**
     * used to indicate an error occured
     * @param msg the message
     * @param e the error
     */
    public static void exitSystemWithError(final String msg, Exception e) {
        e.printStackTrace();
        printExitSystemMsg(msg);
        System.exit(1);
    }

    /**
     * used to indicate an error occured
     * @param e the exception occured
     */
    public static void exitSystemWithError(Exception e) {
        exitSystemWithError(e.getMessage(), e);
    }
}
