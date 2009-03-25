package org.broadinstitute.sting.utils.cmdLine;

import org.apache.log4j.*;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

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
     * Our Argument parser, which handles parsing the command line in GNU format
     */
    protected ArgumentParser m_parser;

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    private static Logger logger = Logger.getRootLogger();// .getLogger(CommandLineProgram.class);

    /**
     * the default log level
     */
    public String logging_level = "ERROR";


    /**
     * where to send the output of our logger
     */
    public String toFile = null;

    /**
     * do we want to silence the command line output
     */
    public Boolean quietMode = false;

    /**
     * do we want to generate debugging information with the logs
     */
    public Boolean debugMode = false;


    /**
     * our logging output patterns
     */
    private static String patternString = "%p %m %n";
    private static String debugPatternString = "%n[level] %p%n[date]\t\t %d{dd MMM yyyy HH:mm:ss,SSS} %n[class]\t\t %C %n[location]\t %l %n[line number]\t %L %n[message]\t %m %n";

    /**
     * the contract for the inheriting class is that they have a setupArgs()
     * function which sets up the args to the specific program.
     */
    protected abstract void setupArgs();

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
    protected Object[] getArgumentSources() { return new Object[] {}; }

    /**
     * this is the function that the inheriting class can expect to have called
     * when all the argument processing is done
     * 
     * @return the return code to exit the program with
     */
    protected abstract int execute();


    /**
     * this is used to indicate if they've asked for help
     */
    public Boolean help = false;


    /**
     * This function is called to start processing the command line, and kick
     * off the execute message of the program.
     *
     * @param clp the command line program to execute
     * @param args the command line arguments passed in
     */
    public static void start(CommandLineProgram clp, String[] args) {

        try {
            // setup a basic log configuration
            BasicConfigurator.configure();

            // setup our log layout
            PatternLayout layout = new PatternLayout();


            // setup the parser
            clp.m_parser = new ArgumentParser(clp.getClass().getName(), clp);

            // setup the default help and logging args controlled by the base class
            clp.setupDefaultArgs();

            // setup the args
            clp.setupArgs();

            // process the args
            if( clp.canAddArgumentsDynamically() ) {
                // if the command-line program can toss in extra args, fetch them and reparse the arguments.
                clp.m_parser.processArgs(args, true);
                Object[] argumentSources = clp.getArgumentSources();
                for( Object argumentSource: argumentSources )
                    clp.addArgumentSource( argumentSource );
                clp.m_parser.processArgs(args, false);
            }
            else {
                clp.m_parser.processArgs(args, false);
            }

            // if we're in debug mode, set the mode up
            if (clp.debugMode) {
                //logger.info("Setting debug");
                layout.setConversionPattern(debugPatternString);
            } else {
                //logger.info("not Setting debug");
                layout.setConversionPattern(patternString);
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

            // they asked for help, give it to them
            if (clp.help) {
                clp.m_parser.printHelp();
                System.exit(1);
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
        catch (org.apache.commons.cli.ParseException e) {
            logger.fatal("Unable to pass command line arguments: " + e.getMessage() );
            clp.m_parser.printHelp();
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
        logger.info("Program Name: " + clp.getClass().getName());

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
    private void setupLoggerLevel() {

        Level par = Level.ERROR;
        if (logging_level.equals("DEBUG")) {
            par = Level.DEBUG;
        }
        if (logging_level.equals("ERROR")) {
            par = Level.ERROR;
        }
        if (logging_level.equals("FATAL")) {
            par = Level.FATAL;
        }
        if (logging_level.equals("INFO")) {
            par = Level.INFO;
        }
        if (logging_level.equals("WARN")) {
            par = Level.WARN;
        }
        if (logging_level.equals("OFF")) {
            par = Level.OFF;
        }

        logger.setLevel(par);
    }

    /**
     * Pass along a new set of valid command line arguments.  In this case,
     * probably a class with @argument or @flag annotations.
     * @param source
     */
    private void addArgumentSource( Object source ) {
        m_parser.addArgumentSource(source);
    }

    /**
     * we have some default options that should always get checked for in the
     * arguments provided to the program
     */
    private void setupDefaultArgs() {
        m_parser.addOptionalFlag("help", "h", "Generate this help message", "help");
        m_parser.addOptionalArg("logging_level", "l", "Set the logging level", "logging_level");
        m_parser.addOptionalArg("log_to_file", "log", "Set the logging location", "toFile");
        m_parser.addOptionalFlag("quiet_output_mode", "quiet", "Set the logging to quiet mode, no output to stdout", "quietMode");
        m_parser.addOptionalFlag("debug_mode", "debug", "Set the logging file string to include a lot of debugging information (SLOW!)", "debugMode");
    }


}
