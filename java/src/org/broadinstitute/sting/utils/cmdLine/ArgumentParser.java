package org.broadinstitute.sting.utils.cmdLine;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import org.broadinstitute.sting.utils.Pair;

/**
 * User: aaron
 * Date: Mar 19, 2009
 * Time: 6:54:15 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class ArgumentParser {

    // what program are we parsing for
    private String programName;

    // our m_options
    private ArrayList<String> m_option_names = new ArrayList<String>();

    // where we eventually want the values to land
    private HashMap<String, Pair<Object,Field>> m_storageLocations = new HashMap<String, Pair<Object,Field>>();

    // create Options object
    protected Options m_options = new Options();

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ArgumentParser.class);

    // the reference to the command line program to fill in
    Object prog;

    public ArgumentParser(String programName, Object prog) {
        this.programName = programName;
        this.prog = prog;
    }


    /**
     * print out the help information
     */
    public void printHelp() {
        // automatically generate the help statement
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(100,
                "java -Xmx4096m -jar dist/GenomeAnalysisTK.jar",
                "",
                m_options,
                "",
                true);
    }


    /**
     * addOptionalArg
     * <p/>
     * Adds an optional argument to check on the command line
     *
     * @param name        the name of the argument, the long name
     * @param letterform  the short form
     * @param description the description of the argument
     * @param fieldname   the field to set when we've parsed this option
     */
    public void addOptionalArg(String name, String letterform, String description, String fieldname) {

        // we always want the help option to be available
        Option opt = OptionBuilder.withLongOpt(name).withArgName(name)
                .hasArg()
                .withDescription(description)
                .create(letterform);

        // add it to the option
        AddToOptionStorage(name, letterform, fieldname, opt);


    }

    /**
     * Used locally to add to the options storage we have, for latter processing
     *
     * @param name       the name of the option
     * @param letterform it's short form
     * @param fieldname  what field it should be stuck into on the calling class
     * @param opt        the option
     */
    private void AddToOptionStorage(String name, String letterform, String fieldname, Option opt) {
        AddToOptionStorage(name, letterform, getField(prog, fieldname), opt);
    }

    /**
     * Used locally to add to the options storage we have, for latter processing
     *
     * @param name       the name of the option
     * @param letterform it's short form
     * @param field      what field it should be stuck into on the calling class
     * @param opt        the option
     */
    private void AddToOptionStorage(String name, String letterform, Pair<Object,Field> field, Option opt) {
        // add to the option list
        m_options.addOption(opt);

        // first check to see if we've already added an option with the same name
        if (m_option_names.contains(letterform)) {
            throw new IllegalArgumentException(letterform + " was already added as an option");
        }

        // add the object with it's name to the storage location
        m_storageLocations.put( name, field );

        // add to the list of m_options
        m_option_names.add(letterform);
    }

    private Pair<Object,Field> getField( Object obj, String fieldName ) {
        try {
            return new Pair<Object,Field>( obj, obj.getClass().getField(fieldName) );
        } catch (NoSuchFieldException e) {
            logger.fatal("Failed to find the field specified by the fieldname parameter.");
            throw new RuntimeException(e.getMessage());
        }
    }
    
    /**
     * addRequiredArg
     * <p/>
     * Adds a required argument to check on the command line
     *
     * @param name        the name of the argument, the long name
     * @param letterform  the short form
     * @param description the description of the argument
     * @param fieldname   what field it should be stuck into on the calling class
     */
    public void addRequiredArg(String name, String letterform, String description, String fieldname) {
        // we always want the help option to be available
        Option opt = OptionBuilder.isRequired()
                .withLongOpt(name)
                .withArgName(name)
                .hasArg()
                .withDescription("(Required Option) " + description)
                .create(letterform);

        // add it to the option
        AddToOptionStorage(name, letterform, fieldname, opt);

    }

    /**
     * addOptionalArg
     * <p/>
     * Adds an optional argument to check on the command line
     *
     * @param name        the name of the argument, the long name
     * @param letterform  the short form
     * @param description the description of the argument
     * @param fieldname   what field it should be stuck into on the calling class
     */
    public void addOptionalArgList(String name, String letterform, String description, String fieldname) {
        // we always want the help option to be available
        Option opt = OptionBuilder.withLongOpt(name).withArgName(name)
                .hasArgs()
                .withDescription(description)
                .create(letterform);
        // add it to the option
        AddToOptionStorage(name, letterform, fieldname, opt);
    }



    /**
     * addRequiredArg
     * <p/>
     * Adds a required argument to check on the command line
     *
     * @param name        the name of the argument, the long name
     * @param letterform  the short form
     * @param description the description of the argument
     * @param fieldname   what field it should be stuck into on the calling class
     */
    public void addRequiredArgList(String name, String letterform, String description, String fieldname) {

        // we always want the help option to be available
        Option opt = OptionBuilder.isRequired()
                .withLongOpt(name)
                .withArgName(name)
                .hasArgs()
                .withDescription("(Required Option) " + description)
                .create(letterform);
        // add it to the option
        AddToOptionStorage(name, letterform, fieldname, opt);

    }

    /**
     * addOptionalFlag
     * <p/>
     * Adds an optional argument to check on the command line
     *
     * @param name        the name of the argument, the long name
     * @param letterform  the short form
     * @param description the description of the argument
     * @param fieldname   what field it should be stuck into on the calling class
     */
    public void addOptionalFlag(String name, String letterform, String description, String fieldname) {

        // if they've passed a non-Boolean as a object, beat them
        try {
            if (!(prog.getClass().getField(fieldname).getType() == Boolean.class)) {
                throw new IllegalArgumentException("Fields to addOptionalFlag must be of type Boolean");
            }
        } catch (NoSuchFieldException e) {
            throw new IllegalArgumentException("Fields to addOptionalFlag must exist!");
        }
        
        // we always want the help option to be available
        Option opt = OptionBuilder.withLongOpt(name)
                .withDescription(description)
                .create(letterform);


        // add it to the option
        AddToOptionStorage(name, letterform, fieldname, opt);

    }


    /**
     * addRequiredFlag
     * <p/>
     * Adds a required argument to check on the command line
     *
     * @param name        the name of the argument, the long name
     * @param letterform  the short form
     * @param description the description of the argument
     * @param fieldname   what field it should be stuck into on the calling class
     */
    public void addRequiredFlag(String name, String letterform, String description, String fieldname) {

        // if they've passed a non-Boolean as a object, beat them
        try {
            if (!(prog.getClass().getField(fieldname).getType() == Boolean.class)) {
                throw new IllegalArgumentException("Fields to addRequiredlFlag must be of type Boolean");
            }
        } catch (NoSuchFieldException e) {
            throw new IllegalArgumentException("Fields to addRequiredlFlag must exist!");
        }

        // we always want the help option to be available
        Option opt = OptionBuilder.isRequired()
                .withLongOpt(name)
                .withDescription("(Required Flag) " + description)
                .create(letterform);

        // add it to the option
        AddToOptionStorage(name, letterform, fieldname, opt);
    }


    /**
     * This function is called to validate all the arguments to the program.
     * If a required Arg isn't found, we generate the help message, and
     * exit the program
     *
     * @param args the command line arguments we recieved
     */
    public void processArgs(String[] args, boolean allowUnrecognized) throws ParseException {
        OurPosixParser parser = new OurPosixParser();
        Collection<Option> opts = m_options.getOptions();

        try {
            parser.parse(m_options, args, !allowUnrecognized);
        }
        catch (UnrecognizedOptionException e) {
            // we don't care about unknown exceptions right now
            logger.warn(e.getMessage());
            if(!allowUnrecognized)
                throw e;
        }

        // Apache CLI can ignore unrecognized arguments with a boolean flag, but
        // you can't get to the unparsed args.  Override PosixParser with a class
        // that can reach in and extract the protected command line.
        // TODO: Holy crap this is wacky.  Find a cleaner way.
        CommandLine cmd = parser.getCmd();

        // logger.info("We have " + opts.size() + " options");
        for (Option opt : opts) {
            if (cmd.hasOption(opt.getOpt())) {
                if (opt.hasArg()) {
                    //logger.info("looking at " + m_storageLocations.get(opt.getLongOpt()));
                    Object obj = m_storageLocations.get(opt.getLongOpt()).first;
                    Field field = m_storageLocations.get(opt.getLongOpt()).second;

                    try {
                        field.set(obj, constructFromString(field, cmd.getOptionValue(opt.getOpt())));
                    } catch (IllegalAccessException e) {
                        logger.fatal("processArgs: cannot convert field " + field.toString());
                        throw new RuntimeException("processArgs: Failed conversion " + e.getMessage());
                    }
                } else {
                    Object obj = m_storageLocations.get(opt.getLongOpt()).first;
                    Field field = m_storageLocations.get(opt.getLongOpt()).second;

                    try {
                        //logger.fatal("about to parse field " + f.getName());
                        field.set(obj, new Boolean(true));
                    } catch (IllegalAccessException e) {
                        logger.fatal("processArgs: cannot convert field " + field.toString());
                        throw new RuntimeException("processArgs: Failed conversion " + e.getMessage());
                    }
                }
            }
        }
    }

    private class OurPosixParser extends PosixParser {
        public CommandLine getCmd() { return cmd; }
    }

    /**
     * Extract arguments stored in annotations from fields of a given class.
     * @param source
     */
    public void addArgumentSource( Object source ) {
        Field[] fields = source.getClass().getFields();
        for(Field field: fields) {
            Argument arg = field.getAnnotation(Argument.class);
            if(arg == null)
                continue;

            String fullName = (arg.fullName().length() != 0) ? arg.fullName() : field.getName().trim().toLowerCase();
            String shortName = (arg.shortName().length() != 0) ? arg.shortName() : fullName.substring(0,1);
            if(shortName.length() != 1)
                throw new IllegalArgumentException("Invalid short name: " + shortName);
            String description = arg.required() ? "(Required Flag) " + arg.doc() : arg.doc();

            // TODO: Handle flags, handle lists
            OptionBuilder ob = OptionBuilder.withLongOpt(fullName).withArgName(fullName).hasArg();
            if( arg.required() ) ob = ob.isRequired();
            if( description.length() != 0 ) ob = ob.withDescription( description );

            Option option = ob.create( shortName );

            AddToOptionStorage(fullName, shortName, new Pair<Object,Field>( source, field ), option );
        }
    }

    private Object constructFromString(Field f, String str) {
        Type type = f.getType();
        // lets go through the types we support
        if (type == Boolean.TYPE) {
            boolean b = false;
            if (str.toLowerCase().equals("true")) {
                b = true;
            }
            Boolean bool = new Boolean(b);
            return bool;
        } else if (type == Integer.TYPE) {
            Integer in = Integer.valueOf(str);
            return in;
        } else if (type == Float.TYPE) {
            Float fl = Float.valueOf(str);
            return fl;
        } else {
            Constructor ctor = null;
            try {
                ctor = f.getType().getConstructor(String.class);
                return ctor.newInstance(str);
            } catch (NoSuchMethodException e) {
                logger.fatal("constructFromString:NoSuchMethodException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:NoSuchMethodException: Failed conversion " + e.getMessage());
            } catch (IllegalAccessException e) {
                logger.fatal("constructFromString:IllegalAccessException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            } catch (InvocationTargetException e) {
                logger.fatal("constructFromString:InvocationTargetException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:InvocationTargetException: Failed conversion " + e.getMessage());
            } catch (InstantiationException e) {
                logger.fatal("constructFromString:InstantiationException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }

        }
    }


}


/**

 public static void main(String[] args) {
 ArgumentParser p = new ArgumentParser("CrapApp");
 p.setupDefaultArgs();
 p.addRequiredArg("Flag","F","a required arg");
 p.addRequiredFlag("Sub","S","a required flag");
 p.addOptionalArg("Boat","T","Maybe you want a boat?");
 String[] str = {"--Flag","rrr","-T","ppp", "--Flag","ttt"};
 p.processArgs(str);
 Iterator<String> r = map.keySet().iterator();
 while (r.hasNext()) {
 String key = r.next();
 String[] q = map.get(key);
 if (q != null) {
 for (String mystr : q) {
 System.err.println("key: " + key + " val: " + mystr);
 }
 }
 }
 }
 */
