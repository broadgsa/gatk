package org.broadinstitute.sting.utils.cmdLine;

import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Modifier;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.HashSet;
import java.util.Set;

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
        AddToOptionStorage(opt, fieldname);


    }

    /**
     * Used locally to add to the options storage we have, for latter processing
     *
     * @param opt        the option
     * @param fieldname  what field it should be stuck into on the calling class
     */
    private void AddToOptionStorage(Option opt, String fieldname) {
        AddToOptionStorage( opt, getField(prog, fieldname) );
    }

    /**
     * Used locally to add to the options storage we have, for latter processing
     *
     * @param opt        the option
     * @param field      what field it should be stuck into on the calling class
     */
    private void AddToOptionStorage(Option opt, Pair<Object,Field> field ) {
        // first check to see if we've already added an option with the same name
        if (m_options.hasOption( opt.getOpt() ))
            throw new IllegalArgumentException(opt.getOpt() + " was already added as an option");

        // add to the option list
        m_options.addOption(opt);

        // add the object with it's name to the storage location
        m_storageLocations.put( opt.getLongOpt(), field );
    }

    /**
     * Used locally to add a group of mutually exclusive options to options storage.
     * @param options A list of pairs of param, field to add.
     */
    private void AddToOptionStorage( List<Pair<Option,Pair<Object,Field>>> options ) {
        OptionGroup optionGroup = new OptionGroup();
        boolean isRequired = true;

        for( Pair<Option,Pair<Object,Field>> option: options ) {
            if (m_options.hasOption(option.first.getOpt()) )
                throw new IllegalArgumentException(option.first.getOpt() + " was already added as an option");

            optionGroup.addOption(option.first);
            m_storageLocations.put( option.first.getLongOpt(), option.second );
            isRequired &= option.first.isRequired();
        }

        optionGroup.setRequired(isRequired);
        m_options.addOptionGroup(optionGroup);
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
        AddToOptionStorage( opt, fieldname );

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
        AddToOptionStorage( opt, fieldname );
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
        AddToOptionStorage(opt, fieldname);

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
        AddToOptionStorage( opt, fieldname );

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
        AddToOptionStorage( opt, fieldname );
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
            if(!allowUnrecognized) {
                logger.warn(e.getMessage());
                throw e;
            }
        }

        // Apache CLI can ignore unrecognized arguments with a boolean flag, but
        // you can't get to the unparsed args.  Override PosixParser with a class
        // that can reach in and extract the protected command line.
        // TODO: Holy crap this is wacky.  Find a cleaner way.
        CommandLine cmd = parser.getCmd();

        // logger.info("We have " + opts.size() + " options");
        for (Option opt : opts) {
            if (cmd.hasOption(opt.getOpt())) {
                //logger.info("looking at " + m_storageLocations.get(opt.getLongOpt()));
                Object obj = m_storageLocations.get(opt.getLongOpt()).first;
                Field field = m_storageLocations.get(opt.getLongOpt()).second;

                try {
                    if (opt.hasArg())
                        field.set(obj, constructFromString(field, cmd.getOptionValues(opt.getOpt())));
                    else
                        field.set(obj, new Boolean(true));
                } catch (IllegalAccessException e) {
                    logger.fatal("processArgs: cannot convert field " + field.toString());
                    throw new RuntimeException("processArgs: Failed conversion " + e.getMessage());
                }
            }
        }
    }

    /**
     * Simple class to wrap the posix parser and get access to the parsed cmd object.
     */
    private class OurPosixParser extends PosixParser {
        public CommandLine getCmd() { return cmd; }
    }

    /**
     * Extract arguments stored in annotations from fields of a given class.
     * @param source Source of arguments, probably provided through Argument annotation.
     */
    public void addArgumentSource( Object source ) {
        Field[] fields = source.getClass().getFields();

        for( Set<Field> optionGroup: groupExclusiveOptions(fields) ) {
            List<Pair<Option,Pair<Object,Field>>> options = new ArrayList<Pair<Option,Pair<Object,Field>>>();
            for( Field field: optionGroup ) {
                Argument argument = field.getAnnotation(Argument.class);
                Option option = createOptionFromField( source, field, argument );
                options.add( new Pair<Option,Pair<Object,Field>>( option, new Pair<Object,Field>( source,field) ) );
            }

            if( options.size() == 1 )
                AddToOptionStorage( options.get(0).first, new Pair<Object,Field>( source, options.get(0).second.second ) );
            else {
                AddToOptionStorage( options );
            }
        }
    }

    /**
     * Group mutually exclusive options together into sets.  Non-exclusive options
     * will be alone a set.  Every option should appear in exactly one set.
     * WARNING: Has no concept of nested dependencies.
     * @param fields list of fields for which to check options.
     * @return
     */
    private List<Set<Field>> groupExclusiveOptions( Field[] fields ) {
        List<Set<Field>> optionGroups = new ArrayList<Set<Field>>();
        for(Field field: fields) {
            Argument argument = field.getAnnotation(Argument.class);
            if(argument == null)
                continue;

            String[] exclusives = argument.exclusive().split(",");
            if( exclusives.length != 0 ) {
                HashSet<Field> matchingFields = new HashSet<Field>();

                // Find the set of all options exclusive to this one.
                matchingFields.add( field );
                for( Field candidate: fields ) {
                    for( String exclusive: exclusives ) {
                        if( candidate.getName().equals( exclusive.trim() ) )
                            matchingFields.add( candidate );
                    }
                }

                // Perform a set intersection to see whether this set of fields intersects
                // with any existing set.
                //
                // If so, add any additional elements to the list.
                boolean setExists = false;
                for( Set<Field> optionGroup: optionGroups ) {
                    Set<Field> working = new HashSet<Field>(optionGroup);
                    working.retainAll( matchingFields );
                    if( working.size() > 0 ) {
                        optionGroup.addAll( matchingFields );
                        setExists = true;
                    }
                }

                // Otherwise, add a new option group.
                if( !setExists )
                    optionGroups.add( matchingFields );
            }
        }

        return optionGroups;
    }

    /**
     * Given a field with some annotations, create a command-line option.  If not enough data is
     * available to create a command-line option, return null.
     * @param source Source class containing the field.
     * @param field Field
     * @return Option representing the field options.
     */
    private Option createOptionFromField( Object source, Field field, Argument argument ) {

        String fullName = (argument.fullName().length() != 0) ? argument.fullName() : field.getName().trim().toLowerCase();
        String shortName = (argument.shortName().length() != 0) ? argument.shortName() : fullName.substring(0,1);
        if(shortName.length() != 1)
            throw new IllegalArgumentException("Invalid short name: " + shortName);
        boolean isFlag = (field.getType() == Boolean.class) || (field.getType() == Boolean.TYPE);
        boolean isCollection = field.getType().isArray() || Collection.class.isAssignableFrom(field.getType());

        if( isFlag && isCollection )
            throw new IllegalArgumentException("Can't have an array of flags.");

        String description = String.format("[%s] %s", source, argument.doc());

        OptionBuilder ob = OptionBuilder.withLongOpt(fullName);
        if( !isFlag ) {
            ob = ob.withArgName(fullName);
            ob = isCollection ? ob.hasArgs() : ob.hasArg();
            if( argument.required() ) {
                ob = ob.isRequired();
                description = String.format("[%s] (Required Option) %s", source, argument.doc());
            }
        }
        if( description.length() != 0 ) ob = ob.withDescription( description );

        Option option = ob.create( shortName );

        return option;
    }

    private Object constructFromString(Field f, String[] strs) {
        Class type = f.getType();

        if( Collection.class.isAssignableFrom(type) ) {
            Collection collection = null;
            Class containedType = null;

            // If this is a parameterized collection, find the contained type.  If blow up if only one type exists.
            if( f.getGenericType() instanceof ParameterizedType ) {
                ParameterizedType parameterizedType = (ParameterizedType)f.getGenericType();
                if( parameterizedType.getActualTypeArguments().length > 1 )
                    throw new IllegalArgumentException("Unable to determine collection type of field: " + f.toString());
                containedType = (Class)parameterizedType.getActualTypeArguments()[0];
            }
            else
                containedType = String.class;

            // If this is a generic interface, pick a concrete implementation to create and pass back.
            // Because of type erasure, don't worry about creating one of exactly the correct type.
            if( Modifier.isInterface(type.getModifiers()) || Modifier.isAbstract(type.getModifiers()) )
            {
                if( java.util.List.class.isAssignableFrom(type) ) type = ArrayList.class;
                else if( java.util.Queue.class.isAssignableFrom(type) ) type = java.util.ArrayDeque.class;
                else if( java.util.Set.class.isAssignableFrom(type) ) type = java.util.TreeSet.class;
            }

            try
            {
                collection = (Collection)type.newInstance();
            }
            catch( Exception ex ) {
                // Runtime exceptions are definitely unexpected parsing simple collection classes.
                throw new IllegalArgumentException(ex);
            }

            for( String str: strs )
                collection.add( constructSingleElement(f,containedType,str) );

            return collection;
        }
        else if( type.isArray() ) {
            Class containedType = type.getComponentType();

            Object arr = Array.newInstance(containedType,strs.length);
            for( int i = 0; i < strs.length; i++ )
                Array.set( arr,i,constructSingleElement(f,containedType,strs[i]) );
            return arr;
        }
        else  {
            if( strs.length != 1 )
                throw new IllegalArgumentException("Passed multiple arguments to an object expecting a single value.");
            return constructSingleElement(f,type,strs[0]);
        }
    }

    /**
     * Builds a single element of the given type.
     * @param f Implies type of data to construct.
     * @param str String representation of data.
     * @return parsed form of String.
     */
    private Object constructSingleElement(Field f, Class type, String str) {
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
        }
        else {
            Constructor ctor = null;
            try {
                ctor = type.getConstructor(String.class);
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
