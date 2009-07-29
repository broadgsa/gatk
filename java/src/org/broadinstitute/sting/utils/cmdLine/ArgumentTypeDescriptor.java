/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.sam.SAMFileWriterBuilder;
import org.broadinstitute.sting.utils.sam.SAMFileReaderBuilder;
import org.apache.log4j.Logger;

import java.lang.reflect.*;
import java.util.*;
import java.io.File;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileReader;

/**
 * An factory capable of providing parsers that can parse any type
 * of supported command-line argument.
 *
 * @author mhanna
 * @version 0.1
 */
public abstract class ArgumentTypeDescriptor {

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ArgumentTypeDescriptor.class);
    
    /**
     * Class reference to the different types of descriptors that the create method can create.
     */
    private static List<ArgumentTypeDescriptor> descriptors = Arrays.asList( new SAMFileReaderArgumentTypeDescriptor(),
                                                                             new SAMFileWriterArgumentTypeDescriptor(),
                                                                             new SimpleArgumentTypeDescriptor(),                                                                                
                                                                             new CompoundArgumentTypeDescriptor() );

    public static ArgumentTypeDescriptor create( Class type ) {
        for( ArgumentTypeDescriptor descriptor: descriptors ) {
            if( descriptor.supports(type) )
                return descriptor;
        }
        throw new StingException("Can't process command-line arguments of type: " + type.getName());
    }

    /**
     * Does this descriptor support classes of the given type?
     * @param type The type to check.
     * @return true if this descriptor supports the given type, false otherwise.
     */
    public abstract boolean supports( Class type );

    /**
     * Given the given argument source and attributes, synthesize argument definitions for command-line arguments.
     * @param source Source class and field for the given argument.
     * @param description Description of the fields that go into a given argument.
     * @return A list of command-line argument definitions supporting this field.
     */
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source, Argument description ) {
        ArgumentDefinition definition =  new ArgumentDefinition( source,
                                                                 getFullName( source, description ),
                                                                 getShortName( source, description ),
                                                                 getDoc( source, description ),
                                                                 isRequired( source, description ),
                                                                 getExclusiveOf( source, description ),
                                                                 getValidationRegex( source, description ) );
        return Collections.singletonList(definition);
    }
    
    public Object parse( ArgumentSource source, ArgumentMatch... values ) {
        return parse( source, source.field.getType(), values ); 
    }

    protected abstract Object parse( ArgumentSource source, Class type, ArgumentMatch... values );

    /**
     * Retrieves the full name of the argument, specifiable with the '--' prefix.  The full name can be
     * either specified explicitly with the fullName annotation parameter or implied by the field name.
     * @return full name of the argument.  Never null.
     */
    protected String getFullName( ArgumentSource source, Argument description ) {
        return description.fullName().trim().length() > 0 ? description.fullName().trim() : source.field.getName().toLowerCase();
    }

    /**
     * Retrieves the short name of the argument, specifiable with the '-' prefix.  The short name can
     * be specified or not; if left unspecified, no short name will be present.
     * @return short name of the argument.  Null if no short name exists.
     */
    protected String getShortName( ArgumentSource source, Argument description ) {
        return description.shortName().trim().length() > 0 ? description.shortName().trim() : null;
    }

    /**
     * Documentation for this argument.  Mandatory field.
     * @return Documentation for this argument.
     */
    protected String getDoc( ArgumentSource source, Argument description ) {
        return description.doc();
    }

    /**
     * Returns whether this field is required.  Note that flag fields are always forced to 'not required'.
     * @return True if the field is mandatory and not a boolean flag.  False otherwise.
     */
    protected boolean isRequired( ArgumentSource source, Argument description ) {
        return description.required() && !source.isFlag();
    }

    /**
     * Specifies other arguments which cannot be used in conjunction with tihs argument.  Comma-separated list.
     * @return A comma-separated list of exclusive arguments, or null if none are present.
     */
    protected String getExclusiveOf( ArgumentSource source, Argument description ) {
        return description.exclusiveOf().trim().length() > 0 ? description.exclusiveOf().trim() : null;
    }

    /**
     * A regular expression which can be used for validation.
     * @return a JVM regex-compatible regular expression, or null to permit any possible value.
     */
    protected String getValidationRegex( ArgumentSource source, Argument description ) {
        return description.validation().trim().length() > 0 ? description.validation().trim() : null;
    }
}

/**
 * Parse simple argument types: java primitives, wrapper classes, and anything that has
 * a simple String constructor.
 */
class SimpleArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    @Override
    public boolean supports( Class type ) {
        if( type.isPrimitive() ) return true;
        if( type.isEnum() ) return true;
        if( primitiveToWrapperMap.containsValue(type) ) return true;

        try {
            type.getConstructor(String.class);
            return true;
        }
        catch( Exception ex ) {
            // An exception thrown above means that the String constructor either doesn't
            // exist or can't be accessed.  In either case, this descriptor doesn't support this type.
            return false;
        }
    }

    @Override
    protected Object parse( ArgumentSource source, Class type, ArgumentMatch... matches ) {
        if( matches.length > 1 || matches[0].values().size() > 1 )
            throw new StingException("Simple argument parser is unable to parse multiple arguments.");
        String value = matches[0].values().get(0);

        // lets go through the types we support
        try {
            if (type.isPrimitive()) {
                Method valueOf = primitiveToWrapperMap.get(type).getMethod("valueOf",String.class);
                return valueOf.invoke(null,value.trim());
            } else if (type.isEnum()) {
                return Enum.valueOf(type,value.toUpperCase().trim());
            } else {
                Constructor ctor = type.getConstructor(String.class);
                return ctor.newInstance(value);
            }
        }
        catch (NoSuchMethodException e) {
            throw new StingException("constructFromString:NoSuchMethodException: Failed conversion " + e.getMessage());
        } catch (IllegalAccessException e) {
            throw new StingException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
        } catch (InvocationTargetException e) {
            throw new StingException("constructFromString:InvocationTargetException: Failed conversion " + e.getMessage());
        } catch (InstantiationException e) {
            throw new StingException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
        }

    }

    /**
     * A mapping of the primitive types to their associated wrapper classes.  Is there really no way to infer
     * this association available in the JRE?
     */
    private static Map<Class,Class> primitiveToWrapperMap = new HashMap<Class,Class>() {
        {
            put( Boolean.TYPE, Boolean.class );
            put( Character.TYPE, Character.class );
            put( Byte.TYPE, Byte.class );
            put( Short.TYPE, Short.class );
            put( Integer.TYPE, Integer.class );
            put( Long.TYPE, Long.class );
            put( Float.TYPE, Float.class );
            put( Double.TYPE, Double.class );
        }
    };    
}

/**
 * Process compound argument types: arrays, and typed and untyped collections.
 */
class CompoundArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    @Override
    public boolean supports( Class type ) {
        return ( Collection.class.isAssignableFrom(type) || type.isArray() );
    }
    
    @Override
    public Object parse( ArgumentSource source, Class type, ArgumentMatch... matches )
    {
        Class componentType = null;

        if( matches.length > 1 )
            throw new StingException("Simple argument parser is unable to combine multiple argument types into a compound argument.");
        ArgumentMatch match = matches[0];

        if( Collection.class.isAssignableFrom(type) ) {

            // If this is a generic interface, pick a concrete implementation to create and pass back.
            // Because of type erasure, don't worry about creating one of exactly the correct type.
            if( Modifier.isInterface(type.getModifiers()) || Modifier.isAbstract(type.getModifiers()) )
            {
                if( java.util.List.class.isAssignableFrom(type) ) type = ArrayList.class;
                else if( java.util.Queue.class.isAssignableFrom(type) ) type = java.util.ArrayDeque.class;
                else if( java.util.Set.class.isAssignableFrom(type) ) type = java.util.TreeSet.class;
            }

            // If this is a parameterized collection, find the contained type.  If blow up if only one type exists.
            if( source.field.getGenericType() instanceof ParameterizedType) {
                ParameterizedType parameterizedType = (ParameterizedType)source.field.getGenericType();
                if( parameterizedType.getActualTypeArguments().length > 1 )
                    throw new IllegalArgumentException("Unable to determine collection type of field: " + source.field.toString());
                componentType = (Class)parameterizedType.getActualTypeArguments()[0];
            }
            else
                componentType = String.class;

            ArgumentTypeDescriptor componentArgumentParser = ArgumentTypeDescriptor.create( componentType );

            Collection collection = null;
            try {
                collection = (Collection)type.newInstance();
            }
            catch (InstantiationException e) {
                logger.fatal("ArgumentParser: InstantiationException: cannot convert field " + source.field.getName());
                throw new StingException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }
            catch (IllegalAccessException e) {
                logger.fatal("ArgumentParser: IllegalAccessException: cannot convert field " + source.field.getName());
                throw new StingException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            }

            for( ArgumentMatch value: match )
                collection.add( componentArgumentParser.parse(source,componentType,value) );

            return collection;

        }
        else if( type.isArray() ) {
            componentType = type.getComponentType();
            ArgumentTypeDescriptor componentArgumentParser = ArgumentTypeDescriptor.create( componentType );
            Object arr = Array.newInstance(componentType,match.values().size());

            int i = 0;
            for( ArgumentMatch value: match )
                Array.set( arr,i++,componentArgumentParser.parse(source,componentType,value));

            return arr;
        }
        else
            throw new StingException("Unsupported compound argument type: " + type);
    }
}

/**
 * Handle SAMFileReaders.
 */
class SAMFileReaderArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    @Override
    public boolean supports( Class type ) {
        return SAMFileReader.class.isAssignableFrom(type);
    }

    @Override
    public Object parse( ArgumentSource source, Class type, ArgumentMatch... matches )  {
        if( matches.length > 1 )
            throw new UnsupportedOperationException("Only an input file name and validation stringency can be supplied when creating a BAM file reader.");

        SAMFileReaderBuilder builder = new SAMFileReaderBuilder();

        ArgumentMatch readerMatch = matches[0];

        if( readerMatch == null )
            throw new StingException("SAM file compression was supplied, but not associated writer was supplied with it.");
        if( readerMatch.values().size() > 1 )
            throw new StingException("Only one filename can be supplied per created BAM file");

        builder.setSAMFile(new File(readerMatch.values().get(0).trim()));

        return builder;
    }
}

/**
 * Handle SAMFileWriters.
 */
class SAMFileWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    private static final String COMPRESSION_FULLNAME = "bam_compression";
    private static final String COMPRESSION_SHORTNAME = "compress";

    @Override
    public boolean supports( Class type ) {
        return SAMFileWriter.class.isAssignableFrom(type);
    }

    @Override
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source, Argument description ) {
        String fullName = description.fullName().trim().length() > 0 ? description.fullName().trim() : "outputBAM";
        String shortName = description.shortName().trim().length() > 0 ? description.shortName().trim() : "ob";

        ArgumentDefinition writerDefinition = new ArgumentDefinition( source,
                                                                      fullName,
                                                                      shortName,
                                                                      getDoc( source, description ),
                                                                      isRequired( source, description ),
                                                                      getExclusiveOf( source, description ),
                                                                      getValidationRegex( source, description ) );
        ArgumentDefinition compressionDefinition = new ArgumentDefinition( source,
                                                                           COMPRESSION_FULLNAME,
                                                                           COMPRESSION_SHORTNAME,
                                                                           "Compression level to use for writing BAM files",
                                                                           false,
                                                                           "",
                                                                           "" );

        return Arrays.asList( writerDefinition, compressionDefinition );
    }

    @Override
    public Object parse( ArgumentSource source, Class type, ArgumentMatch... matches )  {
        if( matches.length > 2 )
            throw new UnsupportedOperationException("Only an input file name and validation stringency can be supplied when creating a BAM file reader.");

        SAMFileWriterBuilder builder = new SAMFileWriterBuilder();

        ArgumentMatch writerMatch = null;
        ArgumentMatch compressionMatch = null;

        for( ArgumentMatch match: matches ) {
            if( match.definition.fullName.equals(COMPRESSION_FULLNAME) )
                compressionMatch = match;
            else
                writerMatch = match;
        }

        if( writerMatch == null )
            throw new StingException("SAM file compression was supplied, but not associated writer was supplied with it.");
        if( writerMatch.values().size() > 1 )
            throw new StingException("Only one filename can be supplied per created BAM file");

        builder.setSAMFile(new File(writerMatch.values().get(0).trim()));

        if( compressionMatch != null ) {
            if( compressionMatch.values().size() > 1 )
                throw new StingException("Only one value can be supplied for BAM compression");
            int compressionLevel = Integer.valueOf(compressionMatch.values().get(0));
            builder.setCompressionLevel(compressionLevel);
        }

        return builder;
    }
    
}
