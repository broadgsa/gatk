package org.broadinstitute.sting.utils;

import java.lang.reflect.Modifier;
import java.lang.reflect.Field;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 30, 2009
 * Time: 5:38:05 PM
 *
 * A set of static utility methods for determining information about this runtime environment.
 * Introspects classes, loads jars, etc.
 */
public class JVMUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private JVMUtils() { }

    /**
     * Determines which location contains the specified class.
     *
     * @return Location (either jar file or directory) of path containing class.
     */
    public static File getLocationFor( Class clazz ) throws IOException {
        try {
            java.net.URI locationURI = clazz.getProtectionDomain().getCodeSource().getLocation().toURI();
            return new File(locationURI);
        }
        catch (java.net.URISyntaxException ex) {
            // a URISyntaxException here must be an IO error; wrap as such.
            throw new IOException(ex);
        }
        catch ( NullPointerException ne ) {
        	throw new IOException("Can not extract code source location for "+clazz.getName());
        }
    }    

    /**
     * Is the specified class a concrete implementation of baseClass?
     * @param clazz Class to check.
     * @return True if clazz is concrete.  False otherwise.
     */
    public static boolean isConcrete( Class clazz ) {
        return !Modifier.isAbstract(clazz.getModifiers()) &&
               !Modifier.isInterface(clazz.getModifiers());
    }

    /**
     * Retrieve all fields available in this object, regardless of where they are declared or
     * whether they're accessible.
     * @param type Type to inspect for fields.
     * @return A list of all available fields.
     */
    public static List<Field> getAllFields(Class type) {
        List<Field> allFields = new ArrayList<Field>();
        while( type != null ) {
            allFields.addAll(Arrays.asList(type.getDeclaredFields()));
            type = type.getSuperclass();
        }
        return allFields;
    }

    /**
     * Find the field with the given name in the class.  Will inspect all fields, independent
     * of access level.
     * @param type Class in which to search for the given field.
     * @param fieldName Name of the field for which to search.
     * @return The field, or null if no such field exists.
     */
    public static Field findField( Class type, String fieldName ) {
        while( type != null ) {
            Field[] fields = type.getDeclaredFields();
            for( Field field: fields ) {
                if( field.getName().equals(fieldName) )
                    return field;
            }
            type = type.getSuperclass();
        }
        return null;
    }

    /**
     * Sets the provided field in the given instance to the given value.  Circumvents access restrictions:
     * a field can be private and still set successfully by this function.
     * @param field Field to set in the given object.
     * @param instance Instance in which to set the field.
     * @param value The value to which to set the given field in the given instance.
     */
    public static void setFieldValue( Field field, Object instance, Object value ) {
        try {
            field.setAccessible(true);
            field.set(instance, value);
        }
        catch( IllegalAccessException ex ) {
            throw new StingException(String.format("Could not set %s in instance %s to %s",field.getName(),instance.getClass().getName(),value.toString()));
        }
    }

    /**
     * Gets the value stored in the provided field in the given instance.
     * @param field Field to set in the given object.
     * @param instance Instance in which to set the field.
     * @return Value currently stored in the given field.
     */
    public static Object getFieldValue( Field field, Object instance ) {
        try {
            field.setAccessible(true);
            return field.get(instance);
        }
        catch( IllegalAccessException ex ) {
            throw new StingException(String.format("Could not retrieve %s in instance %s",field.getName(),instance.getClass().getName()));
        }
    }

}
