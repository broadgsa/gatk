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

package org.broadinstitute.sting.utils.classloader;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.lang.reflect.Modifier;
import java.lang.reflect.Field;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
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
            throw new ReviewedStingException(String.format("Could not set %s in instance %s to %s",field.getName(),instance.getClass().getName(),value.toString()));
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
            throw new ReviewedStingException(String.format("Could not retrieve %s in instance %s",field.getName(),instance.getClass().getName()));
        }
    }

    /**
     * Gets a single object in the list matching or type-compatible with the given type.  Exceptions out if multiple objects match. 
     * @param type The desired type.
     * @param <T> The selected type.
     * @return A collection of the given arguments with the specified type.
     */
    public static <T> T getObjectOfType(Collection<Object> objectsToFilter, Class<T> type) {
        // TODO: Make JVM utils.
        Collection<T> selectedObjects = getObjectsOfType(objectsToFilter,type);
        if(selectedObjects.size() > 1)
            throw new ReviewedStingException("User asked for a single instance of the type, multiple were present");
        if(selectedObjects.size() == 0)
            throw new ReviewedStingException("User asked for a single instance of the type, but none were present");
        return selectedObjects.iterator().next();
    }

    /**
     * Gets a collection of all objects in the list matching or type-compatible with the given type. 
     * @param type The desired type.
     * @param <T> Again, the desired type.  Used so that clients can ignore type safety.
     * @return A collection of the given arguments with the specified type.
     */
    public static <T> Collection<T> getObjectsOfType(Collection<Object> objectsToFilter, Class<T> type) {
        Collection<T> selectedObjects = new ArrayList<T>();
        for(Object object: objectsToFilter) {
            if(type.isAssignableFrom(object.getClass()))
                selectedObjects.add((T)object);
        }
        return selectedObjects;
    }

}
