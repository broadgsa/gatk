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
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.reflections.util.ClasspathHelper;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.*;
import java.net.URL;
import java.util.*;

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
     * @param clazz The specified class.
     * @return Location (either jar file or directory) of path containing class.
     * @throws IOException when the URI cannot be found.
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
     * Is the specified class anonymous?  The plugin manager (for one) generally requires that
     * plugin classes be named so that they can easily be specified from the command line.
     * @param clazz The class on which to perform the anonymous check.
     * @return True if the class is anonymous; false otherwise.
     */
    public static boolean isAnonymous(Class clazz) {
        return clazz.isAnonymousClass();
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
     * @param objectsToFilter objects to filter.
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
     * @param objectsToFilter objects to filter.
     * @param type The desired type.
     * @param <T> Again, the desired type.  Used so that clients can ignore type safety.
     * @return A collection of the given arguments with the specified type.
     */
    @SuppressWarnings("unchecked")
    public static <T> Collection<T> getObjectsOfType(Collection<Object> objectsToFilter, Class<T> type) {
        Collection<T> selectedObjects = new ArrayList<T>();
        for(Object object: objectsToFilter) {
            if(type.isAssignableFrom(object.getClass()))
                selectedObjects.add((T)object);
        }
        return selectedObjects;
    }

    /**
     * Returns the list of class path urls.
     * @return the list of class path urls.
     */
    public static Set<URL> getClasspathURLs() {
        return ClasspathHelper.getUrlsForManifestsCurrentClasspath();
    }

    /**
     * Adds all the generic types from a class definition to the collection.
     * Does not inspect the methods or fields, only the class.
     * @param classes Set to collect the classes.
     * @param type Type to inspect.
     */
    public static void addGenericTypes(Set<Class<?>> classes, Type type) {
        if (type instanceof ParameterizedType) {
            ParameterizedType parameterizedType = (ParameterizedType)type;
            for (Type actualType: parameterizedType.getActualTypeArguments())
                addGenericTypes(classes, actualType);
        } else if (type instanceof GenericArrayType) {
            addGenericTypes(classes, ((GenericArrayType)type).getGenericComponentType());
        } else if (type instanceof WildcardType) {
            WildcardType wildcardType = (WildcardType)type;
            for (Type upperType: wildcardType.getUpperBounds())
                addGenericTypes(classes, upperType);
            for (Type lowerType: wildcardType.getLowerBounds())
                addGenericTypes(classes, lowerType);
        } else if (type instanceof Class<?>) {
            classes.add((Class<?>) type);
        } else {
            throw new StingException("Unknown type: " + type + " (" + type.getClass().getName() + ")");
        }
    }

    public static Class getParameterizedTypeClass(Type t) {
        if ( t instanceof ParameterizedType ) {
            ParameterizedType parameterizedType = (ParameterizedType)t;
            if ( parameterizedType.getActualTypeArguments().length != 1 )
                throw new ReviewedStingException("BUG: more than 1 generic type found on class" + t);
            return (Class)parameterizedType.getActualTypeArguments()[0];
        } else
            throw new ReviewedStingException("BUG: could not find generic type on class " + t);
    }
}
