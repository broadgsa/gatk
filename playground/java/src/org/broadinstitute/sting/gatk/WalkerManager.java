package org.broadinstitute.sting.gatk;

import java.lang.reflect.Modifier;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.jar.JarEntry;
import java.util.jar.JarInputStream;

import org.broadinstitute.sting.gatk.walkers.Walker;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 3:14:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class WalkerManager {
    private Map<String,Walker> walkers = null;

    public WalkerManager() {
        try {
            final File jarFile = getThisJarFile();
            List<Class> walkerClasses = loadClassesOfType( jarFile, Walker.class );

            if(walkerClasses.isEmpty())
                throw new RuntimeException("No walkers were found.");            

            walkers = instantiateWalkers( walkerClasses );
        }
        // IOExceptions here are suspect; they indicate that the WalkerManager can't open its containing jar.
        // Wrap in a RuntimeException.
        catch(IOException ex) {
            throw new RuntimeException(ex);
        }
        // The following two catches are more 'expected'; someone might add a walker that can't be instantiated.
        // TODO: Should these exceptions be handled differently?  Handling them like IOExceptions for the moment.
        catch(InstantiationException ex) {
            throw new RuntimeException(ex);
        }
        catch(IllegalAccessException ex) {
            throw new RuntimeException(ex);
        }
    }

    /**
     * Does a walker with the given name exist?
     * @param walkerName Name of the walker for which to search.
     * @return True if the walker exists, false otherwise.
     */
    public boolean doesWalkerExist(String walkerName) {
        return walkers.containsKey(walkerName);
    }

    /**
     * Gets a walker with the given name, or null if no walker exists.
     * @param walkerName Name of the walker to retrieve.
     * @return The walker object if found; null otherwise.
     */
    public Walker getWalkerByName(String walkerName) {
        return walkers.get(walkerName);
    }

    /**
     * Determines which jar file contains the WalkerManager class.
     * @return Jar file containing the WalkerManager class.
     */
    private File getThisJarFile() throws IOException {
        try {
            java.net.URI jarURI = getClass().getProtectionDomain().getCodeSource().getLocation().toURI();
            return new File( jarURI );
        }
        catch(java.net.URISyntaxException ex) {
            // a URISyntaxException here must be an IO error; wrap as such.
            throw new IOException(ex);
        }
    }

    /**
     * Loads concrete classes from a jar which are both in the same package or 'sub-package' of baseClass,
     * and which extend from baseClass.
     * @param jarFile The jar file to search.
     * @param baseClass The base class / interface for
     * @return A list of classes derived from baseClass.
     */
    private List<Class> loadClassesOfType(final File jarFile, final Class baseClass)
            throws IOException {
        final String PACKAGE_JAR_PREFIX = baseClass.getPackage().getName().replace('.','/');

        List<Class> subclasses = new ArrayList<Class>();

        JarInputStream jarInputStream = new JarInputStream(new FileInputStream( jarFile ) );

        try {
            JarEntry jarEntry = jarInputStream.getNextJarEntry();

            while(jarEntry != null) {
                String jarEntryName = jarEntry.getName();
                if(jarEntryName.startsWith(PACKAGE_JAR_PREFIX) && jarEntryName.endsWith(".class"))
                {
                    String className = jarEntryName.substring(0,jarEntryName.lastIndexOf(".class")).replace('/','.');
                    Class clazz = Class.forName(className);
                    if( baseClass.isAssignableFrom( clazz ) &&
                            !Modifier.isAbstract(clazz.getModifiers()) &&
                            !Modifier.isInterface(clazz.getModifiers()))
                        subclasses.add( clazz );
                }
                jarEntry = jarInputStream.getNextJarEntry();
            }
        }
        catch(ClassNotFoundException ex) {
            // A ClassNotFoundException here must be an IO error; wrap as such.
            throw new IOException(ex);
        }
        finally {
            jarInputStream.close();
        }

        return subclasses;
    }

    /**
     * Instantiate the list of walker classes.  Add them to the walker hashmap.
     * @param walkerClasses Classes to instantiate.
     * @throws InstantiationException
     * @throws IllegalAccessException
     */
    private Map<String,Walker> instantiateWalkers(List<Class> walkerClasses)
            throws InstantiationException, IllegalAccessException {
        Map<String,Walker> walkers = new HashMap<String,Walker>();

        for(Class walkerClass : walkerClasses) {
            Walker walker = (Walker)walkerClass.newInstance();
            String walkerName = walker.getName();

            System.out.printf("* Adding module %s%n", walkerName);            
            walkers.put(walkerName,walker);
        }

        return walkers;
    }
}
