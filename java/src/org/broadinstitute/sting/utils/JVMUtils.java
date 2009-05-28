package org.broadinstitute.sting.utils;

import java.io.File;
import java.io.IOException;
import java.io.FileInputStream;
import java.util.List;
import java.util.ArrayList;
import java.util.jar.JarInputStream;
import java.util.jar.JarEntry;
import java.net.URL;
import java.net.URLClassLoader;
import java.lang.reflect.Modifier;

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
    }

    /**
     * Loads concrete classes from a jar which are both in the same package or 'sub-package' of baseClass,
     * and which extend from baseClass.  Loaded classes must already be on the classpath.
     *
     * @param jarFile The jar file to search.
     * @return A list of classes derived from baseClass.
     */
    public static List<Class> loadInternalClassesFromJar(final File jarFile)
            throws IOException {
        return loadClassesFromJar( jarFile, new InternalLoadingStrategy() );
    }

    /**
     * Loads concrete classes from a jar which are both in the same package or 'sub-package' of baseClass,
     * and which extend from baseClass.  Loaded classes can be outside of the current classpath.
     *
     * @param jarFile The jar file to search.
     * @return A list of classes derived from baseClass.
     */
    public static List<Class> loadExternalClassesFromJar(final File jarFile)
        throws IOException {
        return loadClassesFromJar( jarFile, new ExternalLoadingStrategy(jarFile) );
    }

    /**
     * Loads a list of classes currently on the classpath.
     *
     * @param classFileNames List of files representing classes.
     * @return class objects.
     * @throws IOException Unable to open any of the found classes.
     */
    public static List<Class> loadInternalClasses(List<String> classFileNames)
            throws IOException {
        return loadClasses( classFileNames, new InternalLoadingStrategy() );
    }

    /**
     * Load loose classes, external to the classloader, from the specified directory.
     *
     * @param path source path from which to load classes.
     * @return A list of all loose classes contained in the path directory.
     */
    public static List<Class> loadExternalClasses(final File path, List<String> classFileNames)
            throws IOException {
        return loadClasses( classFileNames, new ExternalLoadingStrategy( path ) );
    }

    /**
     * Convert a filename of the form a/b/c.class to a.b.c.  Makes no assurances about whether the
     * class is valid on any classloader.
     *
     * @param fileName Filename to convert.
     * @return classname represented by that file.
     */
    public static String fileNameToClassName(String fileName) {
        return fileName.substring(0, fileName.lastIndexOf(".class")).replace('/', '.');
    }

    /**
     * Is the specified class a concrete implementation of baseClass?
     * @param clazz Class to check.
     * @param baseClass Base class to check against.
     * @return True if clazz implements baseClass and is not abstract / an interface.  False otherwise.
     */
    public static boolean isConcreteImplementationOf( Class clazz, Class baseClass ) {
        return baseClass.isAssignableFrom(clazz) &&
                !Modifier.isAbstract(clazz.getModifiers()) &&
                !Modifier.isInterface(clazz.getModifiers());
    }

    /**
     * Loads a list of classes from the given jar, using the provided loading strategy.
     * @param jarFile Jar file from which to load.
     * @param loader Dictates how these classes should be loaded.
     * @return A list of loaded classes.
     * @throws IOException In case there's an IO error trying to load the jar.
     */
    private static List<Class> loadClassesFromJar( final File jarFile, final LoadingStrategy loader )
            throws IOException {
        List<Class> classes = new ArrayList<Class>();

        JarInputStream jarInputStream = new JarInputStream(new FileInputStream(jarFile));

        try {
            JarEntry jarEntry = jarInputStream.getNextJarEntry();

            while (jarEntry != null) {
                String jarEntryName = jarEntry.getName();
                if (jarEntryName.endsWith(".class")) {
                    String className = fileNameToClassName(jarEntryName);
                    classes.add( loader.load( className ) );
                }
                jarEntry = jarInputStream.getNextJarEntry();
            }
        }
        catch (ClassNotFoundException ex) {
            // A ClassNotFoundException here must be an IO error; wrap as such.
            throw new IOException(ex);
        }
        finally {
            jarInputStream.close();
        }

        return classes;
    }

    /**
     * Loads a list of classes, using the provided loading strategy.
     * @param classFileNames Which class files to load.
     * @param loader Dictates how these classes should be loaded.
     * @return A list of loaded classes.
     * @throws IOException In case there's an IO error trying to load the jar.
     */
    private static List<Class> loadClasses( List<String> classFileNames, LoadingStrategy loader )
            throws IOException {
        List<Class> classes = new ArrayList<Class>();

        for (String classFileName : classFileNames) {
            String className = fileNameToClassName(classFileName);
            try {
                classes.add( loader.load( className ) );
            }
            catch (ClassNotFoundException ex) {
                // A ClassNotFoundException here must be an IO error; wrap as such.
                throw new IOException(ex);
            }
        }

        return classes;
    }


    /**
     * What mechanism should we use for loading a list of classes?
     */
    private static interface LoadingStrategy {
        Class load( String className ) throws ClassNotFoundException;
    }

    /**
     * An internal loading strategy, for loading classes already on the classpath.
     */
    private static class InternalLoadingStrategy implements LoadingStrategy {
        public Class load( String className )
                throws ClassNotFoundException {
            return Class.forName( className );
        }
    }

    /**
     * An external loading strategy, for loading classes not necessarily already
     * on the classpath.
     */
    private static class ExternalLoadingStrategy implements LoadingStrategy {
        private final ClassLoader classLoader;

        public ExternalLoadingStrategy( final File jarFile ) throws IOException {
            URL pathURL = jarFile.toURI().toURL();
            classLoader = new URLClassLoader(new URL[]{pathURL});
        }

        public Class load( String className ) throws ClassNotFoundException {
            return classLoader.loadClass(className);
        }
    }

}
