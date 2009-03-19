package org.broadinstitute.sting.gatk;

import net.sf.functionalj.reflect.StdReflect;
import net.sf.functionalj.reflect.JdkStdReflect;
import net.sf.functionalj.FunctionN;
import net.sf.functionalj.Functions;

import java.lang.reflect.Modifier;
import java.io.File;
import java.io.FilenameFilter;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLClassLoader;
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

    public WalkerManager( String pluginDirectory ) {
        try {
            final File jarFile = getThisJarFile();

            if(pluginDirectory == null)
                pluginDirectory = jarFile.getParent() + File.separator + "walkers";

            System.out.println("plugin directory: " + pluginDirectory);

            List<Class> walkerClasses = new ArrayList<Class>();

            // Load all classes that live in this jar.
            walkerClasses.addAll( loadClassesFromJar( jarFile ) );

            // Load all classes that live in the extension path.
            File extensionPath = new File( pluginDirectory );
            if(extensionPath.exists())
                walkerClasses.addAll( loadClassesFromPath( extensionPath ) );

            walkerClasses = filterWalkers(walkerClasses);

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
     * @return A list of classes derived from baseClass.
     */
    private List<Class> loadClassesFromJar(final File jarFile)
            throws IOException {
        List<Class> subclasses = new ArrayList<Class>();

        JarInputStream jarInputStream = new JarInputStream(new FileInputStream( jarFile ) );

        try {
            JarEntry jarEntry = jarInputStream.getNextJarEntry();

            while(jarEntry != null) {
                String jarEntryName = jarEntry.getName();
                if(jarEntryName.endsWith(".class"))
                {
                    String className = jarEntryName.substring(0,jarEntryName.lastIndexOf(".class")).replace('/','.');
                    subclasses.add( Class.forName(className) );
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
     * Load loose classes from the specified directory.
     * @param path source path from which to load classes.
     * @return A list of all loose classes contained in the path directory.
     */
    private List<Class> loadClassesFromPath(final File path)
            throws IOException {
        List<Class> subclasses = new ArrayList<Class>();

        URL pathURL = path.toURI().toURL();

        ClassLoader cl = new URLClassLoader(new URL[] { pathURL });

        File[] classFiles = path.listFiles(
                new FilenameFilter() { public boolean accept( File f, String s ) { return s.endsWith(".class"); } });
        for(File classFile: classFiles) {
            // Poor person's way of getting the classname: assume classes are unpackaged.
            // Chop out the classname section of the URL, and assume the class is unpackaged.
            String fileName = classFile.getName();
            String className = fileName.substring(fileName.lastIndexOf(File.separator)+1,fileName.lastIndexOf(".class"));

            try {
                subclasses.add(cl.loadClass(className));
            }
            catch(ClassNotFoundException ex) {
                // Class not found from a list of classes just looked up is an IO error.  Wrap and throw.
                throw new IOException(ex);
            }

        }

        return subclasses;
    }

    /**
     * Given a list of classes, return a list of those classes which extend from the Walker base interface.
     * @param classes Arbitrary list of classes.
     * @return List of classes extending from Walker.
     */
    private List<Class> filterWalkers(List<Class> classes) {
        StdReflect reflect = new JdkStdReflect();
        FunctionN<Boolean> filterFunc = reflect.instanceFunction(new ClassFilter(Walker.class),"filter",Class.class);
        return Functions.findAll(filterFunc.f1(),classes);
    }

    /**
     * A functor returning true for classes which extend from baseClass.
     */
    private class ClassFilter {
        private Class baseClass;

        public ClassFilter(Class baseClass) {
            this.baseClass = baseClass;
        }

        public Boolean filter(Class clazz) {
            return baseClass.isAssignableFrom( clazz ) &&
                    !Modifier.isAbstract( clazz.getModifiers() ) &&
                    !Modifier.isInterface( clazz.getModifiers() );
        }
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
