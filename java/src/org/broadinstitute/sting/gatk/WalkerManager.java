package org.broadinstitute.sting.gatk;

import net.sf.functionalj.reflect.StdReflect;
import net.sf.functionalj.reflect.JdkStdReflect;
import net.sf.functionalj.FunctionN;
import net.sf.functionalj.Functions;

import java.lang.reflect.Field;
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
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.apache.log4j.Logger;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 3:14:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class WalkerManager {

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(WalkerManager.class);

    private Map<String, Class> walkers;

    public WalkerManager(String pluginDirectory) {
        try {
            List<Class> walkerCandidates = new ArrayList<Class>();

            // Load all classes that live in this jar.
            final File location = getThisLocation();
            walkerCandidates.addAll(loadClassesFromLocation(location));

            // Load all classes that live in the extension path.
            if (pluginDirectory == null)
                pluginDirectory = location.getParent() + File.separator + "walkers";
            logger.info("plugin directory: " + pluginDirectory);

            File extensionPath = new File(pluginDirectory);
            if (extensionPath.exists()) {
                List<String> filesInPath = findFilesInPath(extensionPath, "", "class", false);
                walkerCandidates.addAll(loadExternalClasses(extensionPath, filesInPath));
            }

            walkerCandidates = filterWalkers(walkerCandidates);

            if (walkerCandidates.isEmpty())
                throw new RuntimeException("No walkers were found.");

            walkers = createWalkerDatabase(walkerCandidates);
        }
        // IOExceptions here are suspect; they indicate that the WalkerManager can't open its containing jar.
        // Wrap in a RuntimeException.
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    /**
     * Does a walker with the given name exist?
     *
     * @param walkerName Name of the walker for which to search.
     * @return True if the walker exists, false otherwise.
     */
    public boolean doesWalkerExist(String walkerName) {
        return walkers.containsKey(walkerName);
    }

    /**
     * Gets a walker with the given name, or null if no walker exists.
     *
     * @param walkerName Name of the walker to retrieve.
     * @return The walker object if found; null otherwise.
     */
    public Walker createWalkerByName(String walkerName)
            throws InstantiationException, IllegalAccessException {
        Class walker = walkers.get(walkerName);
        return (Walker) walker.newInstance();
    }

    public Class getWalkerClassByName(String walkerName) {
        return walkers.get(walkerName);
    }

    /**
     * Determines which jar file contains the WalkerManager class.
     *
     * @return Jar file containing the WalkerManager class.
     */
    private File getThisLocation() throws IOException {
        try {
            java.net.URI locationURI = getClass().getProtectionDomain().getCodeSource().getLocation().toURI();
            return new File(locationURI);
        }
        catch (java.net.URISyntaxException ex) {
            // a URISyntaxException here must be an IO error; wrap as such.
            throw new IOException(ex);
        }
    }

    /**
     * Load classes internal to the classpath from an arbitrary location.
     *
     * @param location Location from which to load classes.
     * @return List of classes.
     * @throws IOException Problem occurred reading classes.
     */
    private List<Class> loadClassesFromLocation(File location)
            throws IOException {
        if (location.getAbsolutePath().endsWith(".jar"))
            return loadClassesFromJar(location);
        else {
            List<String> classFileNames = findFilesInPath(location, "", "class", true);
            return loadInternalClasses(classFileNames);
        }
    }

    /**
     * Loads concrete classes from a jar which are both in the same package or 'sub-package' of baseClass,
     * and which extend from baseClass.
     *
     * @param jarFile The jar file to search.
     * @return A list of classes derived from baseClass.
     */
    private List<Class> loadClassesFromJar(final File jarFile)
            throws IOException {
        List<Class> subclasses = new ArrayList<Class>();

        JarInputStream jarInputStream = new JarInputStream(new FileInputStream(jarFile));

        try {
            JarEntry jarEntry = jarInputStream.getNextJarEntry();

            while (jarEntry != null) {
                String jarEntryName = jarEntry.getName();
                if (jarEntryName.endsWith(".class")) {
                    String className = fileNameToClassName(jarEntryName);
                    subclasses.add(Class.forName(className));
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

        return subclasses;
    }

    /**
     * Loads a list of classes currently on the classpath.
     *
     * @param classFileNames List of files representing classes.
     * @return class objects.
     * @throws IOException Unable to open any of the found classes.
     */
    private List<Class> loadInternalClasses(List<String> classFileNames)
            throws IOException {
        List<Class> internalClasses = new ArrayList<Class>();

        for (String classFileName : classFileNames) {
            String className = fileNameToClassName(classFileName);
            try {
                internalClasses.add(Class.forName(className));
            }
            catch (ClassNotFoundException ex) {
                // A ClassNotFoundException here must be an IO error; wrap as such.
                throw new IOException(ex);
            }
        }

        return internalClasses;
    }

    /**
     * Load loose classes, external to the classloader, from the specified directory.
     *
     * @param path source path from which to load classes.
     * @return A list of all loose classes contained in the path directory.
     */
    private List<Class> loadExternalClasses(final File path, List<String> classFileNames)
            throws IOException {
        List<Class> subclasses = new ArrayList<Class>();

        URL pathURL = path.toURI().toURL();

        ClassLoader cl = new URLClassLoader(new URL[]{pathURL});

        List<String> filesInPath = findFilesInPath(path, "", "class", false);
        for (String file : filesInPath) {
            String className = fileNameToClassName(file);
            try {
                subclasses.add(cl.loadClass(className));
            }
            catch (ClassNotFoundException ex) {
                // Class not found from a list of classes just looked up is an IO error.  Wrap and throw.
                throw new IOException(ex);
            }
        }

        return subclasses;
    }

    /**
     * Find the files in the given directory matching the given extension.
     *
     * @param basePath       Path to search.
     * @param relativePrefix What directory should the given files be presented relative to?
     * @param extension      Extension for which to search.
     * @param recursive      Search recursively.  Beware of symlinks!
     * @return A list of files matching the specified criteria.
     *         TODO: Move to a utils class.
     *         TODO: Test recursive traversal in the presence of a symlink.
     */
    private List<String> findFilesInPath(final File basePath, final String relativePrefix, final String extension, boolean recursive) {
        List<String> filesInPath = new ArrayList<String>();

        File[] contents = basePath.listFiles(new OrFilenameFilter(new DirectoryFilter(), new ExtensionFilter(extension)));
        for (File content : contents) {
            String relativeFileName = relativePrefix.trim().length() != 0 ?
                    relativePrefix + File.separator + content.getName() :
                    content.getName();
            if (relativeFileName.endsWith(extension))
                filesInPath.add(relativeFileName);
            else if (content.isDirectory() && recursive)
                filesInPath.addAll(findFilesInPath(content, relativeFileName, extension, recursive));
        }

        return filesInPath;
    }

    /**
     * Convert a filename of the form a/b/c.class to a.b.c.  Makes no assurances about whether the
     * class is valid on any classloader.
     *
     * @param fileName Filename to convert.
     * @return classname represented by that file.
     *         TODO: Move to a utils class.
     */
    private String fileNameToClassName(String fileName) {
        return fileName.substring(0, fileName.lastIndexOf(".class")).replace('/', '.');
    }

    /**
     * The following are general-purpose file selection filters.
     * TODO: Move to a utils class.
     */
    private class ExtensionFilter implements FilenameFilter {
        private String extensionName = null;

        public ExtensionFilter(String extensionName) {
            this.extensionName = extensionName;
        }

        public boolean accept(File f, String s) {
            return s.endsWith("." + extensionName);
        }
    }

    private class DirectoryFilter implements FilenameFilter {
        public boolean accept(File f, String s) {
            return new File(f, s).isDirectory();
        }
    }

    private class OrFilenameFilter implements FilenameFilter {
        private FilenameFilter lhs = null, rhs = null;

        public OrFilenameFilter(FilenameFilter lhs, FilenameFilter rhs) {
            this.lhs = lhs;
            this.rhs = rhs;
        }

        public boolean accept(File f, String s) {
            return lhs.accept(f, s) || rhs.accept(f, s);
        }
    }

    /**
     * Given a list of classes, return a list of those classes which extend from the Walker base interface.
     *
     * @param classes Arbitrary list of classes.
     * @return List of classes extending from Walker.
     */
    private List<Class> filterWalkers(List<Class> classes) {
        StdReflect reflect = new JdkStdReflect();
        FunctionN<Boolean> filterFunc = reflect.instanceFunction(new ClassFilter(Walker.class), "filter", Class.class);
        return Functions.findAll(filterFunc.f1(), classes);
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
            return baseClass.isAssignableFrom(clazz) &&
                    !Modifier.isAbstract(clazz.getModifiers()) &&
                    !Modifier.isInterface(clazz.getModifiers());
        }
    }

    /**
     * Instantiate the list of walker classes.  Add them to the walker hashmap.
     *
     * @param walkerClasses Classes to instantiate.
     * @return map of walker name to walker.
     */
    private Map<String, Class> createWalkerDatabase(List<Class> walkerClasses) {
        Map<String, Class> walkers = new HashMap<String, Class>();

        for (Class<Walker> walkerClass : walkerClasses) {
            String walkerName = getWalkerName(walkerClass);
            logger.info(String.format("* Adding module %s%n", walkerName));
            walkers.put(walkerName, walkerClass);
        }

        return walkers;
    }

    /**
     * Create a name for this type of walker.
     *
     * @param walkerType The type of walker.
     * @return A name for this type of walker.
     */
    public static String getWalkerName(Class<Walker> walkerType) {
        String walkerName = "";

        if (walkerType.getAnnotation(WalkerName.class) != null)
            walkerName = walkerType.getAnnotation(WalkerName.class).value().trim();

        if (walkerName.length() == 0) {
            walkerName = walkerType.getSimpleName();
            if (walkerName.endsWith("Walker"))
                walkerName = walkerName.substring(0, walkerName.lastIndexOf("Walker"));
        }

        return walkerName;
    }
}
