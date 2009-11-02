package org.broadinstitute.sting.utils;

import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.AbstractConfiguration;
import org.apache.log4j.Logger;

import java.util.Set;
import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import java.util.jar.Attributes;
import java.util.jar.JarFile;
import java.net.MalformedURLException;
import java.net.URL;
import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;

/**
 * PackageUtils contains some useful methods for package introspection.
 */
public class PackageUtils {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(PackageUtils.class);

    /**
     * A reference into our introspection utility.
     */
    private static Reflections reflections = null;

    static {
        // Initialize general-purpose source tree reflector.
        reflections = new Reflections(new AbstractConfiguration() {
            {
                setUrls(PackageUtils.getUrlsForClasspath(System.getProperty("java.class.path")));
                setScanners(new SubTypesScanner());
            }
        });
    }

    /**
     * Private constructor.  No instantiating this class!
     */
    private PackageUtils() {
    }

    {
    }

    /**
     * Return the classes that implement the specified interface.
     *
     * @param iface the interface which returned classes should implement.
     * @return the list of classes that implement the interface.
     */
    public static <T> List<Class<? extends T>> getClassesImplementingInterface(Class<T> iface) {
        // Load all classes implementing the given interface, then filter out any class that isn't concrete.
        Set<Class<? extends T>> allTypes = reflections.getSubTypesOf(iface);
        List<Class<? extends T>> concreteTypes = new ArrayList<Class<? extends T>>();
        for (Class<? extends T> type : allTypes) {
            if (JVMUtils.isConcrete(type))
                concreteTypes.add(type);
        }

        return concreteTypes;
    }

    /**
     * get a list of URL's, given the current classpath.  This function is a fix for the
     * reflections package, which doesn't correctly expand the classpath from a jar (it doesn't
     * get the manifest and parse out the embedded classpth).
     * @param classpath the current classpath
     * @return
     */
    public static Collection<URL> getUrlsForClasspath(String classpath) {
        // if the classpath is null, we can't really do anything
        if (classpath == null || classpath.equals(""))
            throw new StingException("Classpath cannot be empty or null");
        List<URL> urls = Lists.newArrayList();

        // the current working directory
        String baseDir = System.getProperty("user.dir");

        // our collection of classpath's to expand
        List<JarPath> javaClassPath = new ArrayList<JarPath>();

        // the current classpath can be a list of path's seperated by a semicolon
        String[] classPaths = classpath.split(File.pathSeparator);
        for (String part : classPaths) {
            extractJarClasspath(baseDir, javaClassPath, part);
        }

        // check to make sure each extracted path exists, if so add it to our list
        if (javaClassPath.size() > 0)
            for (JarPath path : javaClassPath) {
                try {
                    if (path.isValid())
                        urls.add(path.toValidFile().toURL());
                } catch (MalformedURLException e) {
                    throw new StingException("could not create url from " + path, e);
                }
            }

        return urls;
    }

    /**
     * extract the classpath from a jar
     *
     * @param baseDir       the base
     * @param javaClassPath the list of jar paths
     * @param part          the current the subsection of the classpath we're processing
     */
    private static void extractJarClasspath(String baseDir, List<JarPath> javaClassPath, String part) {
        try {
            JarFile myJar = new JarFile(part);
            Attributes.Name classPath = new Attributes.Name("Class-Path");
            if (myJar.getManifest().getMainAttributes().containsKey(classPath)) {
                for (String jar : myJar.getManifest().getMainAttributes().getValue(classPath).split(" "))
                    javaClassPath.add(new JarPath(baseDir, new File(part).getParent(), jar));
            }
        } catch (IOException e) {
            logger.warn("could not retreive manifest from " + part);
        }
    }

    /**
     * a simple helper class, to determine the absolute path to a jar file
     */
    static class JarPath {
        private final String mPath;
        private final String mFilename;
        private final String mWorkingDir;

        public JarPath(String workingDir, String path, String filename) {
            this.mPath = path;
            this.mFilename = filename;
            mWorkingDir = workingDir;
        }

        public File toValidFile() {
            if (new File(mFilename).exists())
                return new File(mFilename);
            if (new File(mPath + File.separator + mFilename).exists())
                return new File(mPath + File.separator + mFilename);
            if (new File(mWorkingDir + File.separator + mFilename).exists())
                return new File(mWorkingDir + File.separator + mFilename);
            return null;
        }

        public boolean isValid() {
            return (toValidFile() != null) ? true : false;
        }
    }

}
