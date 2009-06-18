package org.broadinstitute.sting.utils;

import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.io.IOException;

/**
 * PackageUtils contains some useful methods for package introspection.
 */
public class PackageUtils {
    /**
     * Private constructor.  No instantiating this class!
     */
    private PackageUtils() {}

    /**
     * Return the classes that implement the specified interface.
     *
     * @param iface  the interface which returned classes should implement.
     * @return       the list of classes that implement the interface.  How is that not clear by now?!!!!111one!!
     */
    public static ArrayList<Class> getClassesImplementingInterface(Class iface) {
        try {
            final File location = JVMUtils.getLocationFor(iface);

            List<Class> potentialClasses = getClassesFromLocation(location);
            ArrayList<Class> implementingClasses = new ArrayList<Class>();

            for (Class potentialClass : potentialClasses) {
                if (JVMUtils.isConcreteImplementationOf(potentialClass, iface)) {
                    implementingClasses.add(potentialClass);
                }
            }

            return implementingClasses;
        } catch (IOException e) {
            throw new StingException(String.format("Unable to inspect package containing '%s'", iface.getName()));
        }
    }

    /**
     * Return a list of classes at the specified location
     * @param location      the location  where we should start searching (as returned by JVMUtils.getLocationFor(Class))
     * @return              a list of classes at this location
     * @throws IOException  thrown if the jar or directory cannot be inspected
     */
    public static List<Class> getClassesFromLocation(File location) throws IOException {
        if (location.getAbsolutePath().endsWith(".jar"))
            return JVMUtils.loadInternalClassesFromJar(location);
        else {
            List<String> classFileNames = PathUtils.findFilesInPath(location, "", "class", true);
            return JVMUtils.loadInternalClasses(classFileNames);
        }
    }
}
