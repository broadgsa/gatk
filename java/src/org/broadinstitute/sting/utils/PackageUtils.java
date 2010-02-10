package org.broadinstitute.sting.utils;

import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.AbstractConfiguration;
import org.reflections.util.ClasspathHelper;

import java.util.Collections;
import java.util.Set;
import java.util.ArrayList;
import java.util.List;

/**
 * PackageUtils contains some useful methods for package introspection.
 */
public class PackageUtils {
    /**
     * A reference into our introspection utility.
     */
    private static Reflections reflections = null;

    static {
        // Initialize general-purpose source tree reflector.
        reflections = new Reflections( new AbstractConfiguration() {
            {
                setUrls(ClasspathHelper.getUrlsForCurrentClasspath());
                setScanners(new SubTypesScanner());
            }
        });
    }

    /**
     * Private constructor.  No instantiating this class!
     */
    private PackageUtils() {}
    {
    }

    /**
     * Return the classes that implement the specified interface.
     *
     * @param iface  the interface which returned classes should implement.
     * @return       the list of classes that implement the interface.
     */
    public static <T> List<Class<? extends T>> getClassesImplementingInterface(Class<T> iface) {
        // Load all classes implementing the given interface, then filter out any class that isn't concrete.
        Set<Class<? extends T>> allTypes = reflections.getSubTypesOf(iface);
        List<Class<? extends T>> concreteTypes = new ArrayList<Class<? extends T>>();
        for( Class<? extends T> type: allTypes ) {
            if( JVMUtils.isConcrete(type) )
                concreteTypes.add(type);
        }

        return concreteTypes;
    }
}
