package org.broadinstitute.sting.utils;

import java.lang.reflect.Modifier;
import java.io.File;
import java.io.IOException;

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

}
