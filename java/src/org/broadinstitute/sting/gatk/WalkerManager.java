package org.broadinstitute.sting.gatk;

import net.sf.functionalj.reflect.StdReflect;
import net.sf.functionalj.reflect.JdkStdReflect;
import net.sf.functionalj.FunctionN;
import net.sf.functionalj.Functions;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;

import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.Allows;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.utils.JVMUtils;
import org.broadinstitute.sting.utils.PathUtils;
import org.broadinstitute.sting.utils.StingException;
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
            final File location = JVMUtils.getLocationFor( getClass() );
            walkerCandidates.addAll(loadClassesFromLocation(location));

            // Load all classes that live in the extension path.
            if (pluginDirectory == null)
                pluginDirectory = location.getParent() + File.separator + "walkers";
            logger.info("plugin directory: " + pluginDirectory);

            File extensionPath = new File(pluginDirectory);
            if (extensionPath.exists()) {
                List<String> classFilesInPath = PathUtils.findFilesInPath(extensionPath, "", "class", false);
                walkerCandidates.addAll(JVMUtils.loadExternalClasses(extensionPath, classFilesInPath));
                List<String> jarsInPath = PathUtils.findFilesInPath(extensionPath, "", "jar", false);
                for( String jarFileName: jarsInPath ) {
                    File jarFile = new File( extensionPath, jarFileName );
                    walkerCandidates.addAll(JVMUtils.loadExternalClassesFromJar(jarFile) );
                }
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

    /**
     * Retrieves the walker class given a walker name.
     * @param walkerName Name of the walker.
     * @return Class representing the walker.
     */
    public Class getWalkerClassByName(String walkerName) {
        return walkers.get(walkerName);
    }

    /**
     * Gets the data source for the provided walker.
     * @param walker The walker.
     * @return Which type of data source to traverse over...reads or reference?
     */
    public static DataSource getWalkerDataSource(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        By byDataSource = walkerClass.getAnnotation(By.class);
        if( byDataSource == null )
            throw new StingException("Unable to find By annotation for walker class " + walkerClass.getName());
        return byDataSource.value();
    }

    /**
     * Determine whether the given walker supports the given data source.
     * @param walker Walker to query.
     * @param dataSource Source to check for .
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Walker walker, DataSource dataSource) {
        Allows allowsDataSource = getWalkerAllowed(walker);

        // Allows is less restrictive than requires.  If an allows
        // clause is not specified, any kind of data is allowed.
        if( allowsDataSource == null )
            return true;

        return Arrays.asList(allowsDataSource.value()).contains(dataSource);
    }

    /**
     * Determine whether the given walker supports the given reference ordered data.
     * @param walker Walker to query.
     * @param rod Source to check.
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Walker walker, ReferenceOrderedData<? extends ReferenceOrderedDatum> rod) {
        Allows allowsDataSource = getWalkerAllowed(walker);

        // Allows is less restrictive than requires.  If an allows
        // clause is not specified, any kind of data is allowed.
        if( allowsDataSource == null )
            return true;

        // The difference between unspecified RMD and the empty set of metadata can't be detected.
        // Treat an empty 'allows' as 'allow everything'.  Maybe we can have a special RMD flag to account for this
        // case in the future.
        if( allowsDataSource.referenceMetaData().length == 0 )
            return true;

        for( RMD allowed: allowsDataSource.referenceMetaData() ) {
            if( rod.matches(allowed.name(),allowed.type()) )
                return true;
        }
        return false;
    }

    /**
     * Determine whether the given walker requires the given data source.
     * @param walker Walker to query.
     * @param dataSource Source to check for.
     * @return True if the walker allows this data type.  False otherwise.
     */
    public static boolean isRequired(Walker walker, DataSource dataSource) {
        Requires requiresDataSource = getWalkerRequirements(walker);
        return Arrays.asList(requiresDataSource.value()).contains(dataSource);
    }

    /**
     * Get a list of RODs required by the walker.
     * @param walker Walker to query.
     * @return True if the walker allows this data type.  False otherwise.
     */
    public static List<RMD> getRequiredMetaData(Walker walker) {
        Requires requiresDataSource = getWalkerRequirements(walker);
        return Arrays.asList(requiresDataSource.referenceMetaData());
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
            return JVMUtils.loadInternalClassesFromJar(location);
        else {
            List<String> classFileNames = PathUtils.findFilesInPath(location, "", "class", true);
            return JVMUtils.loadInternalClasses(classFileNames);
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
        FunctionN<Boolean> filterFunc = reflect.instanceFunction(new JVMUtils.ClassFilter(Walker.class), "filter", Class.class);
        return Functions.findAll(filterFunc.f1(), classes);
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
            logger.info(String.format("* Adding module %s", walkerName));
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
    public static String getWalkerName(Class<? extends Walker> walkerType) {
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

    /**
     * Utility to get the requires attribute from the walker.
     * Throws an exception if requirements are missing.
     * @param walker Walker to query for required data.
     * @return Required data attribute.
     */
    private static Requires getWalkerRequirements(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        Requires requiresDataSource = walkerClass.getAnnotation(Requires.class);
        if( requiresDataSource == null )
            throw new StingException( "Unable to find data types required by walker class " + walkerClass.getName());
        return requiresDataSource;
    }

    /**
     * Utility to get the forbidden attribute from the walker.
     * @param walker Walker to query for required data.
     * @return Required data attribute.  Null if forbidden info isn't present.
     */
    private static Allows getWalkerAllowed(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        Allows allowsDataSource = walkerClass.getAnnotation(Allows.class);
        return allowsDataSource;
    }
}
