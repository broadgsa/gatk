package org.broadinstitute.sting.gatk;

import java.util.*;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.PackageUtils;
import org.apache.log4j.Logger;
import net.sf.picard.filter.SamRecordFilter;

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

    private Map<String, Class<? extends Walker>> walkersByName;

    public WalkerManager() {
        List<Class<? extends Walker>> walkers = PackageUtils.getClassesImplementingInterface(Walker.class);
        walkersByName = createWalkerDatabase(walkers);
    }

    /**
     * Does a walker with the given name exist?
     *
     * @param walkerName Name of the walker for which to search.
     * @return True if the walker exists, false otherwise.
     */
    public boolean doesWalkerExist(String walkerName) {
        return walkersByName.containsKey(walkerName);
    }

    /**
     * Gets a walker with the given name, or null if no walker exists.
     *
     * @param walkerName Name of the walker to retrieve.
     * @return The walker object if found; null otherwise.
     */
    public Walker createWalkerByName(String walkerName)
            throws InstantiationException, IllegalAccessException {
        Class<? extends Walker> walker = walkersByName.get(walkerName);
        if( walker == null )
            throw new StingException(String.format("Could not find walker with name: %s", walkerName));
        return walker.newInstance();
    }

    /**
     * Retrieves the walker class given a walker name.
     * @param walkerName Name of the walker.
     * @return Class representing the walker.
     */
    public Class getWalkerClassByName(String walkerName) {
        return walkersByName.get(walkerName);
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
     * Extracts filters that the walker has requested be run on the dataset.
     * @param walker Walker to inspect for filtering requests.
     * @return A non-empty list of filters to apply to the reads.
     */
    public static List<SamRecordFilter> getReadFilters(Walker walker) {
        Class<? extends SamRecordFilter>[] filterTypes = getReadFilterTypes(walker);
        List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();

        for( Class<? extends SamRecordFilter> filterType: filterTypes ) {
            try {
                filters.add(filterType.newInstance());
            }
            catch( InstantiationException ex ) {
                throw new StingException("Unable to instantiate filter: " + filterType, ex);
            }
            catch( IllegalAccessException ex ) {
                throw new StingException("Unable to access filter: " + filterType, ex);                
            }
        }

        return filters;
    }

    /**
     * Instantiate the list of walker classes.  Add them to the walker hashmap.
     *
     * @param walkerClasses Classes to instantiate.
     * @return map of walker name to walker.
     */
    private Map<String, Class<? extends Walker>> createWalkerDatabase(List<Class<? extends Walker>> walkerClasses) {
        Map<String, Class<? extends Walker>> walkers = new HashMap<String, Class<? extends Walker>>();

        for (Class<? extends Walker> walkerClass : walkerClasses) {
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

    /**
     * Gets the list of filtering classes specified as walker annotations.
     * @param walker The walker to inspect.
     * @return An array of types extending from SamRecordFilter.  Will never be null.
     */
    private static Class<? extends SamRecordFilter>[] getReadFilterTypes(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        if( !walkerClass.isAnnotationPresent(ReadFilters.class) )
            return new Class[0];
        return walkerClass.getAnnotation(ReadFilters.class).value();
    }
}
