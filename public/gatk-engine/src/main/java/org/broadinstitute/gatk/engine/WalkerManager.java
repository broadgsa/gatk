/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine;

import org.broadinstitute.gatk.engine.filters.DisableableReadFilter;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.gatk.engine.filters.FilterManager;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.ResourceBundleExtractorDoclet;
import org.broadinstitute.gatk.utils.text.TextFormattingUtils;

import java.lang.annotation.Annotation;
import java.util.*;

/**
 * Plugin manager that also provides various utilities for inspecting Walkers.
 */
public class WalkerManager extends PluginManager<Walker> {

    /**
     * A collection of help text for walkers and their enclosing packages.
     */
    private ResourceBundle helpText;

    public WalkerManager() {
        super(Walker.class,"walker","");
        helpText = TextFormattingUtils.GATK_RESOURCE_BUNDLE;
    }

    /**
     * Get the list of walkers currently available to the GATK, organized
     * by package.
     * @param visibleWalkersOnly If true, return only the walker names that aren't hidden.
     * @return Names of currently available walkers.
     */
    public Map<String,Collection<Class<? extends Walker>>> getWalkerNamesByPackage(boolean visibleWalkersOnly) {
        Map<String,Collection<Class<? extends Walker>>> walkersByPackage = new HashMap<String,Collection<Class<? extends Walker>>>();
        for(Class<? extends Walker> walker: getPlugins()) {
            if(visibleWalkersOnly && isHidden(walker))
                continue;

            // Extract the name for the package; if the walker is in the unnamed package, use the empty string
            String walkerPackage = walker.getPackage() != null ? walker.getPackage().getName() : "";
            if(!walkersByPackage.containsKey(walkerPackage))
                walkersByPackage.put(walkerPackage,new ArrayList<Class<? extends Walker>>());
            walkersByPackage.get(walkerPackage).add(walker);
        }
        return Collections.unmodifiableMap(walkersByPackage);
    }

    /**
     * Gets the display name for a given package.
     * @param packageName Fully qualified package name.
     * @return A suitable display name for the package.
     */
    public String getPackageDisplayName(String packageName) {
        // ...try to compute the override from the text of the package name, while accounting for
        // unpackaged walkers.
        String displayName = packageName.substring(packageName.lastIndexOf('.')+1);
        if (displayName.trim().equals("")) displayName = "<unpackaged>";
        return displayName;
    }

    /**
     * Gets the help text associated with a given package name.
     * @param packageName Package for which to search for help text.
     * @return Package help text, or "" if none exists.
     */
    public String getPackageSummaryText(String packageName) {
        String key = String.format("%s.%s",packageName, ResourceBundleExtractorDoclet.SUMMARY_TAGLET_NAME);
        if(!helpText.containsKey(key))
            return "";
        return helpText.getString(key);
    }

    /**
     * Gets the summary help text associated with a given walker type.
     * @param walkerType Type of walker for which to search for help text.
     * @return Walker summary description, or "" if none exists.
     */
    public String getWalkerSummaryText(Class<? extends Walker> walkerType) {
        String walkerSummary = String.format("%s.%s",walkerType.getName(), ResourceBundleExtractorDoclet.SUMMARY_TAGLET_NAME);
        if(!helpText.containsKey(walkerSummary))
            return "";
        return helpText.getString(walkerSummary);
    }

    /**
     * Gets the summary help text associated with a given walker type.
     * @param walker Walker for which to search for help text.
     * @return Walker summary description, or "" if none exists.
     */
    public String getWalkerSummaryText(Walker walker) {
        return getWalkerSummaryText(walker.getClass());
    }

    /**
     * Gets the descriptive help text associated with a given walker type.
     * @param walkerType Type of walker for which to search for help text.
     * @return Walker full description, or "" if none exists.
     */
    public String getWalkerDescriptionText(Class<? extends Walker> walkerType) {
        String walkerDescription = String.format("%s.%s",walkerType.getName(), ResourceBundleExtractorDoclet.DESCRIPTION_TAGLET_NAME);
        if(!helpText.containsKey(walkerDescription))
            return "";
        return helpText.getString(walkerDescription);
    }

    /**
     * Gets the descriptive help text associated with a given walker type.
     * @param walker Walker for which to search for help text.
     * @return Walker full description, or "" if none exists.
     */
    public String getWalkerDescriptionText(Walker walker) {
        return getWalkerDescriptionText(walker.getClass());
    }

    /**
     * Retrieves the walker class given a walker name.
     * @param walkerName Name of the walker.
     * @return Class representing the walker.
     */
    public Class<? extends Walker> getWalkerClassByName(String walkerName) {
        return getPluginsByName().get(walkerName);
    }

    /**
     * Rather than use the default exception, return a MalformedWalkerArgumentsException.
     * @param errorMessage error message from formatErrorMessage()
     * @return - A MalformedWalkerArgumentsException with errorMessage
     */
    @Override
    protected UserException createMalformedArgumentException(final String errorMessage) {
        return new UserException.MalformedWalkerArgumentsException(errorMessage);
    }

    /**
     * Gets the data source for the provided walker.
     * @param walkerClass The class of the walker.
     * @return Which type of data source to traverse over...reads or reference?
     */
    public static DataSource getWalkerDataSource(Class<? extends Walker> walkerClass) {
        By byDataSource = walkerClass.getAnnotation(By.class);
        if( byDataSource == null )
            throw new ReviewedGATKException("Unable to find By annotation for walker class " + walkerClass.getName());
        return byDataSource.value();
    }

    /**
     * Gets the data source for the provided walker.
     * @param walker The walker.
     * @return Which type of data source to traverse over...reads or reference?
     */
    public static DataSource getWalkerDataSource(Walker walker) {
        return getWalkerDataSource(walker.getClass());
    }

    /**
     * Get a list of RODs allowed by the walker.
     * @param walkerClass Class of the walker to query.
     * @return The list of allowed reference meta data.
     */
    public static List<RMD> getAllowsMetaData(Class<? extends Walker> walkerClass) {
        return Collections.<RMD>emptyList();
    }

    /**
     * Determine whether the given walker supports the given data source.
     * @param walkerClass Class of the walker to query.
     * @param dataSource Source to check for .
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Class<? extends Walker> walkerClass, DataSource dataSource) {
        Allows allowsDataSource = getWalkerAllowed(walkerClass);

        // Allows is less restrictive than requires.  If an allows
        // clause is not specified, any kind of data is allowed.
        if( allowsDataSource == null )
            return true;

        return Arrays.asList(allowsDataSource.value()).contains(dataSource);
    }

    /**
     * Determine whether the given walker supports the given data source.
     * @param walker Walker to query.
     * @param dataSource Source to check for .
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Walker walker, DataSource dataSource) {
        return isAllowed(walker.getClass(), dataSource);
    }

    /**
     * Determine whether the given walker supports the given reference ordered data.
     * @param walkerClass Class of the walker to query.
     * @param rod Source to check.
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Class<? extends Walker> walkerClass, ReferenceOrderedDataSource rod) {
        return true;
    }

    /**
     * Determine whether the given walker supports the given reference ordered data.
     * @param walker Walker to query.
     * @param rod Source to check.
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Walker walker, ReferenceOrderedDataSource rod) {
        return isAllowed(walker.getClass(), rod);
    }

    /**
     * Determine whether the given walker requires the given data source.
     * @param walkerClass Class of the walker to query.
     * @param dataSource Source to check for.
     * @return True if the walker allows this data type.  False otherwise.
     */
    public static boolean isRequired(Class<? extends Walker> walkerClass, DataSource dataSource) {
        Requires requiresDataSource = getWalkerRequirements(walkerClass);
        return Arrays.asList(requiresDataSource.value()).contains(dataSource);
    }

    /**
     * Determine whether the given walker requires the given data source.
     * @param walker Walker to query.
     * @param dataSource Source to check for.
     * @return True if the walker allows this data type.  False otherwise.
     */
    public static boolean isRequired(Walker walker, DataSource dataSource) {
        return isRequired(walker.getClass(), dataSource);
    }

    /**
     * Get a list of RODs required by the walker.
     * @param walkerClass Class of the walker to query.
     * @return The list of required reference meta data.
     */
    public static List<RMD> getRequiredMetaData(Class<? extends Walker> walkerClass) {
        return Collections.emptyList();
    }

    /**
     * Get a list of RODs required by the walker.
     * @param walker Walker to query.
     * @return The list of required reference meta data.
     */
    public static List<RMD> getRequiredMetaData(Walker walker) {
        return getRequiredMetaData(walker.getClass());
    }

    /**
     * Reports whether this walker type is hidden -- in other words, whether it'll appear in the help output.
     * @param walkerType Class to test for visibility.
     * @return True if the walker should be hidden.  False otherwise.
     */
    public static boolean isHidden(Class<? extends Walker> walkerType) {
        return walkerType.isAnnotationPresent(Hidden.class);    
    }

    /**
     * Extracts filters that the walker has requested be run on the dataset.
     * @param walkerClass Class of the walker to inspect for filtering requests.
     * @param filterManager Manages the creation of filters.
     * @return A non-empty list of filters to apply to the reads.
     */
    public static List<ReadFilter> getReadFilters(Class<? extends Walker> walkerClass, FilterManager filterManager) {
        List<ReadFilter> filters = new ArrayList<ReadFilter>();
        for(Class<? extends ReadFilter> filterType: getReadFilterTypes(walkerClass))
            filters.add(filterManager.createFilterByType(filterType));
        return filters;
    }

    /**
     * Extracts filters that the walker has requested be run on the dataset.
     * @param walker Walker to inspect for filtering requests.
     * @param filterManager Manages the creation of filters.
     * @return A non-empty list of filters to apply to the reads.
     */
    public static List<ReadFilter> getReadFilters(Walker walker, FilterManager filterManager) {
        return getReadFilters(walker.getClass(), filterManager);
    }

    /**
     * Gets the type of downsampling method requested by the walker.  If an alternative
     * downsampling method is specified on the command-line, the command-line version will
     * be used instead.
     * @param walker The walker to interrogate.
     * @return The downsampling method, as specified by the walker.  Null if none exists.
     */
    public static DownsamplingMethod getDownsamplingMethod( Walker walker ) {
        return getDownsamplingMethod(walker.getClass());
    }

    /**
     * Gets the type of downsampling method requested by the walker.  If an alternative
     * downsampling method is specified on the command-line, the command-line version will
     * be used instead.
     * @param walkerClass The class of the walker to interrogate.
     * @return The downsampling method, as specified by the walker.  Null if none exists.
     */
    public static DownsamplingMethod getDownsamplingMethod( Class<? extends Walker> walkerClass ) {
        DownsamplingMethod downsamplingMethod = null;

        if( walkerClass.isAnnotationPresent(Downsample.class) ) {
            Downsample downsampleParameters = walkerClass.getAnnotation(Downsample.class);
            DownsampleType type = downsampleParameters.by();
            Integer toCoverage = downsampleParameters.toCoverage() >= 0 ? downsampleParameters.toCoverage() : null;
            Double toFraction = downsampleParameters.toFraction() >= 0.0d ? downsampleParameters.toFraction() : null;
            downsamplingMethod = new DownsamplingMethod(type, toCoverage, toFraction);
        }

        return downsamplingMethod;
    }

    public static <T extends Annotation> T getWalkerAnnotation(final Walker walker, final Class<T> clazz) {
        return walker.getClass().getAnnotation(clazz);
    }

    public static ReadTransformer.ApplicationTime getBAQApplicationTime(Walker walker) {
        return walker.getClass().getAnnotation(BAQMode.class).ApplicationTime();
    }    

    /**
     * Create a name for this type of walker.
     *
     * @param walkerType The type of walker.
     * @return A name for this type of walker.
     */
    @Override
    public String getName(Class walkerType) {
        String walkerName = "";

        if (walkerType.getAnnotation(WalkerName.class) != null)
            walkerName = ((WalkerName)walkerType.getAnnotation(WalkerName.class)).value().trim();
        else
            walkerName = super.getName(walkerType);

        return walkerName;
    }

    /**
     * Utility to get the requires attribute from the walker.
     * Throws an exception if requirements are missing.
     * @param walkerClass Class of the walker to query for required data.
     * @return Required data attribute.
     */
    private static Requires getWalkerRequirements(Class<? extends Walker> walkerClass) {
        Requires requiresDataSource = walkerClass.getAnnotation(Requires.class);
        if( requiresDataSource == null )
            throw new ReviewedGATKException( "Unable to find data types required by walker class " + walkerClass.getName());
        return requiresDataSource;
    }

    /**
     * Utility to get the requires attribute from the walker.
     * Throws an exception if requirements are missing.
     * @param walker Walker to query for required data.
     * @return Required data attribute.
     */
    private static Requires getWalkerRequirements(Walker walker) {
        return getWalkerRequirements(walker.getClass());
    }

    /**
     * Utility to get the forbidden attribute from the walker.
     * @param walkerClass Class of the walker to query for required data.
     * @return Required data attribute.  Null if forbidden info isn't present.
     */
    private static Allows getWalkerAllowed(Class<? extends Walker> walkerClass) {
        Allows allowsDataSource = walkerClass.getAnnotation(Allows.class);
        return allowsDataSource;
    }

    /**
     * Utility to get the forbidden attribute from the walker.
     * @param walker Walker to query for required data.
     * @return Required data attribute.  Null if forbidden info isn't present.
     */
    private static Allows getWalkerAllowed(Walker walker) {
        return getWalkerAllowed(walker.getClass());
    }

    /**
     * Gets the list of filtering classes specified as walker annotations.
     * @param walkerClass Class of the walker to inspect.
     * @return An array of types extending from SamRecordFilter.  Will never be null.
     */
    public static Collection<Class<? extends ReadFilter>> getReadFilterTypes(Class<?> walkerClass) {
        List<Class<? extends ReadFilter>> filterTypes = new ArrayList<Class<? extends ReadFilter>>();
        while(walkerClass != null) {
            // Add the read filters in the ReadFilters annotation
            if(walkerClass.isAnnotationPresent(ReadFilters.class)) {
                for ( Class c : walkerClass.getAnnotation(ReadFilters.class).value() ) {
                    if( !filterTypes.contains(c) )
                        filterTypes.add(c);
                }
            }
            // Remove read filters in the DisabledReadFilters annotation
            if(walkerClass.isAnnotationPresent(DisabledReadFilters.class)) {
                for ( Class c : walkerClass.getAnnotation(DisabledReadFilters.class).value() ) {
                    if ( filterTypes.contains(c) )
                        filterTypes.remove(c);
                }
            }
            walkerClass = walkerClass.getSuperclass();
        }
        return filterTypes;
    }

    /**
     * Gets the list of filtering classes specified as walker annotations.
     * @param walker The walker to inspect.
     * @return An array of types extending from SamRecordFilter.  Will never be null.
     */
    public static Collection<Class<? extends ReadFilter>> getReadFilterTypes(Walker walker) {
        return getReadFilterTypes(walker.getClass());
    }
}
