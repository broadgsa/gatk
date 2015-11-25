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

package org.broadinstitute.gatk.tools.walkers.help;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.StringEscapeUtils;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.help.GATKDocUtils;
import org.broadinstitute.gatk.utils.help.GenericDocumentationHandler;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.lang.annotation.Annotation;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class WalkerDocumentationHandler extends GenericDocumentationHandler {
    private final static String CMDLINE_GATK_URL = HelpConstants.GATK_DOCS_URL + "org_broadinstitute_gatk_engine_CommandLineGATK.php";

    @Override
    protected CommandLineProgram createCommandLineProgram() {
        return new CommandLineGATK();
    }

    /**
     * Umbrella function that groups the collection of values for specific annotations applied to an
     * instance of class c. Lists of collected values are added directly to the "toProcess" object.
     * Requires being able to instantiate the class.
     *
     * @param classToProcess the object to instantiate and query for the annotation
     * @param root the root of the document handler, to which we'll store collected annotations
     */
    @Override
    protected void getClazzAnnotations(Class classToProcess, Map<String, Object> root) {
        //
        // attempt to instantiate the class
        final Object instance = makeInstanceIfPossible(classToProcess);
        if (instance != null) {
            final Class myClass = instance.getClass();
            // Get parallelism options
            final HashSet<HashMap<String, Object>> parallelOptions = getParallelism(myClass, new HashSet<HashMap<String, Object>>());
            root.put("parallel", parallelOptions);
            // Get annotation info (what type of annotation, standard etc.)
            final HashSet<String> annotInfo = getAnnotInfo(myClass, new HashSet<String>());
            root.put("annotinfo", StringUtils.join(annotInfo, ", "));
            // Get annotation field (whether it goes in INFO or FORMAT)
            root.put("annotfield", getAnnotField(myClass));
            // Get walker type if applicable
            root.put("walkertype", getWalkerType(myClass));
            // Get partition type if applicable
            root.put("partitiontype", getPartitionType(myClass));
            // Get read filter annotations (ReadFilters) if applicable
            final HashSet<HashMap<String, Object>> bucket= getReadFilters(myClass, new HashSet<HashMap<String, Object>>());
            root.put("readfilters", bucket);
            // Get default downsampling settings
            final HashMap<String, Object> dsSettings = getDownSamplingSettings(myClass, new HashMap<String, Object>());
            root.put("downsampling", dsSettings);
            // Get reference window size settings
            final HashMap<String, Object> refwindow = getRefWindow(myClass, new HashMap<String, Object>());
            root.put("refwindow", refwindow);
            // Get ActiveRegion size settings
            final HashMap<String, Object> activeRegion = getActiveRegion(myClass, new HashMap<String, Object>());
            root.put("activeregion", activeRegion);
            // Get annotation header line description if applicable
            final Object annotDescriptLines = getAnnotDescript(instance, myClass);
            root.put("annotdescript", annotDescriptLines);

            // anything else?
        } else {
            // put empty items to avoid blowups
            root.put("parallel", new HashSet<String>());
            root.put("annotinfo", "");
            root.put("annotfield", "");
            root.put("walkertype", "");
            root.put("partitiontype", "");
            root.put("readfilters", new HashSet<HashMap<String, Object>>());
            root.put("downsampling", new HashMap<String, Object>());
            root.put("refwindow", new HashMap<String, Object>());
            root.put("activeregion", new HashMap<String, Object>());
            root.put("annotdescript", new ArrayList<HashMap<String, Object>>());
        }
    }

    /**
     * Utility function that looks up annotation descriptions if applicable.
     *
     * @param myClass the class to query
     * @return a hash map of descriptions, otherwise an empty map
     */
    private Object getAnnotDescript(Object instance, Class myClass) {
        //
        // Check if the class has the method we want
        for (Method classMethod : myClass.getMethods()) {
            if (classMethod.toString().contains("getDescriptions") && classMethod.toString().contains("annotator")) {
                try {
                    String headerLine = (classMethod.invoke(instance)).toString();
                    Pattern p = Pattern.compile("(INFO=<.*?>|FORMAT=<.*?>)");
                    Matcher m = p.matcher(headerLine);
                    List<String> annotLines = new ArrayList<>();
                    while (m.find()) {
                        annotLines.add(StringEscapeUtils.escapeHtml(m.group()));
                        System.out.println("found "+m.group());
                    }
                    return annotLines;
                } catch (IllegalArgumentException e) {
                } catch (IllegalAccessException e) {
                } catch (InvocationTargetException e) {
                }
            }
        }
        return null;
    }

    /**
     * Utility function that checks which parallelism options are available for an instance of class c.
     *
     * @param myClass the class to query for the interfaces
     * @param parallelOptions an empty HashSet in which to collect the info
     * @return a hash set of parallelism options, otherwise an empty set
     */
    private HashSet<HashMap<String, Object>> getParallelism(Class myClass, HashSet<HashMap<String, Object>> parallelOptions) {
        //
        // Retrieve interfaces
        Class[] implementedInterfaces = myClass.getInterfaces();
        for (Class intfClass : implementedInterfaces) {
            final HashMap<String, Object> nugget = new HashMap<String, Object>();
            if (intfClass.getSimpleName().equals("TreeReducible")) {
                nugget.put("name", intfClass.getSimpleName());
                nugget.put("arg", HelpConstants.ARG_TREEREDUCIBLE);
                nugget.put("link", CMDLINE_GATK_URL + "#" + HelpConstants.ARG_TREEREDUCIBLE);
            } else if (intfClass.getSimpleName().equals("NanoSchedulable")) {
                nugget.put("name", intfClass.getSimpleName());
                nugget.put("arg", HelpConstants.ARG_NANOSCHEDULABLE);
                nugget.put("link", CMDLINE_GATK_URL + "#" + HelpConstants.ARG_NANOSCHEDULABLE);
            } else {
                continue;
            }
            parallelOptions.add(nugget);
        }
        // Look up superclasses recursively
        final Class mySuperClass = myClass.getSuperclass();
        if (mySuperClass.getSimpleName().equals("Object")) {
            return parallelOptions;
        }
        return getParallelism(mySuperClass, parallelOptions);
    }

    /**
     * Utility function that looks up whether the annotation goes in INFO or FORMAT field.
     *
     * @param myClass the class to query for the interfaces
     * @return a String specifying the annotation field
     */
    private final String getAnnotField(Class myClass) {
        //
        // Look up superclasses recursively until we find either
        // GenotypeAnnotation or InfoFieldAnnotation
        final Class mySuperClass = myClass.getSuperclass();
        if (mySuperClass == InfoFieldAnnotation.class) {
            return "INFO (variant-level)";
        } else if (mySuperClass == GenotypeAnnotation.class) {
            return "FORMAT (sample genotype-level)";
        } else if (mySuperClass.getSimpleName().equals("Object")) {
            return "";
        }
        return getAnnotField(mySuperClass);
    }

    /**
     * Utility function that determines the annotation type for an instance of class c.
     *
     * @param myClass the class to query for the interfaces
     * @param annotInfo an empty HashSet in which to collect the info
     * @return a hash set of the annotation types, otherwise an empty set
     */
    private HashSet<String> getAnnotInfo(Class myClass, HashSet<String> annotInfo) {
        //
        // Retrieve interfaces
        Class[] implementedInterfaces = myClass.getInterfaces();
        for (Class intfClass : implementedInterfaces) {
            if (intfClass.getName().contains("Annotation")) {
                annotInfo.add(intfClass.getSimpleName());
            }
        }
        // Look up superclasses recursively
        final Class mySuperClass = myClass.getSuperclass();
        if (mySuperClass.getSimpleName().equals("Object")) {
            return annotInfo;
        }
        return getAnnotInfo(mySuperClass, annotInfo);
    }

    /**
     * Utility function that determines the default downsampling settings for an instance of class c.
     *
     * @param myClass the class to query for the settings
     * @param dsSettings an empty HashMap in which to collect the info
     * @return a hash set of the downsampling settings, otherwise an empty set
     */
    private HashMap<String, Object> getDownSamplingSettings(Class myClass, HashMap<String, Object> dsSettings) {
        //
        // Check for RODWalker first
        if (!checkForRODWalker(myClass).equals("yes")) {
            //
            // Retrieve annotation
            if (myClass.isAnnotationPresent(Downsample.class)) {
                final Annotation thisAnnotation = myClass.getAnnotation(Downsample.class);
                if(thisAnnotation instanceof Downsample) {
                    final Downsample dsAnnotation = (Downsample) thisAnnotation;
                    dsSettings.put("by", dsAnnotation.by().toString());
                    dsSettings.put("to_cov", dsAnnotation.toCoverage());
                }
            }
        }
        return dsSettings;
    }

    /**
     * Utility function that determines the reference window size for an instance of class c.
     *
     * @param myClass the class to query for the settings
     * @param refWindow an empty HashMap in which to collect the info
     * @return a HashMap of the window start and stop, otherwise an empty HashMap
     */
    private HashMap<String, Object> getRefWindow(Class myClass, HashMap<String, Object> refWindow) {
        //
        // Retrieve annotation
        if (myClass.isAnnotationPresent(Reference.class)) {
            final Annotation thisAnnotation = myClass.getAnnotation(Reference.class);
            if(thisAnnotation instanceof Reference) {
                final Reference refAnnotation = (Reference) thisAnnotation;
                refWindow.put("start", refAnnotation.window().start());
                refWindow.put("stop", refAnnotation.window().stop());
            }
        }
        return refWindow;
    }

    /**
     * Utility function that determines the ActiveRegion settings for an instance of class c.
     *
     * @param myClass the class to query for the settings
     * @param activeRegion an empty HashMap in which to collect the info
     * @return a HashMap of the ActiveRegion parameters, otherwise an empty HashMap
     */
    private HashMap<String, Object> getActiveRegion(Class myClass, HashMap<String, Object> activeRegion) {
        //
        // Retrieve annotation
        if (myClass.isAnnotationPresent(ActiveRegionTraversalParameters.class)) {
            final Annotation thisAnnotation = myClass.getAnnotation(ActiveRegionTraversalParameters.class);
            if(thisAnnotation instanceof ActiveRegionTraversalParameters) {
                final ActiveRegionTraversalParameters arAnnotation = (ActiveRegionTraversalParameters) thisAnnotation;
                activeRegion.put("ext", arAnnotation.extension());
                activeRegion.put("max", arAnnotation.maxRegion());
                activeRegion.put("min", arAnnotation.minRegion());
            }
        }
        return activeRegion;
    }

    /**
     * Utility function that determines the partition type of an instance of class c.
     *
     * @param myClass the class to query for the annotation
     * @return the partition type if applicable, otherwise an empty string
     */
    private String getPartitionType(Class myClass) {
        //
        // Retrieve annotation
        if (myClass.isAnnotationPresent(PartitionBy.class)) {
            final Annotation thisAnnotation = myClass.getAnnotation(PartitionBy.class);
            if(thisAnnotation instanceof PartitionBy) {
                final PartitionBy partAnnotation = (PartitionBy) thisAnnotation;
                return partAnnotation.value().toString();
            }
        }
        return "";
    }

    /**
     * Utility function that determines the type of walker subclassed by an instance of class c.
     *
     * @param myClass the class to query for the annotation
     * @return the type of walker if applicable, otherwise an empty string
     */
    private String getWalkerType(Class myClass) {
        //
        // Look up superclasses recursively until we find either Walker or Object
        final Class mySuperClass = myClass.getSuperclass();
        if (mySuperClass.getSimpleName().equals("Walker")) {
            return myClass.getSimpleName();
        } else if (mySuperClass.getSimpleName().equals("Object")) {
            return "";
        }
        return getWalkerType(mySuperClass);
    }

    /**
     * Utility function that checks whether an instance of class c is a subclass of RODWalker.
     *
     * @param myClass the class to query for the annotation
     * @return "yes" or "no" (can't use a Boolean because of the recursion)
     */
    private String checkForRODWalker(Class myClass) {
        //
        // Look up superclasses recursively until we find either RODWalker or (Walker or Object)
        final Class mySuperClass = myClass.getSuperclass();
        if (mySuperClass.getSimpleName().equals("RodWalker")) {
            return "yes";
        } else if (mySuperClass.getSimpleName().equals("Object") || mySuperClass.getSimpleName().equals("Walker")) {
            return "";
        }
        return checkForRODWalker(mySuperClass);
    }

    /**
     * Utility function that finds the values of ReadFilters annotation applied to an instance of class c.
     *
     * @param myClass the class to query for the annotation
     * @param bucket a container in which we store the annotations collected
     * @return a hash set of values, otherwise an empty set
     */
    private HashSet<HashMap<String, Object>> getReadFilters(Class myClass, HashSet<HashMap<String, Object>> bucket) {
        //
        // Retrieve annotation
        if (myClass.isAnnotationPresent(ReadFilters.class)) {
            final Annotation thisAnnotation = myClass.getAnnotation(ReadFilters.class);
            if(thisAnnotation instanceof ReadFilters) {
                final ReadFilters rfAnnotation = (ReadFilters) thisAnnotation;
                for (Class<?> filter : rfAnnotation.value()) {
                    // make hashmap of simplename and url
                    final HashMap<String, Object> nugget = new HashMap<String, Object>();
                    nugget.put("name", filter.getSimpleName());
                    nugget.put("filename", GATKDocUtils.phpFilenameForClass(filter));
                    bucket.add(nugget);
                }
            }
        }
        // Look up superclasses recursively
        final Class mySuperClass = myClass.getSuperclass();
        if (mySuperClass.getSimpleName().equals("Object")) {
            return bucket;
        }
        return getReadFilters(mySuperClass, bucket);
    }
}
