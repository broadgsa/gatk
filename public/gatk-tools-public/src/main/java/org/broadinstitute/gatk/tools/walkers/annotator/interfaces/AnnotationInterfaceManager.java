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

package org.broadinstitute.gatk.tools.walkers.annotator.interfaces;

import org.broadinstitute.gatk.utils.DeprecatedToolChecks;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.*;

public class AnnotationInterfaceManager {
    private static final String NULL_ANNOTATION_NAME = "none";
    private static final String NULL_ANNOTATION_GROUP_NAME = "none";

    private static PluginManager<InfoFieldAnnotation> infoFieldAnnotationPluginManager = new PluginManager<InfoFieldAnnotation>(InfoFieldAnnotation.class);
    private static PluginManager<GenotypeAnnotation> genotypeAnnotationPluginManager = new PluginManager<GenotypeAnnotation>(GenotypeAnnotation.class);
    private static PluginManager<AnnotationType> annotationTypePluginManager = new PluginManager<AnnotationType>(AnnotationType.class);

    public static List<InfoFieldAnnotation> createAllInfoFieldAnnotations() {
        return infoFieldAnnotationPluginManager.createAllTypes();
    }

    public static List<GenotypeAnnotation> createAllGenotypeAnnotations() {
        return genotypeAnnotationPluginManager.createAllTypes();
    }

    public static void validateAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse) {
        HashMap<String, Class> classMap = new HashMap<String, Class>();
        for ( Class c : infoFieldAnnotationPluginManager.getPlugins() )
            classMap.put(c.getSimpleName(), c);
        for ( Class c : genotypeAnnotationPluginManager.getPlugins() )
            classMap.put(c.getSimpleName(), c);
        for ( Class c : annotationTypePluginManager.getInterfaces() )
            classMap.put(c.getSimpleName(), c);

        if ( annotationGroupsToUse.size() != 1 || !NULL_ANNOTATION_GROUP_NAME.equals(annotationGroupsToUse.get(0)) ) {
            for ( String group : annotationGroupsToUse ) {
                Class interfaceClass = classMap.get(group);
                if ( interfaceClass == null )
                    interfaceClass = classMap.get(group + "Annotation");
                if ( interfaceClass == null )
                    throw new UserException.BadArgumentValue("group", "Annotation group " + group + " was not found; please check that you have specified the group name correctly");
            }
        }

        // validate the specific classes provided
        if ( annotationsToUse.size() != 1 || !NULL_ANNOTATION_NAME.equals(annotationsToUse.get(0)) ) {
            for (String annotation : annotationsToUse) {
                Class annotationClass = classMap.get(annotation);
                if (annotationClass == null)
                    annotationClass = classMap.get(annotation + "Annotation");
                if (annotationClass == null) {
                    if (DeprecatedToolChecks.isDeprecatedAnnotation(annotation)) {
                        throw new UserException.DeprecatedAnnotation(annotation, DeprecatedToolChecks.getAnnotationDeprecationInfo(annotation));
                    } else {
                        throw new UserException.BadArgumentValue("annotation", "Annotation " + annotation + " was not found; please check that you have specified the annotation name correctly");
                    }
                }
            }
        }
    }

    public static List<InfoFieldAnnotation> createInfoFieldAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse) {
        return createAnnotations(infoFieldAnnotationPluginManager, annotationGroupsToUse, annotationsToUse);
    }

    public static List<GenotypeAnnotation> createGenotypeAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse) {
        return createAnnotations(genotypeAnnotationPluginManager, annotationGroupsToUse, annotationsToUse);
    }

    private static <T> List<T> createAnnotations(PluginManager<T> pluginManager, List<String> annotationGroupsToUse, List<String> annotationsToUse) {
        // get the instances
        List<T> annotations = new ArrayList<T>();

        // get the classes from the provided groups (interfaces)
        // create a map for all annotation classes which implement our top-level interfaces
        HashMap<String, Class> classMap = new HashMap<String, Class>();
        for ( Class c : pluginManager.getPlugins() )
            classMap.put(c.getSimpleName(), c);
        for ( Class c : annotationTypePluginManager.getInterfaces() )
            classMap.put(c.getSimpleName(), c);

        // use a TreeSet so that classes are returned deterministically (the plugin manager apparently isn't deterministic)
        TreeSet<Class> classes = new TreeSet<Class>(new Comparator<Class>() {
            public int compare(Class o1, Class o2) {
                return o1.getSimpleName().compareTo(o2.getSimpleName());
            }
        });

        if ( annotationGroupsToUse.size() != 1 || !NULL_ANNOTATION_GROUP_NAME.equals(annotationGroupsToUse.get(0)) ) {
            for ( String group : annotationGroupsToUse ) {
                Class interfaceClass = classMap.get(group);
                if ( interfaceClass == null )
                    interfaceClass = classMap.get(group + "Annotation");
                if ( interfaceClass != null )
                    classes.addAll(pluginManager.getPluginsImplementing(interfaceClass));
            }
        }

        // get the specific classes provided
        if ( annotationsToUse.size() != 1 || !NULL_ANNOTATION_NAME.equals(annotationsToUse.get(0)) ) {
            for (String annotation : annotationsToUse) {
                Class annotationClass = classMap.get(annotation);
                if (annotationClass == null)
                    annotationClass = classMap.get(annotation + "Annotation");
                if (annotationClass != null)
                    classes.add(annotationClass);
            }
        }

        // note that technically an annotation can work on both the INFO and FORMAT fields
        for ( Class c : classes )
            annotations.add((T)pluginManager.createByType(c));

        return annotations;
    }
}
