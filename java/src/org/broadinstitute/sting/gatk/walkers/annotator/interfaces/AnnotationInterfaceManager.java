package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.*;

public class AnnotationInterfaceManager {
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

        if ( annotationGroupsToUse.size() != 1 || !"none".equals(annotationGroupsToUse.get(0)) ) {
            for ( String group : annotationGroupsToUse ) {
                Class interfaceClass = classMap.get(group);
                if ( interfaceClass == null )
                    interfaceClass = classMap.get(group + "Annotation");
                if ( interfaceClass == null )
                    throw new UserException.BadArgumentValue("group", "Class " + group + " is not found; please check that you have specified the class name correctly");
            }
        }

        // validate the specific classes provided
        for ( String annotation : annotationsToUse ) {
            Class annotationClass = classMap.get(annotation);
            if ( annotationClass == null )
                annotationClass = classMap.get(annotation + "Annotation");
            if ( annotationClass == null )
                throw new UserException.BadArgumentValue("annotation", "Class " + annotation + " is not found; please check that you have specified the class name correctly");
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

        if ( annotationGroupsToUse.size() != 1 || !"none".equals(annotationGroupsToUse.get(0)) ) {
            for ( String group : annotationGroupsToUse ) {
                Class interfaceClass = classMap.get(group);
                if ( interfaceClass == null )
                    interfaceClass = classMap.get(group + "Annotation");
                if ( interfaceClass != null )
                    classes.addAll(pluginManager.getPluginsImplementing(interfaceClass));
            }
        }

        // get the specific classes provided
        for ( String annotation : annotationsToUse ) {
            Class annotationClass = classMap.get(annotation);
            if ( annotationClass == null )
                annotationClass = classMap.get(annotation + "Annotation");
            if ( annotationClass != null )
                classes.add(annotationClass);
        }

        // note that technically an annotation can work on both the INFO and FORMAT fields
        for ( Class c : classes )
            annotations.add(pluginManager.createByType(c));

        return annotations;
    }
}
