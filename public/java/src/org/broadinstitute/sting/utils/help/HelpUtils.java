/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.help;

import com.sun.javadoc.FieldDoc;
import com.sun.javadoc.PackageDoc;
import com.sun.javadoc.ProgramElementDoc;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotationType;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;

import java.lang.reflect.Field;
import java.util.List;

public class HelpUtils {

    protected static boolean assignableToClass(ProgramElementDoc classDoc, Class lhsClass, boolean requireConcrete) {
        try {
            Class type = getClassForDoc(classDoc);
            return lhsClass.isAssignableFrom(type) && (!requireConcrete || JVMUtils.isConcrete(type));
        } catch (Throwable t) {
            // Ignore errors.
            return false;
        }
    }

    protected static Class getClassForDoc(ProgramElementDoc doc) throws ClassNotFoundException {
        return Class.forName(getClassName(doc));
    }

    protected static Field getFieldForFieldDoc(FieldDoc fieldDoc) {
        try {
            Class clazz = getClassForDoc(fieldDoc.containingClass());
            return JVMUtils.findField(clazz, fieldDoc.name());
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Reconstitute the class name from the given class JavaDoc object.
     *
     * @param doc the Javadoc model for the given class.
     * @return The (string) class name of the given class.
     */
    protected static String getClassName(ProgramElementDoc doc) {
        PackageDoc containingPackage = doc.containingPackage();
        return containingPackage.name().length() > 0 ?
                String.format("%s.%s", containingPackage.name(), doc.name()) :
                String.format("%s", doc.name());
    }

    /**
     * Simple method to print a list of available annotations.
     */
    public static void listAnnotations() {
        System.out.println("\nThis is a list of available Variant Annotations for use with tools such as UnifiedGenotyper, HaplotypeCaller and VariantAnnotator. Please see the Technical Documentation for more details about these annotations:");
        System.out.println("http://www.broadinstitute.org/gatk/gatkdocs/");
        System.out.println("\nStandard annotations in the list below are marked with a '*'.");
        List<Class<? extends InfoFieldAnnotation>> infoAnnotationClasses = new PluginManager<InfoFieldAnnotation>(InfoFieldAnnotation.class).getPlugins();
        System.out.println("\nAvailable annotations for the VCF INFO field:");
        for (int i = 0; i < infoAnnotationClasses.size(); i++)
            System.out.println("\t" + (StandardAnnotation.class.isAssignableFrom(infoAnnotationClasses.get(i)) ? "*" : "") + infoAnnotationClasses.get(i).getSimpleName());
        System.out.println();
        List<Class<? extends GenotypeAnnotation>> genotypeAnnotationClasses = new PluginManager<GenotypeAnnotation>(GenotypeAnnotation.class).getPlugins();
        System.out.println("\nAvailable annotations for the VCF FORMAT field:");
        for (int i = 0; i < genotypeAnnotationClasses.size(); i++)
            System.out.println("\t" + (StandardAnnotation.class.isAssignableFrom(genotypeAnnotationClasses.get(i)) ? "*" : "") + genotypeAnnotationClasses.get(i).getSimpleName());
        System.out.println();
        System.out.println("\nAvailable classes/groups of annotations:");
        for ( Class c : new PluginManager<AnnotationType>(AnnotationType.class).getInterfaces() )
            System.out.println("\t" + c.getSimpleName());
        System.out.println();
    }

}