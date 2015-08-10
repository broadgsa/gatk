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

import org.broadinstitute.gatk.utils.classloader.PluginManager;

import java.util.List;

public class AnnotationHelpUtils {

    /**
     * Simple method to print a list of available annotations.
     */
    public static void listAnnotations() {
        System.out.println("\nThis is a list of available Variant Annotations for use with tools such as UnifiedGenotyper, HaplotypeCaller and VariantAnnotator. Please see the Technical Documentation for more details about these annotations:");
        System.out.println("http://www.broadinstitute.org/gatk/tooldocs/");
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
