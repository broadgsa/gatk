/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotationType;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.playground.gatk.walkers.annotator.GenomicAnnotation;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.classloader.PackageUtils;

import java.util.*;

public class AnnotationByAlleleFrequencyWalker  extends RodWalker<Integer, Integer> {

    /* Get selected annotations from a given VCF File. If these annotations match a site reference ROD, output the annotation value.
    Usage example:
        java -jar dist/GenomeAnalysisTK.jar -R reference.fasta \
        -T AnnotationByAlleleFrequency \
        -A QualByDepth  \
        -B eval,VCF,eval.vcf \
        -B HapMap,VCF,ref.vcf


     */
    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////

    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected String[] annotationsToUse = {};

    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { };

    private ArrayList<InfoFieldAnnotation> requestedInfoAnnotations;
    private ArrayList<GenotypeAnnotation> requestedGenotypeAnnotations;


    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

        // create a map for all annotation classes which implement our top-level interfaces
        HashMap<String, Class> classMap = new HashMap<String, Class>();
        for ( Class c : PackageUtils.getClassesImplementingInterface(InfoFieldAnnotation.class) )
            classMap.put(c.getSimpleName(), c);
        for ( Class c : PackageUtils.getClassesImplementingInterface(GenotypeAnnotation.class) )
            classMap.put(c.getSimpleName(), c);
        for ( Class c : PackageUtils.getInterfacesExtendingInterface(AnnotationType.class) )
            classMap.put(c.getSimpleName(), c);

        HashSet<Class> classes = new HashSet<Class>();
        // get the classes from the provided groups (interfaces)
        for ( String group : annotationClassesToUse ) {
            Class interfaceClass = classMap.get(group);
            if ( interfaceClass == null )
                interfaceClass = classMap.get(group + "Annotation");
            if ( interfaceClass == null )
                throw new StingException("Class " + group + " is not found; please check that you have specified the class name correctly");
            classes.addAll(PackageUtils.getClassesImplementingInterface(interfaceClass));
        }

        // get the specific classes provided
        for ( String annotation : annotationsToUse ) {
            Class annotationClass = classMap.get(annotation);
            if ( annotationClass == null )
                annotationClass = classMap.get(annotation + "Annotation");
            if ( annotationClass == null )
                throw new StingException("Class " + annotation + " is not found; please check that you have specified the class name correctly");
            classes.add(annotationClass);
        }

        // get the instances
        requestedInfoAnnotations = new ArrayList<InfoFieldAnnotation>();
        requestedGenotypeAnnotations = new ArrayList<GenotypeAnnotation>();

        for ( Class c : classes ) {
            // note that technically an annotation can work on both the INFO and FORMAT fields
            if ( InfoFieldAnnotation.class.isAssignableFrom(c) )
                requestedInfoAnnotations.add((InfoFieldAnnotation)getInstance(c));
            if ( GenotypeAnnotation.class.isAssignableFrom(c) )
                requestedGenotypeAnnotations.add((GenotypeAnnotation)getInstance(c));
        }

    }

    private static <T> ArrayList<T> getInstances(List<Class<? extends T>> classes) {
        ArrayList<T> objects = new ArrayList<T>();
        for ( Class c : classes )
            objects.add((T)getInstance(c));
        return objects;
    }

    private static <T> T getInstance(Class<T> c) {
        try {
            return c.newInstance();
        } catch (InstantiationException e) {
            throw new StingException(String.format("Cannot instantiate annotation class '%s': must be concrete class", c.getSimpleName()));
        } catch (IllegalAccessException e) {
            throw new StingException(String.format("Cannot instantiate annotation class '%s': must have no-arg constructor", c.getSimpleName()));
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------


    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if( tracker != null ) {
            EnumSet<VariantContext.Type> vc = EnumSet.of(VariantContext.Type.SNP);
            GenomeLoc loc = context.getLocation();
            VariantContext vc_eval;
            VariantContext vc_ref;


            try {
                vc_eval = tracker.getVariantContext("eval", vc, loc, true);
                vc_ref = tracker.getVariantContext("HapMap", vc, loc, true);
            } catch (java.util.NoSuchElementException e) {
                return 0;
            }

            if (vc_ref == null || vc_eval == null)  {
                return 0;
            }

            // Get Allele frequency for reference ROD
            double af_ref = Double.valueOf(2*vc_ref.getHomVarCount()+vc_ref.getHetCount());
            af_ref = af_ref / (vc_ref.getChromosomeCount());
            System.out.format("AF_Ref: %5.4f ", af_ref);

            // Get Allele frequency for eval ROD
            double af_eval = Double.valueOf(2*vc_eval.getHomVarCount()+vc_eval.getHetCount());
            af_eval = af_eval / (vc_eval.getChromosomeCount());
            System.out.format("AF_Eval: %5.4f ", af_eval);


/*
            String qq = vc_eval.getAttributeAsString("AF");
            System.out.format("AF_EvalVCF: %s ", qq);
            */
            
            //go through all the requested info annotationTypes
            for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations )
            {

                String key = annotationType.getKeyName();
                String value_str = vc_eval.getAttributeAsString(key);
                System.out.format("%s: %s ", key, value_str);
            }
            System.out.println();

        }

        return 1; // This value isn't actually used for anything
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer reduceInit() {
        return 0; // Nothing to do here
    }

    public Integer reduce( Integer value, Integer sum ) {
        return 0; // Nothing to do here
    }

    public void onTraversalDone( Integer sum ) {

   }
}