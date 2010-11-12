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

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotationInterfaceManager;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;

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

    private List<InfoFieldAnnotation> requestedInfoAnnotations;
    private List<GenotypeAnnotation> requestedGenotypeAnnotations;


    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        List<String> annotationClasses = Arrays.asList(annotationClassesToUse);
        List<String> annotations = Arrays.asList(annotationsToUse);
        AnnotationInterfaceManager.validateAnnotations(annotationClasses, annotations);
        requestedInfoAnnotations = AnnotationInterfaceManager.createInfoFieldAnnotations(annotationClasses, annotations);
        requestedGenotypeAnnotations = AnnotationInterfaceManager.createGenotypeAnnotations(annotationClasses, annotations);
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
        } catch (Exception e) {
            throw new DynamicClassResolutionException(c, e);
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
                vc_eval = tracker.getVariantContext(ref,"eval", vc, loc, true);
                vc_ref = tracker.getVariantContext(ref,"HapMap", vc, loc, true);
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

                String key = annotationType.getKeyNames().get(0);
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