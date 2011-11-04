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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;


public class VariantAnnotatorEngine {

    private List<InfoFieldAnnotation> requestedInfoAnnotations;
    private List<GenotypeAnnotation> requestedGenotypeAnnotations;
    private List<VAExpression> requestedExpressions = new ArrayList<VAExpression>();

    private HashMap<RodBinding<VariantContext>, String> dbAnnotations = new HashMap<RodBinding<VariantContext>, String>();
    private AnnotatorCompatibleWalker walker;
    private GenomeAnalysisEngine toolkit;

    private static class VAExpression {

        public String fullName, fieldName;
        public RodBinding<VariantContext> binding;

        public VAExpression(String fullEpression, List<RodBinding<VariantContext>> bindings) {
            int indexOfDot = fullEpression.lastIndexOf(".");
            if ( indexOfDot == -1 )
                throw new UserException.BadArgumentValue(fullEpression, "it should be in rodname.value format");

            fullName = fullEpression;
            fieldName = fullEpression.substring(indexOfDot+1);

            String bindingName = fullEpression.substring(0, indexOfDot);
            for ( RodBinding<VariantContext> rod : bindings ) {
                if ( rod.getName().equals(bindingName) ) {
                    binding = rod;
                    break;
                }
            }
        }
    }

    // use this constructor if you want all possible annotations
    public VariantAnnotatorEngine(List<String> annotationsToExclude, AnnotatorCompatibleWalker walker, GenomeAnalysisEngine toolkit) {
        this.walker = walker;
        this.toolkit = toolkit;
        requestedInfoAnnotations = AnnotationInterfaceManager.createAllInfoFieldAnnotations();
        requestedGenotypeAnnotations = AnnotationInterfaceManager.createAllGenotypeAnnotations();
        excludeAnnotations(annotationsToExclude);
        initializeDBs();
    }

    // use this constructor if you want to select specific annotations (and/or interfaces)
    public VariantAnnotatorEngine(List<String> annotationGroupsToUse, List<String> annotationsToUse, List<String> annotationsToExclude, AnnotatorCompatibleWalker walker, GenomeAnalysisEngine toolkit) {
        this.walker = walker;
        this.toolkit = toolkit;
        initializeAnnotations(annotationGroupsToUse, annotationsToUse, annotationsToExclude);
        initializeDBs();
    }

    // select specific expressions to use
    public void initializeExpressions(List<String> expressionsToUse) {
        // set up the expressions
        for ( String expression : expressionsToUse )
            requestedExpressions.add(new VAExpression(expression, walker.getResourceRodBindings()));
    }

    private void initializeAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse, List<String> annotationsToExclude) {
        AnnotationInterfaceManager.validateAnnotations(annotationGroupsToUse, annotationsToUse);
        requestedInfoAnnotations = AnnotationInterfaceManager.createInfoFieldAnnotations(annotationGroupsToUse, annotationsToUse);
        requestedGenotypeAnnotations = AnnotationInterfaceManager.createGenotypeAnnotations(annotationGroupsToUse, annotationsToUse);
        excludeAnnotations(annotationsToExclude);
    }

    private void excludeAnnotations(List<String> annotationsToExclude) {
        if ( annotationsToExclude.size() == 0 )
            return;

        List<InfoFieldAnnotation> tempRequestedInfoAnnotations = new ArrayList<InfoFieldAnnotation>(requestedInfoAnnotations.size());
        for ( InfoFieldAnnotation annotation : requestedInfoAnnotations ) {
            if ( !annotationsToExclude.contains(annotation.getClass().getSimpleName()) )
                tempRequestedInfoAnnotations.add(annotation);
        }
        requestedInfoAnnotations = tempRequestedInfoAnnotations;

        List<GenotypeAnnotation> tempRequestedGenotypeAnnotations = new ArrayList<GenotypeAnnotation>(requestedGenotypeAnnotations.size());
        for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
            if ( !annotationsToExclude.contains(annotation.getClass().getSimpleName()) )
                tempRequestedGenotypeAnnotations.add(annotation);
        }
        requestedGenotypeAnnotations = tempRequestedGenotypeAnnotations;
    }

    private void initializeDBs() {

        // check to see whether comp rods were included
        final RodBinding<VariantContext> dbsnp = walker.getDbsnpRodBinding();
        if ( dbsnp != null &&  dbsnp.isBound() )
            dbAnnotations.put(dbsnp, VCFConstants.DBSNP_KEY);

        final List<RodBinding<VariantContext>> comps = walker.getCompRodBindings();
        for ( RodBinding<VariantContext> rod : comps )
            dbAnnotations.put(rod, rod.getName());
    }

    public void invokeAnnotationInitializationMethods( Set<VCFHeaderLine> headerLines ) {
        for ( VariantAnnotatorAnnotation annotation : requestedInfoAnnotations ) {
            annotation.initialize(walker, toolkit, headerLines);
        }

        for ( VariantAnnotatorAnnotation annotation : requestedGenotypeAnnotations ) {
            annotation.initialize(walker, toolkit, headerLines);
        }
    }

    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {

        Set<VCFHeaderLine> descriptions = new HashSet<VCFHeaderLine>();

        for ( InfoFieldAnnotation annotation : requestedInfoAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( String db : dbAnnotations.values() )
            descriptions.add(new VCFInfoHeaderLine(db, 0, VCFHeaderLineType.Flag, (db.equals(VCFConstants.DBSNP_KEY) ? "dbSNP" : db) + " Membership"));

        return descriptions;
    }

    public VariantContext annotateContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {

        Map<String, Object> infoAnnotations = new LinkedHashMap<String, Object>(vc.getAttributes());

        // annotate db occurrences
        annotateDBs(tracker, ref, vc, infoAnnotations);

        // annotate expressions where available
        annotateExpressions(tracker, ref, infoAnnotations);

        // go through all the requested info annotationTypes
        for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, walker, ref, stratifiedContexts, vc);
            if ( annotationsFromCurrentType != null )
                infoAnnotations.putAll(annotationsFromCurrentType);
        }

        // generate a new annotated VC
        final VariantContext annotatedVC = VariantContext.modifyAttributes(vc, infoAnnotations);

        // annotate genotypes, creating another new VC in the process
        return VariantContext.modifyGenotypes(annotatedVC, annotateGenotypes(tracker, ref, stratifiedContexts, vc));
    }

    private void annotateDBs(RefMetaDataTracker tracker, ReferenceContext ref, VariantContext vc, Map<String, Object> infoAnnotations) {
        for ( Map.Entry<RodBinding<VariantContext>, String> dbSet : dbAnnotations.entrySet() ) {
            if ( dbSet.getValue().equals(VCFConstants.DBSNP_KEY) ) {
                String rsID = VCFUtils.rsIDOfFirstRealVariant(tracker.getValues(dbSet.getKey(), ref.getLocus()), vc.getType());
                infoAnnotations.put(VCFConstants.DBSNP_KEY, rsID != null);
                // annotate dbsnp id if available and not already there
                if ( rsID != null && (!vc.hasID() || vc.getID().equals(VCFConstants.EMPTY_ID_FIELD)) )
                    infoAnnotations.put(VariantContext.ID_KEY, rsID);
            } else {
                boolean overlapsComp = false;
                for ( VariantContext comp : tracker.getValues(dbSet.getKey(), ref.getLocus()) ) {
                    if ( !comp.isFiltered() ) {
                        overlapsComp = true;
                        break;
                    }
                }
                infoAnnotations.put(dbSet.getValue(), overlapsComp);
            }
        }
    }

    private void annotateExpressions(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, Object> infoAnnotations) {
        for ( VAExpression expression : requestedExpressions ) {
            Collection<VariantContext> VCs = tracker.getValues(expression.binding, ref.getLocus());
            if ( VCs.size() == 0 )
                continue;

            VariantContext vc = VCs.iterator().next();
            if ( vc.hasAttribute(expression.fieldName) )
                infoAnnotations.put(expression.fullName, vc.getAttribute(expression.fieldName));
        }
    }

    private Map<String, Genotype> annotateGenotypes(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( requestedGenotypeAnnotations.size() == 0 )
            return vc.getGenotypes();

        Map<String, Genotype> genotypes = new HashMap<String, Genotype>(vc.getNSamples());
        for ( Map.Entry<String, Genotype> g : vc.getGenotypes().entrySet() ) {
            Genotype genotype = g.getValue();
            AlignmentContext context = stratifiedContexts.get(g.getKey());
            if ( context == null ) {
                genotypes.put(g.getKey(), genotype);
                continue;
            }

            Map<String, Object> genotypeAnnotations = new HashMap<String, Object>(genotype.getAttributes());
            for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
                Map<String, Object> result = annotation.annotate(tracker, walker, ref, context, vc, genotype);
                if ( result != null )
                    genotypeAnnotations.putAll(result);
            }
            genotypes.put(g.getKey(), new Genotype(g.getKey(), genotype.getAlleles(), genotype.getNegLog10PError(), genotype.getFilters(), genotypeAnnotations, genotype.isPhased()));
        }

        return genotypes;
    }
}
