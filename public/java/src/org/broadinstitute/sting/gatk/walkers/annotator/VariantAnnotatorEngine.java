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
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;


public class VariantAnnotatorEngine {

    private List<InfoFieldAnnotation> requestedInfoAnnotations;
    private List<GenotypeAnnotation> requestedGenotypeAnnotations;
    private List<VAExpression> requestedExpressions = new ArrayList<VAExpression>();

    private final HashMap<RodBinding<VariantContext>, String> dbAnnotations = new HashMap<RodBinding<VariantContext>, String>();
    private final AnnotatorCompatibleWalker walker;
    private final GenomeAnalysisEngine toolkit;

    private boolean requireStrictAlleleMatch = false;

    protected static class VAExpression {

        public String fullName, fieldName;
        public RodBinding<VariantContext> binding;

        public VAExpression(String fullExpression, List<RodBinding<VariantContext>> bindings) {
            int indexOfDot = fullExpression.lastIndexOf(".");
            if ( indexOfDot == -1 )
                throw new UserException.BadArgumentValue(fullExpression, "it should be in rodname.value format");

            fullName = fullExpression;
            fieldName = fullExpression.substring(indexOfDot+1);

            String bindingName = fullExpression.substring(0, indexOfDot);
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

    // experimental constructor for active region traversal
    public VariantAnnotatorEngine(GenomeAnalysisEngine toolkit) {
        this.walker = null;
        this.toolkit = toolkit;
        requestedInfoAnnotations = AnnotationInterfaceManager.createInfoFieldAnnotations(Arrays.asList("ActiveRegionBasedAnnotation"), Collections.<String>emptyList());
    }

    // select specific expressions to use
    public void initializeExpressions(List<String> expressionsToUse) {
        // set up the expressions
        for ( String expression : expressionsToUse )
            requestedExpressions.add(new VAExpression(expression, walker.getResourceRodBindings()));
    }

    protected List<VAExpression> getRequestedExpressions() { return requestedExpressions; }

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

    public void setRequireStrictAlleleMatch( final boolean requireStrictAlleleMatch ) {
        this.requireStrictAlleleMatch = requireStrictAlleleMatch;
    }

    public VariantContext annotateContext(final RefMetaDataTracker tracker, final ReferenceContext ref, final Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        Map<String, Object> infoAnnotations = new LinkedHashMap<String, Object>(vc.getAttributes());

        // annotate db occurrences
        vc = annotateDBs(tracker, ref, vc, infoAnnotations);

        // annotate expressions where available
        annotateExpressions(tracker, ref, infoAnnotations);

        // go through all the requested info annotationTypes
        for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, walker, ref, stratifiedContexts, vc);
            if ( annotationsFromCurrentType != null )
                infoAnnotations.putAll(annotationsFromCurrentType);
        }

        // generate a new annotated VC
        VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(infoAnnotations);

        // annotate genotypes, creating another new VC in the process
        return builder.genotypes(annotateGenotypes(tracker, ref, stratifiedContexts, vc)).make();
    }

    public VariantContext annotateContext(final Map<String, Map<Allele, List<GATKSAMRecord>>> stratifiedContexts, VariantContext vc) {
        Map<String, Object> infoAnnotations = new LinkedHashMap<String, Object>(vc.getAttributes());

        // go through all the requested info annotationTypes
        for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            Map<String, Object> annotationsFromCurrentType = ((ActiveRegionBasedAnnotation)annotationType).annotate(stratifiedContexts, vc);
            if ( annotationsFromCurrentType != null ) {
                infoAnnotations.putAll(annotationsFromCurrentType);
            }
        }

        // generate a new annotated VC
        return new VariantContextBuilder(vc).attributes(infoAnnotations).make();
    }

    private VariantContext annotateDBs(RefMetaDataTracker tracker, ReferenceContext ref, VariantContext vc, Map<String, Object> infoAnnotations) {
        for ( Map.Entry<RodBinding<VariantContext>, String> dbSet : dbAnnotations.entrySet() ) {
            if ( dbSet.getValue().equals(VCFConstants.DBSNP_KEY) ) {
                final String rsID = VCFUtils.rsIDOfFirstRealVariant(tracker.getValues(dbSet.getKey(), ref.getLocus()), vc.getType());
                
                // add the ID if appropriate
                if ( rsID != null ) {
                    // put the DB key into the INFO field
                    infoAnnotations.put(VCFConstants.DBSNP_KEY, true);

                    if ( vc.emptyID() ) {
                        vc = new VariantContextBuilder(vc).id(rsID).make();
                    } else if ( walker.alwaysAppendDbsnpId() && vc.getID().indexOf(rsID) == -1 ) {
                        final String newRsID = vc.getID() + VCFConstants.ID_FIELD_SEPARATOR + rsID;
                        vc = new VariantContextBuilder(vc).id(newRsID).make();
                    }
                }
            } else {
                boolean overlapsComp = false;
                for ( VariantContext comp : tracker.getValues(dbSet.getKey(), ref.getLocus()) ) {
                    if ( !comp.isFiltered() && ( !requireStrictAlleleMatch || comp.getAlleles().equals(vc.getAlleles()) ) ) {
                        overlapsComp = true;
                        break;
                    }
                }
                if ( overlapsComp )
                    infoAnnotations.put(dbSet.getValue(), overlapsComp);
            }
        }

        return vc;
    }

    private void annotateExpressions(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, Object> infoAnnotations) {
        for ( VAExpression expression : requestedExpressions ) {
            Collection<VariantContext> VCs = tracker.getValues(expression.binding, ref.getLocus());
            if ( VCs.size() == 0 )
                continue;

            VariantContext vc = VCs.iterator().next();
            // special-case the ID field
            if ( expression.fieldName.equals("ID") ) {
                if ( vc.hasID() )
                    infoAnnotations.put(expression.fullName, vc.getID());
            } else if ( vc.hasAttribute(expression.fieldName) ) {
                infoAnnotations.put(expression.fullName, vc.getAttribute(expression.fieldName));
            }
        }
    }

    private GenotypesContext annotateGenotypes(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( requestedGenotypeAnnotations.isEmpty() )
            return vc.getGenotypes();

        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            AlignmentContext context = stratifiedContexts.get(genotype.getSampleName());

            if ( context == null ) {
                genotypes.add(genotype);
            } else {
                final GenotypeBuilder gb = new GenotypeBuilder(genotype);
                for ( final GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
                    annotation.annotate(tracker, walker, ref, context, vc, genotype, gb);
                }
                genotypes.add(gb.make());
            }
        }

        return genotypes;
    }
}
