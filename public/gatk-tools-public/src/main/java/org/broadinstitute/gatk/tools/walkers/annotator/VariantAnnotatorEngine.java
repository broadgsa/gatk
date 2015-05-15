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

package org.broadinstitute.gatk.tools.walkers.annotator;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections.ListUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.*;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;


public class VariantAnnotatorEngine {
    private final static Logger logger = Logger.getLogger(VariantAnnotatorEngine.class);
    private List<InfoFieldAnnotation> requestedInfoAnnotations = Collections.emptyList();
    private List<GenotypeAnnotation> requestedGenotypeAnnotations = Collections.emptyList();
    private List<VAExpression> requestedExpressions = new ArrayList<>();

    private final AnnotatorCompatible walker;
    private final GenomeAnalysisEngine toolkit;

    VariantOverlapAnnotator variantOverlapAnnotator = null;

    // Map of info field name to info field
    private final Map<String, VCFInfoHeaderLine> hInfoMap = new HashMap<>();

    protected static class VAExpression {

        public String fullName, fieldName;
        public RodBinding<VariantContext> binding;

        public VAExpression(String fullExpression, List<RodBinding<VariantContext>> bindings){
            final int indexOfDot = fullExpression.lastIndexOf(".");
            if ( indexOfDot == -1 )
                throw new UserException.BadArgumentValue(fullExpression, "it should be in rodname.value format");

            fullName = fullExpression;
            fieldName = fullExpression.substring(indexOfDot+1);

            final String bindingName = fullExpression.substring(0, indexOfDot);
            for ( final RodBinding<VariantContext> rod : bindings ) {
                if ( rod.getName().equals(bindingName) ) {
                    binding = rod;
                    break;
                }
            }
        }
    }

    // use this constructor if you want all possible annotations
    public VariantAnnotatorEngine(List<String> annotationsToExclude, AnnotatorCompatible walker, GenomeAnalysisEngine toolkit) {
        this.walker = walker;
        this.toolkit = toolkit;
        requestedInfoAnnotations = AnnotationInterfaceManager.createAllInfoFieldAnnotations();
        requestedGenotypeAnnotations = AnnotationInterfaceManager.createAllGenotypeAnnotations();
        excludeAnnotations(annotationsToExclude);
        initializeDBs(toolkit);
    }

    // use this constructor if you want to select specific annotations (and/or interfaces)
    public VariantAnnotatorEngine(List<String> annotationGroupsToUse, List<String> annotationsToUse, List<String> annotationsToExclude, AnnotatorCompatible walker, GenomeAnalysisEngine toolkit) {
        this.walker = walker;
        this.toolkit = toolkit;
        initializeAnnotations(annotationGroupsToUse, annotationsToUse, annotationsToExclude);
        initializeDBs(toolkit);
    }

    public void makeHeaderInfoMap(final Set<VCFHeaderLine> hInfo ){
        for ( VCFHeaderLine hLine : hInfo ) {
            if ( hLine instanceof VCFInfoHeaderLine )
                hInfoMap.put( ((VCFInfoHeaderLine)hLine).getID(), (VCFInfoHeaderLine)hLine);
        }
    }

    // select specific expressions to use
    public void initializeExpressions(Set<String> expressionsToUse) {
        // set up the expressions
        for ( final String expression : expressionsToUse )
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
        if ( annotationsToExclude.isEmpty() )
            return;

        final List<InfoFieldAnnotation> tempRequestedInfoAnnotations = new ArrayList<>(requestedInfoAnnotations.size());
        for ( final InfoFieldAnnotation annotation : requestedInfoAnnotations ) {
            if ( !annotationsToExclude.contains(annotation.getClass().getSimpleName()) )
                tempRequestedInfoAnnotations.add(annotation);
        }
        requestedInfoAnnotations = tempRequestedInfoAnnotations;

        final List<GenotypeAnnotation> tempRequestedGenotypeAnnotations = new ArrayList<>(requestedGenotypeAnnotations.size());
        for ( final GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
            if ( !annotationsToExclude.contains(annotation.getClass().getSimpleName()) )
                tempRequestedGenotypeAnnotations.add(annotation);
        }
        requestedGenotypeAnnotations = tempRequestedGenotypeAnnotations;
    }

    private void initializeDBs(final GenomeAnalysisEngine engine) {
        // check to see whether comp rods were included
        RodBinding<VariantContext> dbSNPBinding = walker.getDbsnpRodBinding();
        if ( dbSNPBinding != null && ! dbSNPBinding.isBound() )
            dbSNPBinding = null;

        final Map<RodBinding<VariantContext>, String> overlapBindings = new LinkedHashMap<>();
        for ( final RodBinding<VariantContext> b : walker.getCompRodBindings())
            if ( b.isBound() ) overlapBindings.put(b, b.getName());
        if ( dbSNPBinding != null && ! overlapBindings.keySet().contains(VCFConstants.DBSNP_KEY) )
            overlapBindings.put(dbSNPBinding, VCFConstants.DBSNP_KEY); // add overlap detection with DBSNP by default

        variantOverlapAnnotator = new VariantOverlapAnnotator(dbSNPBinding, overlapBindings, engine.getGenomeLocParser());
    }

    public void invokeAnnotationInitializationMethods( final Set<VCFHeaderLine> headerLines ) {
        for ( final VariantAnnotatorAnnotation annotation : requestedInfoAnnotations ) {
            annotation.initialize(walker, toolkit, headerLines);
        }

        for ( final VariantAnnotatorAnnotation annotation : requestedGenotypeAnnotations ) {
            annotation.initialize(walker, toolkit, headerLines);
        }
    }

    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {
        final Set<VCFHeaderLine> descriptions = new HashSet<>();

        for ( final InfoFieldAnnotation annotation : requestedInfoAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( final GenotypeAnnotation annotation : requestedGenotypeAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( final String db : variantOverlapAnnotator.getOverlapNames() ) {
            if ( VCFStandardHeaderLines.getInfoLine(db, false) != null )
                descriptions.add(VCFStandardHeaderLines.getInfoLine(db));
            else
                descriptions.add(new VCFInfoHeaderLine(db, 0, VCFHeaderLineType.Flag, db + " Membership"));
        }

        return descriptions;
    }

    public VariantContext annotateContext(final RefMetaDataTracker tracker,
                                          final ReferenceContext ref,
                                          final Map<String, AlignmentContext> stratifiedContexts,
                                          final VariantContext vc) {
        return annotateContext(tracker, ref, stratifiedContexts, vc, null);
    }

    public VariantContext annotateContext(final RefMetaDataTracker tracker,
                                          final ReferenceContext ref,
                                          final Map<String, AlignmentContext> stratifiedContexts,
                                          final VariantContext vc,
                                          final Map<String,PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final Map<String, Object> infoAnnotations = new LinkedHashMap<>(vc.getAttributes());

        // annotate expressions where available
        annotateExpressions(tracker, ref.getLocus(), vc, infoAnnotations);

        // go through all the requested info annotationTypes
        for ( final InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, walker, ref, stratifiedContexts, vc, perReadAlleleLikelihoodMap);
            if ( annotationsFromCurrentType != null )
                infoAnnotations.putAll(annotationsFromCurrentType);
        }

        // generate a new annotated VC
        final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(infoAnnotations);

        // annotate genotypes, creating another new VC in the process
        final VariantContext annotated = builder.genotypes(annotateGenotypes(tracker, ref, stratifiedContexts, vc, perReadAlleleLikelihoodMap)).make();

        // annotate db occurrences
        return annotateDBs(tracker, annotated);
    }

    public VariantContext annotateContextForActiveRegion(final ReferenceContext referenceContext,
                                                         final RefMetaDataTracker tracker,
                                                         final ReadLikelihoods<Allele> readLikelihoods,
                                                         final VariantContext vc) {
        //TODO we transform the read-likelihood into the Map^2 previous version for the sake of not changing of not changing annotation interface.
        //TODO should we change those interfaces?

        final Map<String, PerReadAlleleLikelihoodMap> annotationLikelihoods = readLikelihoods.toPerReadAlleleLikelihoodMap();
        return annotateContextForActiveRegion(referenceContext, tracker, annotationLikelihoods, vc);
    }

    public VariantContext annotateContextForActiveRegion(final ReferenceContext referenceContext,
                                                         final RefMetaDataTracker tracker,
                                                         final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                                         final VariantContext vc) {
        final Map<String, Object> infoAnnotations = new LinkedHashMap<>(vc.getAttributes());

        // go through all the requested info annotationTypes
        for ( final InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            if ( !(annotationType instanceof ActiveRegionBasedAnnotation) )
                continue;

            final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(referenceContext, perReadAlleleLikelihoodMap, vc);
            if ( annotationsFromCurrentType != null ) {
                infoAnnotations.putAll(annotationsFromCurrentType);
            }
        }

        // generate a new annotated VC
        final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(infoAnnotations);

        // annotate genotypes, creating another new VC in the process
        final VariantContext annotated = builder.genotypes(annotateGenotypes(null, null, null, vc, perReadAlleleLikelihoodMap)).make();

        // annotate db occurrences
        return annotateDBs(tracker, annotated);
    }

    /**
     * Annotate the ID field and other DBs for the given Variant Context
     *
     * @param tracker  ref meta data tracker (cannot be null)
     * @param vc       variant context to annotate
     * @return non-null annotated version of vc
     */
    @Requires({"tracker != null && loc != null && vc != null && infoAnnotations != null"})
    @Ensures("result != null")
    private VariantContext annotateDBs(final RefMetaDataTracker tracker, VariantContext vc) {
        return variantOverlapAnnotator.annotateOverlaps(tracker, variantOverlapAnnotator.annotateRsID(tracker, vc));
    }

    /**
     * Annotate the requested expressions
     *
     * @param tracker   ref meta data tracker (cannot be null)
     * @param loc       the location on the genome
     * @param vc        variant context to annotate
     * @param infoAnnotations the annotations for the requested expressions
     */
    @Requires({"tracker != null && loc != null && vc != null"})
    private void annotateExpressions(final RefMetaDataTracker tracker, final GenomeLoc loc, final VariantContext vc, final Map<String, Object> infoAnnotations){

        // each requested expression
        for ( final VAExpression expression : requestedExpressions ) {

            // get the variant contexts for all the expressions at the location
            final Collection<VariantContext> expressionVCs = tracker.getValues(expression.binding, loc);
            if ( expressionVCs.isEmpty() )
                continue;

            // get the expression's variant context
            final VariantContext expressionVC = expressionVCs.iterator().next();

            // special-case the ID field
            if ( expression.fieldName.equals("ID") ) {
                if ( expressionVC.hasID() )
                    infoAnnotations.put(expression.fullName, expressionVC.getID());
            } else if (expression.fieldName.equals("ALT")) {
                infoAnnotations.put(expression.fullName, expressionVC.getAlternateAllele(0).getDisplayString());
            } else if ( expressionVC.hasAttribute(expression.fieldName) ) {
                // find the info field
                final VCFInfoHeaderLine  hInfo = hInfoMap.get(expression.fullName);
                if ( hInfo == null ){
                    throw new UserException("Cannot annotate expression " + expression.fullName + " at " + loc + " for variant allele(s) " + vc.getAlleles() + ", missing header info");
                }

                // can not annotate if more variant than expression alleles
                if ( expressionVC.getNAlleles() < vc.getNAlleles() ) {
                    logger.warn("Skipping expression " + expression.fullName + " at " + loc + ", can not match " + expressionVC.getNAlleles() + " in the expression to " +
                            vc.getNAlleles() + " in the variant");
                    continue;
                }

                //
                // Add the info field annotations
                //

                final boolean isMultiAllelic = expressionVC.getNAlleles() > 2;
                final boolean useRefAndAltAlleles = VCFHeaderLineCount.R == hInfo.getCountType();
                final boolean useAltAlleles = VCFHeaderLineCount.A == hInfo.getCountType();
                List<Allele> usedExpressionAlleles = null;

                // Multiallelic and count of A or R
                if ( isMultiAllelic && (useAltAlleles || useRefAndAltAlleles) ){

                    // remove brackets and spaces from expression attribute
                    final String cleanedExpression = expressionVC.getAttribute(expression.fieldName).toString().replaceAll("[\\[\\]\\s]", "");

                    // map where key = expression allele string value = expression value corresponding to the allele
                    final Map<String, String> mapAlleleToExpressionValue = new HashMap<String, String>();

                    // get comma separated expression values
                    ArrayList<String> expressionValuesList = new ArrayList<String>(Arrays.asList(cleanedExpression.split(",")));

                    if ( vc.isSNP() && expressionVC.isMixed() ){
                        final VariantContextBuilder builder = new VariantContextBuilder(expressionVC);
                        List<Allele> sameLengthAlleles = new ArrayList<Allele>();

                        // get alt alleles that are the same length as the ref allele
                        Iterator<String> expressionValuesIterator = expressionValuesList.iterator();
                        for ( Allele allele : expressionVC.getAlleles() ){
                            if ( allele.isNonReference() ){
                                if ( !expressionValuesIterator.hasNext() ){
                                    logger.warn("Cannot annotate expression " + expression.fullName + " at " + loc + " for expression allele): " + allele);
                                    break;
                                }
                                expressionValuesIterator.next();
                                if ( allele.length() == expressionVC.getReference().length() ) {
                                    sameLengthAlleles.add(allele);
                                }
                                else {
                                    // remove unused expression values
                                    expressionValuesIterator.remove();
                                }
                            } else {
                                if ( useRefAndAltAlleles )
                                    expressionValuesIterator.remove();
                            }
                        }

                        if (!sameLengthAlleles.isEmpty()) {
                            sameLengthAlleles.add(0, expressionVC.getReference());
                            VariantContext variantContext = builder.alleles(sameLengthAlleles).make();
                            // extract the SNPs
                            VariantContext variantContextTrimmed = GATKVariantContextUtils.trimAlleles(variantContext, true, true);
                            usedExpressionAlleles = useRefAndAltAlleles ? variantContextTrimmed.getAlleles() : variantContextTrimmed.getAlternateAlleles();
                        }
                    } else {
                        // get the alleles common to the expression and variant
                        usedExpressionAlleles = useRefAndAltAlleles ? expressionVC.getAlleles() : expressionVC.getAlternateAlleles();
                    }

                    final List<Allele> commonAlleles = ListUtils.intersection(usedExpressionAlleles, vc.getAlleles());

                    // the number of expression values must be the same as the number of alleles
                    if ( expressionValuesList.size() != usedExpressionAlleles.size() ) {
                        logger.warn("Cannot annotate expression " + expression.fullName + " at " + loc + " for variant allele(s): " + vc.getAlleles() + ", " +
                                expressionValuesList.size() + " expression values is not equal to " + usedExpressionAlleles.size() + " expression alleles");
                        continue;
                    }

                    // map the used expression alleles to it's value
                    for (int i = 0; i != expressionValuesList.size(); i++)
                        mapAlleleToExpressionValue.put(usedExpressionAlleles.get(i).getBaseString(), expressionValuesList.get(i));

                    // add the variants expression values to the annotation
                    final List<String> annotationValues = new ArrayList<String>();
                    for (final Allele commonAllele : commonAlleles) {
                        annotationValues.add(mapAlleleToExpressionValue.get(commonAllele.getBaseString()));
                    }

                    infoAnnotations.put(expression.fullName, annotationValues);
                } else {
                    infoAnnotations.put(expression.fullName, expressionVC.getAttribute(expression.fieldName));
                }
            }
        }
    }

    private GenotypesContext annotateGenotypes(final RefMetaDataTracker tracker,
                                               final ReferenceContext ref, final Map<String, AlignmentContext> stratifiedContexts,
                                               final VariantContext vc,
                                               final Map<String,PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( requestedGenotypeAnnotations.isEmpty() )
            return vc.getGenotypes();

        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            AlignmentContext context = null;
            PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = null;
            if (stratifiedContexts != null)
                context = stratifiedContexts.get(genotype.getSampleName());
            if (stratifiedPerReadAlleleLikelihoodMap != null)
                perReadAlleleLikelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(genotype.getSampleName());


            final GenotypeBuilder gb = new GenotypeBuilder(genotype);
            for ( final GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
                annotation.annotate(tracker, walker, ref, context, vc, genotype, gb, perReadAlleleLikelihoodMap);
            }
            genotypes.add(gb.make());
        }

        return genotypes;
    }
}
