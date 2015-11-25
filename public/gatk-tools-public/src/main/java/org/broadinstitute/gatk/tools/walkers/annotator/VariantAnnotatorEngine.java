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

package org.broadinstitute.gatk.tools.walkers.annotator;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.*;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;


public class VariantAnnotatorEngine {
    private final static Logger logger = Logger.getLogger(VariantAnnotatorEngine.class);
    private List<InfoFieldAnnotation> requestedInfoAnnotations = Collections.emptyList();
    private List<InfoFieldAnnotation> requestedReducibleInfoAnnotations = new ArrayList<>();
    private List<InfoFieldAnnotation> requestedNonReducibleInfoAnnotations = new ArrayList<>();
    private List<GenotypeAnnotation> requestedGenotypeAnnotations = Collections.emptyList();
    private List<VAExpression> requestedExpressions = new ArrayList<>();
    private boolean expressionAlleleConcordance = false;

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
        setReducibleAnnotations();
        initializeDBs(toolkit);
    }

    // use this constructor if you want to select specific annotations (and/or interfaces)
    public VariantAnnotatorEngine(List<String> annotationGroupsToUse, List<String> annotationsToUse, List<String> annotationsToExclude, AnnotatorCompatible walker, GenomeAnalysisEngine toolkit) {
        this.walker = walker;
        this.toolkit = toolkit;
        initializeAnnotations(annotationGroupsToUse, annotationsToUse, annotationsToExclude);
        setReducibleAnnotations();
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

    // set whether enforing allele concordance for expression
    public void setExpressionAlleleConcordance(Boolean expressionAlleleConcordance){
        this.expressionAlleleConcordance = expressionAlleleConcordance;
    }

    protected List<VAExpression> getRequestedExpressions() { return requestedExpressions; }

    public List<InfoFieldAnnotation> getRequestedReducibleInfoAnnotations() { return Collections.unmodifiableList(requestedReducibleInfoAnnotations); }

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
        // annotate genotypes
        final VariantContextBuilder builder = new VariantContextBuilder(vc).genotypes(annotateGenotypes(tracker, ref, stratifiedContexts, vc, perReadAlleleLikelihoodMap));
        VariantContext newGenotypeAnnotatedVC = builder.make();

        // annotate expressions where available
        final Map<String, Object> infoAnnotations = new LinkedHashMap<>(newGenotypeAnnotatedVC.getAttributes());
        annotateExpressions(tracker, ref.getLocus(), newGenotypeAnnotatedVC, infoAnnotations);

        // go through all the requested info annotationTypes
        for ( final InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, walker, ref, stratifiedContexts, newGenotypeAnnotatedVC, perReadAlleleLikelihoodMap);
            if ( annotationsFromCurrentType != null )
                infoAnnotations.putAll(annotationsFromCurrentType);
        }

        // create a new VC in the with info and genotype annotations
        final VariantContext annotated = builder.attributes(infoAnnotations).make();

        // annotate db occurrences
        return annotateDBs(tracker, annotated);
    }

    /**
     *
     * @param referenceContext
     * @param tracker
     * @param readLikelihoods
     * @param vc
     * @param useRaw    output annotation data as raw data? (Yes in the case of gVCF mode for HaplotypeCaller)
     * @return
     */
    public VariantContext annotateContextForActiveRegion(final ReferenceContext referenceContext,
                                                         final RefMetaDataTracker tracker,
                                                         final ReadLikelihoods<Allele> readLikelihoods,
                                                         final VariantContext vc,
                                                         final boolean useRaw) {
        //TODO we transform the read-likelihood into the Map^2 previous version for the sake of not changing of not changing annotation interface.
        //TODO should we change those interfaces?

        final Map<String, PerReadAlleleLikelihoodMap> annotationLikelihoods = readLikelihoods.toPerReadAlleleLikelihoodMap();
        return annotateContextForActiveRegion(referenceContext, tracker, annotationLikelihoods, vc, useRaw);
    }

    /**
     *
     * @param referenceContext
     * @param tracker
     * @param perReadAlleleLikelihoodMap
     * @param vc
     * @param useRaw    output annotation data as raw data? (Yes in the case of gVCF mode for HaplotypeCaller)
     * @return
     */
    public VariantContext annotateContextForActiveRegion(final ReferenceContext referenceContext,
                                                         final RefMetaDataTracker tracker,
                                                         final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                                         final VariantContext vc,
                                                         final boolean useRaw) {
        // annotate genotypes
        final VariantContextBuilder builder = new VariantContextBuilder(vc).genotypes(annotateGenotypes(null, null, null, vc, perReadAlleleLikelihoodMap));
        VariantContext newGenotypeAnnotatedVC = builder.make();

        final Map<String, Object> infoAnnotations = new LinkedHashMap<>(newGenotypeAnnotatedVC.getAttributes());

        // go through all the requested info annotationTypes that are reducible
        if (useRaw) {
            for (final InfoFieldAnnotation annotationType : requestedReducibleInfoAnnotations) {
                if (!(annotationType instanceof ActiveRegionBasedAnnotation))
                    continue;


                    ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;
                    final Map<String, Object> annotationsFromCurrentType = currentASannotation.annotateRawData(null, null, referenceContext, null, newGenotypeAnnotatedVC, perReadAlleleLikelihoodMap);
                    if (annotationsFromCurrentType != null) {
                        infoAnnotations.putAll(annotationsFromCurrentType);
                    }
            }
        }
        //if not in reference-confidence mode, do annotate with reducible annotations, but skip the raw data and go straight to the finalized values
        else {
            for (final InfoFieldAnnotation annotationType : requestedReducibleInfoAnnotations) {
                if (!(annotationType instanceof ActiveRegionBasedAnnotation))
                    continue;

                final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(null, null, referenceContext, null, newGenotypeAnnotatedVC, perReadAlleleLikelihoodMap);
                if (annotationsFromCurrentType != null) {
                    infoAnnotations.putAll(annotationsFromCurrentType);
                }
            }
        }
        //leave this in or else the median will overwrite until we do truly allele-specific
        //// for now do both allele-specific and not
        for ( final InfoFieldAnnotation annotationType : requestedNonReducibleInfoAnnotations ) {
            if ( !(annotationType instanceof ActiveRegionBasedAnnotation) )
                continue;

                final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(referenceContext, perReadAlleleLikelihoodMap, newGenotypeAnnotatedVC);
                if (annotationsFromCurrentType != null) {
                    infoAnnotations.putAll(annotationsFromCurrentType);
                }
        }

        // create a new VC with info and genotype annotations
        final VariantContext annotated = builder.attributes(infoAnnotations).make();

        // annotate db occurrences
        return annotateDBs(tracker, annotated);
    }

    /**
     * Combine (raw) data for reducible annotations (those that use raw data in gVCFs)
     * Mutates annotationMap by removing the annotations that were combined
     * @param allelesList   the list of merged alleles across all variants being combined
     * @param annotationMap attributes of merged variant contexts -- is modifying by removing successfully combined annotations
     * @return  a map containing the keys and raw values for the combined annotations
     */
    public Map<String, Object> combineAnnotations(final List<Allele> allelesList, Map<String, List<ReducibleAnnotationData>> annotationMap) {
        Map<String, Object> combinedAnnotations = new HashMap<>();

        // go through all the requested reducible info annotationTypes
        for (final InfoFieldAnnotation annotationType : requestedReducibleInfoAnnotations) {
                ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;
                if (annotationMap.containsKey(currentASannotation.getRawKeyName())) {
                    final List<ReducibleAnnotationData> annotationValue = annotationMap.get(currentASannotation.getRawKeyName());
                    final Map<String, Object> annotationsFromCurrentType = currentASannotation.combineRawData(allelesList, annotationValue);
                    combinedAnnotations.putAll(annotationsFromCurrentType);
                    //remove the combined annotations so that the next method only processes the non-reducible ones
                    annotationMap.remove(currentASannotation.getRawKeyName());
                }
        }
        return combinedAnnotations;
    }

    /**
     * Finalize reducible annotations (those that use raw data in gVCFs)
     * @param vc    the merged VC with the final set of alleles, possibly subset to the number of maxAltAlleles for genotyping
     * @param originalVC    the merged but non-subset VC that contains the full list of merged alleles
     * @return  a VariantContext with the final annotation values for reducible annotations
     */
    public VariantContext finalizeAnnotations(VariantContext vc, VariantContext originalVC) {
        final Map<String, Object> infoAnnotations = new LinkedHashMap<>(vc.getAttributes());

        // go through all the requested info annotationTypes
        for ( final InfoFieldAnnotation annotationType : requestedReducibleInfoAnnotations ) {

            ReducibleAnnotation currentASannotation = (ReducibleAnnotation)annotationType;

            final Map<String, Object> annotationsFromCurrentType = currentASannotation.finalizeRawData(vc, originalVC);
            if ( annotationsFromCurrentType != null ) {
                infoAnnotations.putAll(annotationsFromCurrentType);
                //clean up raw annotation data after annotations are finalized
                infoAnnotations.remove(currentASannotation.getRawKeyName());
            }
        }

        // generate a new annotated VC
        final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(infoAnnotations);

        // annotate genotypes, creating another new VC in the process
        final VariantContext annotated = builder.make();
        return annotated;
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
            } else if (expression.fieldName.equals("FILTER")) {
                if ( expressionVC.isFiltered() ) {
                    infoAnnotations.put(expression.fullName, expressionVC.getFilters().toString().replace("[", "").replace("]", "").replace(" ", ""));
                } else {
                    infoAnnotations.put(expression.fullName, "PASS");
                }
            } else if ( expressionVC.hasAttribute(expression.fieldName) ) {
                // find the info field
                final VCFInfoHeaderLine hInfo = hInfoMap.get(expression.fullName);
                if ( hInfo == null ){
                    throw new UserException("Cannot annotate expression " + expression.fullName + " at " + loc + " for variant allele(s) " + vc.getAlleles() + ", missing header info");
                }

                //
                // Add the info field annotations
                //
                final boolean useRefAndAltAlleles = VCFHeaderLineCount.R == hInfo.getCountType();
                final boolean useAltAlleles = VCFHeaderLineCount.A == hInfo.getCountType();

                // Annotation uses ref and/or alt alleles or enforce allele concordance
                if ( (useAltAlleles || useRefAndAltAlleles) || expressionAlleleConcordance ){

                    // remove brackets and spaces from expression value
                    final String cleanedExpressionValue = expressionVC.getAttribute(expression.fieldName).toString().replaceAll("[\\[\\]\\s]", "");

                    // get comma separated expression values
                    final ArrayList<String> expressionValuesList = new ArrayList<String>(Arrays.asList(cleanedExpressionValue.split(",")));

                    // get the minimum biallelics without genotypes
                    final List<VariantContext> minBiallelicVCs = getMinRepresentationBiallelics(vc);
                    final List<VariantContext> minBiallelicExprVCs = getMinRepresentationBiallelics(expressionVC);

                    // check concordance
                    final List<String> annotationValues = new ArrayList<>();
                    boolean canAnnotate = false;
                    for ( final VariantContext biallelicVC : minBiallelicVCs ) {
                        // check that ref and alt alleles are the same
                        List<Allele> exprAlleles = biallelicVC.getAlleles();
                        boolean isAlleleConcordant = false;
                        int i = 0;
                        for ( final VariantContext biallelicExprVC : minBiallelicExprVCs ){
                            List<Allele> alleles = biallelicExprVC.getAlleles();
                            // concordant
                            if ( alleles.equals(exprAlleles) ){
                                // get the value for the reference if needed.
                                if ( i == 0 && useRefAndAltAlleles )
                                    annotationValues.add(expressionValuesList.get(i++));
                                // use annotation expression and add to vc
                                annotationValues.add(expressionValuesList.get(i));
                                isAlleleConcordant = true;
                                canAnnotate = true;
                                break;
                            }
                            i++;
                        }

                        // can not find allele match so set to annotation value to zero
                        if ( !isAlleleConcordant )
                            annotationValues.add("0");
                    }

                    // no allele matches so can not annotate
                    if ( !canAnnotate )
                        continue;

                    // add the annotation values
                    infoAnnotations.put(expression.fullName, annotationValues);
                } else {
                    // use all of the expression values
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

    /**
     * Break the variant context into bialleles (reference and alternate alleles) and trim to a minimum representation
     *
     * @param vc variant context to annotate
     * @return list of biallelics trimmed to a minimum representation
     */
    private List<VariantContext> getMinRepresentationBiallelics(final VariantContext vc) {
        final List<VariantContext> minRepresentationBiallelicVCs = new ArrayList<VariantContext>();
        final boolean isMultiAllelic = vc.getNAlleles() > 2;
        if (isMultiAllelic) {
            final List<VariantContext> vcList = GATKVariantContextUtils.splitVariantContextToBiallelics(vc);
            for (final VariantContext biallelicVC : vcList) {
                if (!biallelicVC.isSNP())
                    minRepresentationBiallelicVCs.add(GATKVariantContextUtils.trimAlleles(biallelicVC, true, true));
                else
                    minRepresentationBiallelicVCs.add(biallelicVC);
            }
        } else {
            minRepresentationBiallelicVCs.add(vc);
        }

        return minRepresentationBiallelicVCs;
    }

    private void setReducibleAnnotations() {
        for(final InfoFieldAnnotation annotationType : requestedInfoAnnotations) {
            if (annotationType instanceof ReducibleAnnotation)
                requestedReducibleInfoAnnotations.add(annotationType);
            else
                requestedNonReducibleInfoAnnotations.add(annotationType);
        }
    }
}
