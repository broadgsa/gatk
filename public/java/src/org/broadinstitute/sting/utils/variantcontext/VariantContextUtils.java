/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlEngine;
import org.apache.log4j.Logger;
import org.broad.tribble.util.popgen.HardyWeinbergCalculation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.AbstractVCFCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.Serializable;
import java.util.*;

public class VariantContextUtils {
    private static Logger logger = Logger.getLogger(VariantContextUtils.class);
    public final static String MERGE_INTERSECTION = "Intersection";
    public final static String MERGE_FILTER_IN_ALL = "FilteredInAll";
    public final static String MERGE_REF_IN_ALL = "ReferenceInAll";
    public final static String MERGE_FILTER_PREFIX = "filterIn";

    final public static JexlEngine engine = new JexlEngine();
    static {
        engine.setSilent(false); // will throw errors now for selects that don't evaluate properly
        engine.setLenient(false);
    }

    /**
     * Update the attributes of the attributes map given the VariantContext to reflect the
     * proper chromosome-based VCF tags
     *
     * @param vc          the VariantContext
     * @param attributes  the attributes map to populate; must not be null; may contain old values
     * @param removeStaleValues should we remove stale values from the mapping?
     * @return the attributes map provided as input, returned for programming convenience
     */
    public static Map<String, Object> calculateChromosomeCounts(VariantContext vc, Map<String, Object> attributes, boolean removeStaleValues) {
        // if everyone is a no-call, remove the old attributes if requested
        if ( vc.getCalledChrCount() == 0 && removeStaleValues ) {
            if ( attributes.containsKey(VCFConstants.ALLELE_COUNT_KEY) )
                attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
            if ( attributes.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY) )
                attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
            if ( attributes.containsKey(VCFConstants.ALLELE_NUMBER_KEY) )
                attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
            return attributes;
        }

        if ( vc.hasGenotypes() ) {
            attributes.put(VCFConstants.ALLELE_NUMBER_KEY, vc.getCalledChrCount());

            // if there are alternate alleles, record the relevant tags
            if ( vc.getAlternateAlleles().size() > 0 ) {
                ArrayList<String> alleleFreqs = new ArrayList<String>();
                ArrayList<Integer> alleleCounts = new ArrayList<Integer>();
                double totalChromosomes = (double)vc.getCalledChrCount();
                for ( Allele allele : vc.getAlternateAlleles() ) {
                    int altChromosomes = vc.getCalledChrCount(allele);
                    alleleCounts.add(altChromosomes);
                    String freq = String.format(makePrecisionFormatStringFromDenominatorValue(totalChromosomes), ((double)altChromosomes / totalChromosomes));
                    alleleFreqs.add(freq);
                }

                attributes.put(VCFConstants.ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
                attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs.size() == 1 ? alleleFreqs.get(0) : alleleFreqs);
            }
            else {
                attributes.put(VCFConstants.ALLELE_COUNT_KEY, 0);
                attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0);
            }
        }

        return attributes;
    }

    /**
     * Update the attributes of the attributes map in the VariantContextBuilder to reflect the proper
     * chromosome-based VCF tags based on the current VC produced by builder.make()
     *
     * @param builder     the VariantContextBuilder we are updating
     * @param removeStaleValues should we remove stale values from the mapping?
     */
    public static void calculateChromosomeCounts(VariantContextBuilder builder, boolean removeStaleValues) {
        final VariantContext vc = builder.make();

        // if everyone is a no-call, remove the old attributes if requested
        if ( vc.getCalledChrCount() == 0 && removeStaleValues ) {
            if ( vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) )
                builder.rmAttribute(VCFConstants.ALLELE_COUNT_KEY);
            if ( vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) )
                builder.rmAttribute(VCFConstants.ALLELE_FREQUENCY_KEY);
            if ( vc.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) )
                builder.rmAttribute(VCFConstants.ALLELE_NUMBER_KEY);
            return;
        }

        if ( vc.hasGenotypes() ) {
            builder.attribute(VCFConstants.ALLELE_NUMBER_KEY, vc.getCalledChrCount());

            // if there are alternate alleles, record the relevant tags
            if ( vc.getAlternateAlleles().size() > 0 ) {
                ArrayList<String> alleleFreqs = new ArrayList<String>();
                ArrayList<Integer> alleleCounts = new ArrayList<Integer>();
                double totalChromosomes = (double)vc.getCalledChrCount();
                for ( Allele allele : vc.getAlternateAlleles() ) {
                    int altChromosomes = vc.getCalledChrCount(allele);
                    alleleCounts.add(altChromosomes);
                    String freq = String.format(makePrecisionFormatStringFromDenominatorValue(totalChromosomes), ((double)altChromosomes / totalChromosomes));
                    alleleFreqs.add(freq);
                }

                builder.attribute(VCFConstants.ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
                builder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs.size() == 1 ? alleleFreqs.get(0) : alleleFreqs);
            }
            else {
                builder.attribute(VCFConstants.ALLELE_COUNT_KEY, 0);
                builder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0);
            }
        }
    }

    private static String makePrecisionFormatStringFromDenominatorValue(double maxValue) {
        int precision = 1;

        while ( maxValue > 1 ) {
            precision++;
            maxValue /= 10.0;
        }

        return "%." + precision + "f";
    }

    public static Genotype removePLs(Genotype g) {
        Map<String, Object> attrs = new HashMap<String, Object>(g.getAttributes());
        attrs.remove(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY);
        attrs.remove(VCFConstants.GENOTYPE_LIKELIHOODS_KEY);
        return new Genotype(g.getSampleName(), g.getAlleles(), g.getLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, attrs, g.isPhased());
    }

    public static VariantContext createVariantContextWithPaddedAlleles(VariantContext inputVC, boolean refBaseShouldBeAppliedToEndOfAlleles) {
        // see if we need to pad common reference base from all alleles
        boolean padVC;

        // We need to pad a VC with a common base if the length of the reference allele is less than the length of the VariantContext.
        // This happens because the position of e.g. an indel is always one before the actual event (as per VCF convention).
        long locLength = (inputVC.getEnd() - inputVC.getStart()) + 1;
        if (inputVC.hasSymbolicAlleles())
            padVC = true;
        else if (inputVC.getReference().length() == locLength)
            padVC = false;
        else if (inputVC.getReference().length() == locLength-1)
            padVC = true;
        else throw new IllegalArgumentException("Badly formed variant context at location " + String.valueOf(inputVC.getStart()) +
                    " in contig " + inputVC.getChr() + ". Reference length must be at most one base shorter than location size");

        // nothing to do if we don't need to pad bases
        if (padVC) {
            if ( !inputVC.hasReferenceBaseForIndel() )
                throw new ReviewedStingException("Badly formed variant context at location " + inputVC.getChr() + ":" + inputVC.getStart() + "; no padded reference base is available.");

            Byte refByte = inputVC.getReferenceBaseForIndel();

            List<Allele> alleles = new ArrayList<Allele>();

            for (Allele a : inputVC.getAlleles()) {
                // get bases for current allele and create a new one with trimmed bases
                if (a.isSymbolic()) {
                    alleles.add(a);
                } else {
                    String newBases;
                    if ( refBaseShouldBeAppliedToEndOfAlleles )
                        newBases = a.getBaseString() + new String(new byte[]{refByte});
                    else
                        newBases = new String(new byte[]{refByte}) + a.getBaseString();
                    alleles.add(Allele.create(newBases,a.isReference()));
                }
            }

            // now we can recreate new genotypes with trimmed alleles
            GenotypesContext genotypes = GenotypesContext.create(inputVC.getNSamples());
            for (final Genotype g : inputVC.getGenotypes() ) {
                List<Allele> inAlleles = g.getAlleles();
                List<Allele> newGenotypeAlleles = new ArrayList<Allele>(g.getAlleles().size());
                for (Allele a : inAlleles) {
                    if (a.isCalled()) {
                        if (a.isSymbolic()) {
                            newGenotypeAlleles.add(a);
                        } else {
                            String newBases;
                            if ( refBaseShouldBeAppliedToEndOfAlleles )
                                newBases = a.getBaseString() + new String(new byte[]{refByte});
                            else
                                newBases = new String(new byte[]{refByte}) + a.getBaseString();
                            newGenotypeAlleles.add(Allele.create(newBases,a.isReference()));
                        }
                    }
                    else {
                        // add no-call allele
                        newGenotypeAlleles.add(Allele.NO_CALL);
                    }
                }
                genotypes.add(new Genotype(g.getSampleName(), newGenotypeAlleles, g.getLog10PError(),
                        g.getFilters(), g.getAttributes(), g.isPhased()));

            }

            return new VariantContextBuilder(inputVC).alleles(alleles).genotypes(genotypes).make();
        }
        else
            return inputVC;

    }

    /**
     * A simple but common wrapper for matching VariantContext objects using JEXL expressions
     */
    public static class JexlVCMatchExp {
        public String name;
        public Expression exp;

        /**
         * Create a new matcher expression with name and JEXL expression exp
         * @param name name
         * @param exp  expression
         */
        public JexlVCMatchExp(String name, Expression exp) {
            this.name = name;
            this.exp = exp;
        }
    }

    /**
     * Method for creating JexlVCMatchExp from input walker arguments names and exps.  These two arrays contain
     * the name associated with each JEXL expression. initializeMatchExps will parse each expression and return
     * a list of JexlVCMatchExp, in order, that correspond to the names and exps.  These are suitable input to
     * match() below.
     *
     * @param names names
     * @param exps  expressions
     * @return list of matches
     */
    public static List<JexlVCMatchExp> initializeMatchExps(String[] names, String[] exps) {
        if ( names == null || exps == null )
            throw new ReviewedStingException("BUG: neither names nor exps can be null: names " + Arrays.toString(names) + " exps=" + Arrays.toString(exps) );

        if ( names.length != exps.length )
            throw new UserException("Inconsistent number of provided filter names and expressions: names=" + Arrays.toString(names) + " exps=" + Arrays.toString(exps));

        Map<String, String> map = new HashMap<String, String>();
        for ( int i = 0; i < names.length; i++ ) { map.put(names[i], exps[i]); }

        return VariantContextUtils.initializeMatchExps(map);
    }

    public static List<JexlVCMatchExp> initializeMatchExps(ArrayList<String> names, ArrayList<String> exps) {
        String[] nameArray = new String[names.size()];
        String[] expArray = new String[exps.size()];
        return initializeMatchExps(names.toArray(nameArray), exps.toArray(expArray));
    }


    /**
     * Method for creating JexlVCMatchExp from input walker arguments mapping from names to exps.  These two arrays contain
     * the name associated with each JEXL expression. initializeMatchExps will parse each expression and return
     * a list of JexlVCMatchExp, in order, that correspond to the names and exps.  These are suitable input to
     * match() below.
     *
     * @param names_and_exps mapping of names to expressions
     * @return list of matches
     */
    public static List<JexlVCMatchExp> initializeMatchExps(Map<String, String> names_and_exps) {
        List<JexlVCMatchExp> exps = new ArrayList<JexlVCMatchExp>();

        for ( Map.Entry<String, String> elt : names_and_exps.entrySet() ) {
            String name = elt.getKey();
            String expStr = elt.getValue();

            if ( name == null || expStr == null ) throw new IllegalArgumentException("Cannot create null expressions : " + name +  " " + expStr);
            try {
                Expression exp = engine.createExpression(expStr);
                exps.add(new JexlVCMatchExp(name, exp));
            } catch (Exception e) {
                throw new UserException.BadArgumentValue(name, "Invalid expression used (" + expStr + "). Please see the JEXL docs for correct syntax.") ;
            }
        }

        return exps;
    }

    /**
     * Returns true if exp match VC.  See collection<> version for full docs.
     * @param vc    variant context
     * @param exp   expression
     * @return true if there is a match
     */
    public static boolean match(VariantContext vc, JexlVCMatchExp exp) {
        return match(vc,Arrays.asList(exp)).get(exp);
    }

    /**
     * Matches each JexlVCMatchExp exp against the data contained in vc, and returns a map from these
     * expressions to true (if they matched) or false (if they didn't).  This the best way to apply JEXL
     * expressions to VariantContext records.  Use initializeMatchExps() to create the list of JexlVCMatchExp
     * expressions.
     *
     * @param vc   variant context
     * @param exps expressions
     * @return true if there is a match
     */
    public static Map<JexlVCMatchExp, Boolean> match(VariantContext vc, Collection<JexlVCMatchExp> exps) {
        return new JEXLMap(exps,vc);

    }

    /**
     * Returns true if exp match VC/g.  See collection<> version for full docs.
     * @param vc   variant context
     * @param g    genotype
     * @param exp   expression
     * @return true if there is a match
     */
    public static boolean match(VariantContext vc, Genotype g, JexlVCMatchExp exp) {
        return match(vc,g,Arrays.asList(exp)).get(exp);
    }

    /**
     * Matches each JexlVCMatchExp exp against the data contained in vc/g, and returns a map from these
     * expressions to true (if they matched) or false (if they didn't).  This the best way to apply JEXL
     * expressions to VariantContext records/genotypes.  Use initializeMatchExps() to create the list of JexlVCMatchExp
     * expressions.
     *
     * @param vc   variant context
     * @param g    genotype
     * @param exps expressions
     * @return true if there is a match
     */
    public static Map<JexlVCMatchExp, Boolean> match(VariantContext vc, Genotype g, Collection<JexlVCMatchExp> exps) {
        return new JEXLMap(exps,vc,g);
    }

    public static double computeHardyWeinbergPvalue(VariantContext vc) {
        if ( vc.getCalledChrCount() == 0 )
            return 0.0;
        return HardyWeinbergCalculation.hwCalculate(vc.getHomRefCount(), vc.getHetCount(), vc.getHomVarCount());
    }

    /**
     * Returns a newly allocated VC that is the same as VC, but without genotypes
     * @param vc  variant context
     * @return  new VC without genotypes
     */
    @Requires("vc != null")
    @Ensures("result != null")
    public static VariantContext sitesOnlyVariantContext(VariantContext vc) {
        return new VariantContextBuilder(vc).noGenotypes().make();
    }

    /**
     * Returns a newly allocated list of VC, where each VC is the same as the input VCs, but without genotypes
     * @param vcs  collection of VCs
     * @return new VCs without genotypes
     */
    @Requires("vcs != null")
    @Ensures("result != null")
    public static Collection<VariantContext> sitesOnlyVariantContexts(Collection<VariantContext> vcs) {
        List<VariantContext> r = new ArrayList<VariantContext>();
        for ( VariantContext vc : vcs )
            r.add(sitesOnlyVariantContext(vc));
        return r;
    }

    private final static Map<String, Object> subsetAttributes(final CommonInfo igc, final Collection<String> keysToPreserve) {
        Map<String, Object> attributes = new HashMap<String, Object>(keysToPreserve.size());
        for ( final String key : keysToPreserve  ) {
            if ( igc.hasAttribute(key) )
                attributes.put(key, igc.getAttribute(key));
        }
        return attributes;
    }

    /**
     * @deprecated use variant context builder version instead
     * @param vc
     * @param keysToPreserve
     * @return
     */
    @Deprecated
    public static VariantContext pruneVariantContext(final VariantContext vc, Collection<String> keysToPreserve ) {
        return pruneVariantContext(new VariantContextBuilder(vc), keysToPreserve).make();
    }

    public static VariantContextBuilder pruneVariantContext(final VariantContextBuilder builder, Collection<String> keysToPreserve ) {
        final VariantContext vc = builder.make();
        if ( keysToPreserve == null ) keysToPreserve = Collections.emptyList();

        // VC info
        final Map<String, Object> attributes = subsetAttributes(vc.commonInfo, keysToPreserve);

        // Genotypes
        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype g : vc.getGenotypes() ) {
            Map<String, Object> genotypeAttributes = subsetAttributes(g.commonInfo, keysToPreserve);
            genotypes.add(new Genotype(g.getSampleName(), g.getAlleles(), g.getLog10PError(), g.getFilters(),
                    genotypeAttributes, g.isPhased()));
        }

        return builder.genotypes(genotypes).attributes(attributes);
    }

    public enum GenotypeMergeType {
        /**
         * Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
         */
        UNIQUIFY,
        /**
         * Take genotypes in priority order (see the priority argument).
         */
        PRIORITIZE,
        /**
         * Take the genotypes in any order.
         */
        UNSORTED,
        /**
         * Require that all samples/genotypes be unique between all inputs.
         */
        REQUIRE_UNIQUE
    }

    public enum FilteredRecordMergeType {
        /**
         * Union - leaves the record if any record is unfiltered.
         */
        KEEP_IF_ANY_UNFILTERED,
        /**
         * Requires all records present at site to be unfiltered. VCF files that don't contain the record don't influence this.
         */
        KEEP_IF_ALL_UNFILTERED
    }

    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniqifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name
     *
     * @param genomeLocParser           loc parser
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param annotateOrigin            should we annotate the set it came from?
     * @param printMessages             should we print messages?
     * @param setKey                    the key name of the set
     * @param filteredAreUncalled       are filtered records uncalled?
     * @param mergeInfoWithMaxAC        should we merge in info from the VC with maximum allele count?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final GenomeLocParser genomeLocParser,
                                             final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean annotateOrigin,
                                             final boolean printMessages,
                                             final String setKey,
                                             final boolean filteredAreUncalled,
                                             final boolean mergeInfoWithMaxAC ) {
        if ( unsortedVCs == null || unsortedVCs.size() == 0 )
            return null;

        if ( annotateOrigin && priorityListOfVCs == null )
            throw new IllegalArgumentException("Cannot merge calls and annotate their origins without a complete priority list of VariantContexts");

        if ( genotypeMergeOptions == GenotypeMergeType.REQUIRE_UNIQUE )
            verifyUniqueSampleNames(unsortedVCs);

        List<VariantContext> prepaddedVCs = sortVariantContextsByPriority(unsortedVCs, priorityListOfVCs, genotypeMergeOptions);
        // Make sure all variant contexts are padded with reference base in case of indels if necessary
        List<VariantContext> VCs = new ArrayList<VariantContext>();

        for (VariantContext vc : prepaddedVCs) {
            // also a reasonable place to remove filtered calls, if needed
            if ( ! filteredAreUncalled || vc.isNotFiltered() )
                VCs.add(createVariantContextWithPaddedAlleles(vc, false));
        }
        if ( VCs.size() == 0 ) // everything is filtered out and we're filteredAreUncalled
            return null;

        // establish the baseline info from the first VC
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();
        final Allele refAllele = determineReferenceAllele(VCs);

        final Set<Allele> alleles = new LinkedHashSet<Allele>();
        final Set<String> filters = new TreeSet<String>();
        final Map<String, Object> attributes = new TreeMap<String, Object>();
        final Set<String> inconsistentAttributes = new HashSet<String>();
        final Set<String> variantSources = new HashSet<String>(); // contains the set of sources we found in our set of VCs that are variant
        final Set<String> rsIDs = new LinkedHashSet<String>(1); // most of the time there's one id

        GenomeLoc loc = getLocation(genomeLocParser,first);
        int depth = 0;
        int maxAC = -1;
        final Map<String, Object> attributesWithMaxAC = new TreeMap<String, Object>();
        double log10PError = 1;
        VariantContext vcWithMaxAC = null;
        GenotypesContext genotypes = GenotypesContext.create();

        // counting the number of filtered and variant VCs
        int nFiltered = 0;

        boolean remapped = false;

        // cycle through and add info from the other VCs, making sure the loc/reference matches

        for ( VariantContext vc : VCs ) {
            if ( loc.getStart() != vc.getStart() ) // || !first.getReference().equals(vc.getReference()) )
                throw new ReviewedStingException("BUG: attempting to merge VariantContexts with different start sites: first="+ first.toString() + " second=" + vc.toString());

            if ( getLocation(genomeLocParser,vc).size() > loc.size() )
                loc = getLocation(genomeLocParser,vc); // get the longest location

            nFiltered += vc.isFiltered() ? 1 : 0;
            if ( vc.isVariant() ) variantSources.add(vc.getSource());

            AlleleMapper alleleMapping = resolveIncompatibleAlleles(refAllele, vc, alleles);
            remapped = remapped || alleleMapping.needsRemapping();

            alleles.addAll(alleleMapping.values());

            mergeGenotypes(genotypes, vc, alleleMapping, genotypeMergeOptions == GenotypeMergeType.UNIQUIFY);

            log10PError = Math.min(log10PError, vc.isVariant() ? vc.getLog10PError() : 1);

            filters.addAll(vc.getFilters());

            //
            // add attributes
            //
            // special case DP (add it up) and ID (just preserve it)
            //
            if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            if ( vc.hasID() ) rsIDs.add(vc.getID());
            if (mergeInfoWithMaxAC && vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
                String rawAlleleCounts = vc.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY, null);
                // lets see if the string contains a , separator
                if (rawAlleleCounts.contains(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) {
                    List<String> alleleCountArray = Arrays.asList(rawAlleleCounts.substring(1, rawAlleleCounts.length() - 1).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
                    for (String alleleCount : alleleCountArray) {
                        final int ac = Integer.valueOf(alleleCount.trim());
                        if (ac > maxAC) {
                            maxAC = ac;
                            vcWithMaxAC = vc;
                        }
                    }
                } else {
                    final int ac = Integer.valueOf(rawAlleleCounts);
                    if (ac > maxAC) {
                        maxAC = ac;
                        vcWithMaxAC = vc;
                    }
                }
            }

            for (Map.Entry<String, Object> p : vc.getAttributes().entrySet()) {
                String key = p.getKey();
                // if we don't like the key already, don't go anywhere
                if ( ! inconsistentAttributes.contains(key) ) {
                    boolean alreadyFound = attributes.containsKey(key);
                    Object boundValue = attributes.get(key);
                    boolean boundIsMissingValue = alreadyFound && boundValue.equals(VCFConstants.MISSING_VALUE_v4);

                    if ( alreadyFound && ! boundValue.equals(p.getValue()) && ! boundIsMissingValue ) {
                        // we found the value but we're inconsistent, put it in the exclude list
                        //System.out.printf("Inconsistent INFO values: %s => %s and %s%n", key, boundValue, p.getValue());
                        inconsistentAttributes.add(key);
                        attributes.remove(key);
                    } else if ( ! alreadyFound || boundIsMissingValue )  { // no value
                        //if ( vc != first ) System.out.printf("Adding key %s => %s%n", p.getKey(), p.getValue());
                        attributes.put(key, p.getValue());
                    }
                }
            }
        }

        // if we have more alternate alleles in the merged VC than in one or more of the
        // original VCs, we need to strip out the GL/PLs (because they are no longer accurate), as well as allele-dependent attributes like AC,AF
        for ( VariantContext vc : VCs ) {
            if (vc.alleles.size() == 1)
                continue;
            if ( hasPLIncompatibleAlleles(alleles, vc.alleles)) {
                logger.warn(String.format("Stripping PLs at %s due incompatible alleles merged=%s vs. single=%s",
                        genomeLocParser.createGenomeLoc(vc), alleles, vc.alleles));
                genotypes = stripPLs(genotypes);
                // this will remove stale AC,AF attributed from vc
                calculateChromosomeCounts(vc, attributes, true);
                break;
            }
        }

        // take the VC with the maxAC and pull the attributes into a modifiable map
        if ( mergeInfoWithMaxAC && vcWithMaxAC != null ) {
            attributesWithMaxAC.putAll(vcWithMaxAC.getAttributes());
        }

        // if at least one record was unfiltered and we want a union, clear all of the filters
        if ( filteredRecordMergeType == FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED && nFiltered != VCs.size() )
            filters.clear();


        if ( annotateOrigin ) { // we care about where the call came from
            String setValue;
            if ( nFiltered == 0 && variantSources.size() == priorityListOfVCs.size() ) // nothing was unfiltered
                setValue = MERGE_INTERSECTION;
            else if ( nFiltered == VCs.size() )     // everything was filtered out
                setValue = MERGE_FILTER_IN_ALL;
            else if ( variantSources.isEmpty() )               // everyone was reference
                setValue = MERGE_REF_IN_ALL;
            else {
                LinkedHashSet<String> s = new LinkedHashSet<String>();
                for ( VariantContext vc : VCs )
                    if ( vc.isVariant() )
                        s.add( vc.isFiltered() ? MERGE_FILTER_PREFIX + vc.getSource() : vc.getSource() );
                setValue = Utils.join("-", s);
            }

            if ( setKey != null ) {
                attributes.put(setKey, setValue);
                if( mergeInfoWithMaxAC && vcWithMaxAC != null ) { attributesWithMaxAC.put(setKey, vcWithMaxAC.getSource()); }
            }
        }

        if ( depth > 0 )
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID);
        builder.loc(loc.getContig(), loc.getStart(), loc.getStop());
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        builder.log10PError(log10PError);
        builder.filters(filters).attributes(mergeInfoWithMaxAC ? attributesWithMaxAC : attributes);

        // Trim the padded bases of all alleles if necessary
        VariantContext merged = createVariantContextWithTrimmedAlleles(builder.make());
        if ( printMessages && remapped ) System.out.printf("Remapped => %s%n", merged);
        return merged;
    }

    private static final boolean hasPLIncompatibleAlleles(final Collection<Allele> alleleSet1, final Collection<Allele> alleleSet2) {
        final Iterator<Allele> it1 = alleleSet1.iterator();
        final Iterator<Allele> it2 = alleleSet2.iterator();

        while ( it1.hasNext() && it2.hasNext() ) {
            final Allele a1 = it1.next();
            final Allele a2 = it2.next();
            if ( ! a1.equals(a2) )
                return true;
        }

        // by this point, at least one of the iterators is empty.  All of the elements
        // we've compared are equal up until this point.  But it's possible that the
        // sets aren't the same size, which is indicated by the test below.  If they
        // are of the same size, though, the sets are compatible
        return it1.hasNext() || it2.hasNext();
    }

    public static boolean allelesAreSubset(VariantContext vc1, VariantContext vc2) {
        // if all alleles of vc1 are a contained in alleles of vc2, return true
        if (!vc1.getReference().equals(vc2.getReference()))
            return false;

        for (Allele a :vc1.getAlternateAlleles()) {
            if (!vc2.getAlternateAlleles().contains(a))
                return false;
        }

        return true;
    }

    public static VariantContext createVariantContextWithTrimmedAlleles(VariantContext inputVC) {
        // see if we need to trim common reference base from all alleles
        boolean trimVC;

        // We need to trim common reference base from all alleles in all genotypes if a ref base is common to all alleles
        Allele refAllele = inputVC.getReference();
        if (!inputVC.isVariant())
            trimVC = false;
        else if (refAllele.isNull())
            trimVC = false;
        else {
            trimVC = (AbstractVCFCodec.computeForwardClipping(new ArrayList<Allele>(inputVC.getAlternateAlleles()),
                    inputVC.getReference().getDisplayString()) > 0);
         }

        // nothing to do if we don't need to trim bases
        if (trimVC) {
            List<Allele> alleles = new ArrayList<Allele>();
            GenotypesContext genotypes = GenotypesContext.create();

            // set the reference base for indels in the attributes
            Map<String,Object> attributes = new TreeMap<String,Object>(inputVC.getAttributes());

            Map<Allele, Allele> originalToTrimmedAlleleMap = new HashMap<Allele, Allele>();

            for (Allele a : inputVC.getAlleles()) {
                if (a.isSymbolic()) {
                    alleles.add(a);
                    originalToTrimmedAlleleMap.put(a, a);
                } else {
                    // get bases for current allele and create a new one with trimmed bases
                    byte[] newBases = Arrays.copyOfRange(a.getBases(), 1, a.length());
                    Allele trimmedAllele = Allele.create(newBases, a.isReference());
                    alleles.add(trimmedAllele);
                    originalToTrimmedAlleleMap.put(a, trimmedAllele);
                }
            }

            // detect case where we're trimming bases but resulting vc doesn't have any null allele. In that case, we keep original representation
            // example: mixed records such as {TA*,TGA,TG}
            boolean hasNullAlleles = false;

            for (Allele a: originalToTrimmedAlleleMap.values()) {
                if (a.isNull())
                    hasNullAlleles = true;
                if (a.isReference())
                    refAllele = a;
             }

             if (!hasNullAlleles)
               return inputVC;
           // now we can recreate new genotypes with trimmed alleles
            for ( final Genotype genotype : inputVC.getGenotypes() ) {

                List<Allele> originalAlleles = genotype.getAlleles();
                List<Allele> trimmedAlleles = new ArrayList<Allele>();
                for ( Allele a : originalAlleles ) {
                    if ( a.isCalled() )
                        trimmedAlleles.add(originalToTrimmedAlleleMap.get(a));
                    else
                        trimmedAlleles.add(Allele.NO_CALL);
                }
                genotypes.add(Genotype.modifyAlleles(genotype, trimmedAlleles));

            }

            final VariantContextBuilder builder = new VariantContextBuilder(inputVC);
            return builder.alleles(alleles).genotypes(genotypes).attributes(attributes).referenceBaseForIndel(new Byte(inputVC.getReference().getBases()[0])).make();
        }

        return inputVC;
    }

    public static GenotypesContext stripPLs(GenotypesContext genotypes) {
        GenotypesContext newGs = GenotypesContext.create(genotypes.size());

        for ( final Genotype g : genotypes ) {
            newGs.add(g.hasLikelihoods() ? removePLs(g) : g);
        }

        return newGs;
    }

    public static Map<VariantContext.Type, List<VariantContext>> separateVariantContextsByType(Collection<VariantContext> VCs) {
        HashMap<VariantContext.Type, List<VariantContext>> mappedVCs = new HashMap<VariantContext.Type, List<VariantContext>>();
        for ( VariantContext vc : VCs ) {

            // look at previous variant contexts of different type. If:
            // a) otherVC has alleles which are subset of vc, remove otherVC from its list and add otherVC to  vc's list
            // b) vc has alleles which are subset of otherVC. Then, add vc to otherVC's type list (rather, do nothing since vc will be added automatically to its list)
            // c) neither: do nothing, just add vc to its own list
            boolean addtoOwnList = true;
            for (VariantContext.Type type : VariantContext.Type.values()) {
                if (type.equals(vc.getType()))
                    continue;

                if (!mappedVCs.containsKey(type))
                    continue;

                List<VariantContext> vcList = mappedVCs.get(type);
                for (int k=0; k <  vcList.size(); k++) {
                    VariantContext otherVC = vcList.get(k);
                    if (allelesAreSubset(otherVC,vc)) {
                        // otherVC has a type different than vc and its alleles are a subset of vc: remove otherVC from its list and add it to vc's type list
                        vcList.remove(k);
                        // avoid having empty lists
                        if (vcList.size() == 0)
                            mappedVCs.remove(vcList);
                        if ( !mappedVCs.containsKey(vc.getType()) )
                            mappedVCs.put(vc.getType(), new ArrayList<VariantContext>());
                        mappedVCs.get(vc.getType()).add(otherVC);
                        break;
                    }
                    else if (allelesAreSubset(vc,otherVC)) {
                        // vc has a type different than otherVC and its alleles are a subset of VC: add vc to otherVC's type list and don't add to its own
                        mappedVCs.get(type).add(vc);
                        addtoOwnList = false;
                        break;
                    }
                }
            }
            if (addtoOwnList) {
                if ( !mappedVCs.containsKey(vc.getType()) )
                    mappedVCs.put(vc.getType(), new ArrayList<VariantContext>());
                mappedVCs.get(vc.getType()).add(vc);
                }
        }

        return mappedVCs;
    }

    private static class AlleleMapper {
        private VariantContext vc = null;
        private Map<Allele, Allele> map = null;
        public AlleleMapper(VariantContext vc)          { this.vc = vc; }
        public AlleleMapper(Map<Allele, Allele> map)    { this.map = map; }
        public boolean needsRemapping()                 { return this.map != null; }
        public Collection<Allele> values()              { return map != null ? map.values() : vc.getAlleles(); }

        public Allele remap(Allele a)                   { return map != null && map.containsKey(a) ? map.get(a) : a; }

        public List<Allele> remap(List<Allele> as) {
            List<Allele> newAs = new ArrayList<Allele>();
            for ( Allele a : as ) {
                //System.out.printf("  Remapping %s => %s%n", a, remap(a));
                newAs.add(remap(a));
            }
            return newAs;
        }
    }

    static private void verifyUniqueSampleNames(Collection<VariantContext> unsortedVCs) {
        Set<String> names = new HashSet<String>();
        for ( VariantContext vc : unsortedVCs ) {
            for ( String name : vc.getSampleNames() ) {
                //System.out.printf("Checking %s %b%n", name, names.contains(name));
                if ( names.contains(name) )
                    throw new UserException("REQUIRE_UNIQUE sample names is true but duplicate names were discovered " + name);
            }

            names.addAll(vc.getSampleNames());
        }
    }


    static private Allele determineReferenceAllele(List<VariantContext> VCs) {
        Allele ref = null;

        for ( VariantContext vc : VCs ) {
            Allele myRef = vc.getReference();
            if ( ref == null || ref.length() < myRef.length() )
                ref = myRef;
            else if ( ref.length() == myRef.length() && ! ref.equals(myRef) )
                throw new UserException.BadInput(String.format("The provided variant file(s) have inconsistent references for the same position(s) at %s:%d, %s vs. %s", vc.getChr(), vc.getStart(), ref, myRef));
        }

        return ref;
    }

    static private AlleleMapper resolveIncompatibleAlleles(Allele refAllele, VariantContext vc, Set<Allele> allAlleles) {
        if ( refAllele.equals(vc.getReference()) )
            return new AlleleMapper(vc);
        else {
            // we really need to do some work.  The refAllele is the longest reference allele seen at this
            // start site.  So imagine it is:
            //
            // refAllele: ACGTGA
            // myRef:     ACGT
            // myAlt:     -
            //
            // We need to remap all of the alleles in vc to include the extra GA so that
            // myRef => refAllele and myAlt => GA
            //

            Allele myRef = vc.getReference();
            if ( refAllele.length() <= myRef.length() ) throw new ReviewedStingException("BUG: myRef="+myRef+" is longer than refAllele="+refAllele);
            byte[] extraBases = Arrays.copyOfRange(refAllele.getBases(), myRef.length(), refAllele.length());

//            System.out.printf("Remapping allele at %s%n", vc);
//            System.out.printf("ref   %s%n", refAllele);
//            System.out.printf("myref %s%n", myRef );
//            System.out.printf("extrabases %s%n", new String(extraBases));

            Map<Allele, Allele> map = new HashMap<Allele, Allele>();
            for ( Allele a : vc.getAlleles() ) {
                if ( a.isReference() )
                    map.put(a, refAllele);
                else {
                    Allele extended = Allele.extend(a, extraBases);
                    for ( Allele b : allAlleles )
                        if ( extended.equals(b) )
                            extended = b;
//                    System.out.printf("  Extending %s => %s%n", a, extended);
                    map.put(a, extended);
                }
            }

            // debugging
//            System.out.printf("mapping %s%n", map);

            return new AlleleMapper(map);
        }
    }

    static class CompareByPriority implements Comparator<VariantContext>, Serializable {
        List<String> priorityListOfVCs;
        public CompareByPriority(List<String> priorityListOfVCs) {
            this.priorityListOfVCs = priorityListOfVCs;
        }

        private int getIndex(VariantContext vc) {
            int i = priorityListOfVCs.indexOf(vc.getSource());
            if ( i == -1 ) throw new UserException.BadArgumentValue(Utils.join(",", priorityListOfVCs), "Priority list " + priorityListOfVCs + " doesn't contain variant context " + vc.getSource());
            return i;
        }

        public int compare(VariantContext vc1, VariantContext vc2) {
            return Integer.valueOf(getIndex(vc1)).compareTo(getIndex(vc2));
        }
    }

    public static List<VariantContext> sortVariantContextsByPriority(Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs, GenotypeMergeType mergeOption ) {
        if ( mergeOption == GenotypeMergeType.PRIORITIZE && priorityListOfVCs == null )
            throw new IllegalArgumentException("Cannot merge calls by priority with a null priority list");

        if ( priorityListOfVCs == null || mergeOption == GenotypeMergeType.UNSORTED )
            return new ArrayList<VariantContext>(unsortedVCs);
        else {
            ArrayList<VariantContext> sorted = new ArrayList<VariantContext>(unsortedVCs);
            Collections.sort(sorted, new CompareByPriority(priorityListOfVCs));
            return sorted;
        }
    }

    private static void mergeGenotypes(GenotypesContext mergedGenotypes, VariantContext oneVC, AlleleMapper alleleMapping, boolean uniqifySamples) {
        for ( Genotype g : oneVC.getGenotypes() ) {
            String name = mergedSampleName(oneVC.getSource(), g.getSampleName(), uniqifySamples);
            if ( mergedGenotypes.containsSample(name) ) {
                // only add if the name is new
                Genotype newG = g;

                if ( uniqifySamples || alleleMapping.needsRemapping() ) {
                    final List<Allele> alleles = alleleMapping.needsRemapping() ? alleleMapping.remap(g.getAlleles()) : g.getAlleles();
                    newG = new Genotype(name, alleles, g.getLog10PError(), g.getFilters(), g.getAttributes(), g.isPhased());
                }

                mergedGenotypes.add(newG);
            }
        }
    }

    public static String mergedSampleName(String trackName, String sampleName, boolean uniqify ) {
        return uniqify ? sampleName + "." + trackName : sampleName;
    }

    /**
     * Returns a context identical to this with the REF and ALT alleles reverse complemented.
     *
     * @param vc        variant context
     * @return new vc
     */
    public static VariantContext reverseComplement(VariantContext vc) {
        // create a mapping from original allele to reverse complemented allele
        HashMap<Allele, Allele> alleleMap = new HashMap<Allele, Allele>(vc.getAlleles().size());
        for ( Allele originalAllele : vc.getAlleles() ) {
            Allele newAllele;
            if ( originalAllele.isNoCall() || originalAllele.isNull() )
                newAllele = originalAllele;
            else
                newAllele = Allele.create(BaseUtils.simpleReverseComplement(originalAllele.getBases()), originalAllele.isReference());
            alleleMap.put(originalAllele, newAllele);
        }

        // create new Genotype objects
        GenotypesContext newGenotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            List<Allele> newAlleles = new ArrayList<Allele>();
            for ( Allele allele : genotype.getAlleles() ) {
                Allele newAllele = alleleMap.get(allele);
                if ( newAllele == null )
                    newAllele = Allele.NO_CALL;
                newAlleles.add(newAllele);
            }
            newGenotypes.add(Genotype.modifyAlleles(genotype, newAlleles));
        }

        return new VariantContextBuilder(vc).alleles(alleleMap.values()).genotypes(newGenotypes).make();
    }

    public static VariantContext purgeUnallowedGenotypeAttributes(VariantContext vc, Set<String> allowedAttributes) {
        if ( allowedAttributes == null )
            return vc;

        GenotypesContext newGenotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            Map<String, Object> attrs = new HashMap<String, Object>();
            for ( Map.Entry<String, Object> attr : genotype.getAttributes().entrySet() ) {
                if ( allowedAttributes.contains(attr.getKey()) )
                    attrs.put(attr.getKey(), attr.getValue());
            }
            newGenotypes.add(Genotype.modifyAttributes(genotype, attrs));
        }

        return new VariantContextBuilder(vc).genotypes(newGenotypes).make();
    }

    public static BaseUtils.BaseSubstitutionType getSNPSubstitutionType(VariantContext context) {
        if (!context.isSNP() || !context.isBiallelic())
            throw new IllegalStateException("Requested SNP substitution type for bialleic non-SNP " + context);
        return BaseUtils.SNPSubstitutionType(context.getReference().getBases()[0], context.getAlternateAllele(0).getBases()[0]);
    }

    /**
     * If this is a BiAlleic SNP, is it a transition?
     */
    public static boolean isTransition(VariantContext context) {
        return getSNPSubstitutionType(context) == BaseUtils.BaseSubstitutionType.TRANSITION;
    }

    /**
     * If this is a BiAlleic SNP, is it a transversion?
     */
    public static boolean isTransversion(VariantContext context) {
        return getSNPSubstitutionType(context) == BaseUtils.BaseSubstitutionType.TRANSVERSION;
    }

    /**
     * create a genome location, given a variant context
     * @param genomeLocParser parser
     * @param vc the variant context
     * @return the genomeLoc
     */
    public static final GenomeLoc getLocation(GenomeLocParser genomeLocParser,VariantContext vc) {
        return genomeLocParser.createGenomeLoc(vc.getChr(), vc.getStart(), vc.getEnd(), true);
    }

    public static final Set<String> genotypeNames(final Collection<Genotype> genotypes) {
        final Set<String> names = new HashSet<String>(genotypes.size());
        for ( final Genotype g : genotypes )
            names.add(g.getSampleName());
        return names;
    }
}
