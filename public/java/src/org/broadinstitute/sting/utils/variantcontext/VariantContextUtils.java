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
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.StringUtil;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlEngine;
import org.apache.log4j.Logger;
import org.broad.tribble.util.popgen.HardyWeinbergCalculation;
import org.broadinstitute.sting.gatk.walkers.phasing.ReadBackedPhasingWalker;
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
     * Create a new VariantContext
     *
     * @param name            name
     * @param loc             location
     * @param alleles         alleles
     * @param genotypes       genotypes set
     * @param negLog10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     * @return VariantContext object
     */
    public static VariantContext toVC(String name, GenomeLoc loc, Collection<Allele> alleles, Collection<Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        return new VariantContext(name, loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes != null ? VariantContext.genotypeCollectionToMap(new TreeMap<String, Genotype>(), genotypes) : null, negLog10PError, filters, attributes);
    }

    /**
     * Create a new variant context without genotypes and no Perror, no filters, and no attributes
     * @param name            name
     * @param loc             location
     * @param alleles         alleles
     * @return VariantContext object
     */
    public static VariantContext toVC(String name, GenomeLoc loc, Collection<Allele> alleles) {
        return new VariantContext (name, loc.getContig(), loc.getStart(), loc.getStop(), alleles, VariantContext.NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    /**
     * Create a new variant context without genotypes and no Perror, no filters, and no attributes
     * @param name            name
     * @param loc             location
     * @param alleles         alleles
     * @param genotypes       genotypes
     * @return VariantContext object
     */
    public static VariantContext toVC(String name, GenomeLoc loc, Collection<Allele> alleles, Collection<Genotype> genotypes) {
        return new VariantContext(name, loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    /**
     * Copy constructor
     *
     * @param other the VariantContext to copy
     * @return VariantContext object
     */
    public static VariantContext toVC(VariantContext other) {
        return new VariantContext(other.getSource(), other.getChr(), other.getStart(), other.getEnd(), other.getAlleles(), other.getGenotypes(), other.getNegLog10PError(), other.getFilters(), other.getAttributes());
    }

    /**
     * Update the attributes of the attributes map given the VariantContext to reflect the proper chromosome-based VCF tags
     *
     * @param vc          the VariantContext
     * @param attributes  the attributes map to populate; must not be null; may contain old values
     * @param removeStaleValues should we remove stale values from the mapping?
     */
    public static void calculateChromosomeCounts(VariantContext vc, Map<String, Object> attributes, boolean removeStaleValues) {
        // if everyone is a no-call, remove the old attributes if requested
        if ( vc.getChromosomeCount() == 0 && removeStaleValues ) {
            if ( attributes.containsKey(VCFConstants.ALLELE_COUNT_KEY) )
                attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
            if ( attributes.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY) )
                attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
            if ( attributes.containsKey(VCFConstants.ALLELE_NUMBER_KEY) )
                attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
            return;
        }

        if ( vc.hasGenotypes() ) {
            attributes.put(VCFConstants.ALLELE_NUMBER_KEY, vc.getChromosomeCount());

            // if there are alternate alleles, record the relevant tags
            if ( vc.getAlternateAlleles().size() > 0 ) {
                ArrayList<String> alleleFreqs = new ArrayList<String>();
                ArrayList<Integer> alleleCounts = new ArrayList<Integer>();
                double totalChromosomes = (double)vc.getChromosomeCount();
                for ( Allele allele : vc.getAlternateAlleles() ) {
                    int altChromosomes = vc.getChromosomeCount(allele);
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
        return new Genotype(g.getSampleName(), g.getAlleles(), g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, attrs, g.isPhased());
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
        if ( vc.getChromosomeCount() == 0 )
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
        return VariantContext.modifyGenotypes(vc, null);
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

    public static VariantContext pruneVariantContext(VariantContext vc) {
        return pruneVariantContext(vc, null);
    }

    public static VariantContext pruneVariantContext(final VariantContext vc, final Collection<String> keysToPreserve ) {
        final MutableVariantContext mvc = new MutableVariantContext(vc);

        if ( keysToPreserve == null || keysToPreserve.size() == 0 )
            mvc.clearAttributes();
        else {
            final Map<String, Object> d = mvc.getAttributes();
            mvc.clearAttributes();
            for ( String key : keysToPreserve )
                if ( d.containsKey(key) )
                    mvc.putAttribute(key, d.get(key));
        }

        // this must be done as the ID is stored in the attributes field
        if ( vc.hasID() ) mvc.setID(vc.getID());

        Collection<Genotype> gs = mvc.getGenotypes().values();
        mvc.clearGenotypes();
        for ( Genotype g : gs ) {
            MutableGenotype mg = new MutableGenotype(g);
            mg.clearAttributes();
            if ( keysToPreserve != null )
                for ( String key : keysToPreserve )
                    if ( g.hasAttribute(key) )
                        mg.putAttribute(key, g.getAttribute(key));
            mvc.addGenotype(mg);
        }

        return mvc;
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
     * Performs a master merge on the VCs.  Here there is a master input [contains all of the information] and many
     * VCs containing partial, extra genotype information which should be added to the master.  For example,
     * we scatter out the phasing algorithm over some samples in the master, producing a minimal VCF with phasing
     * information per genotype.  The master merge will add the PQ information from each genotype record, where
     * appropriate, to the master VC.
     *
     * @param unsortedVCs   collection of VCs
     * @param masterName    name of master VC
     * @return  master-merged VC
     */
    public static VariantContext masterMerge(Collection<VariantContext> unsortedVCs, String masterName) {
        VariantContext master = findMaster(unsortedVCs, masterName);
        Map<String, Genotype> genotypes = master.getGenotypes();
        for (Genotype g : genotypes.values()) {
            genotypes.put(g.getSampleName(), new MutableGenotype(g));
        }

        Map<String, Object> masterAttributes = new HashMap<String, Object>(master.getAttributes());

        for (VariantContext vc : unsortedVCs) {
            if (!vc.getSource().equals(masterName)) {
                for (Genotype g : vc.getGenotypes().values()) {
                    MutableGenotype masterG = (MutableGenotype) genotypes.get(g.getSampleName());
                    for (Map.Entry<String, Object> attr : g.getAttributes().entrySet()) {
                        if (!masterG.hasAttribute(attr.getKey())) {
                            //System.out.printf("Adding GT attribute %s to masterG %s, new %s%n", attr, masterG, g);
                            masterG.putAttribute(attr.getKey(), attr.getValue());
                        }
                    }

                    if (masterG.isPhased() != g.isPhased()) {
                        if (masterG.sameGenotype(g)) {
                            // System.out.printf("Updating phasing %s to masterG %s, new %s%n", g.isPhased(), masterG, g);
                            masterG.setAlleles(g.getAlleles());
                            masterG.setPhase(g.isPhased());
                        }
                        //else System.out.println("WARNING: Not updating phase, since genotypes differ between master file and auxiliary info file!");
                    }

//                    if ( MathUtils.compareDoubles(masterG.getNegLog10PError(), g.getNegLog10PError()) != 0 ) {
//                        System.out.printf("Updating GQ %s to masterG %s, new %s%n", g.getNegLog10PError(), masterG, g);
//                        masterG.setNegLog10PError(g.getNegLog10PError());
//                    }

                }

                for (Map.Entry<String, Object> attr : vc.getAttributes().entrySet()) {
                    if (!masterAttributes.containsKey(attr.getKey())) {
                        //System.out.printf("Adding VC attribute %s to master %s, new %s%n", attr, master, vc);
                        masterAttributes.put(attr.getKey(), attr.getValue());
                    }
                }
            }
        }

        return new VariantContext(master.getSource(), master.getChr(), master.getStart(), master.getEnd(), master.getAlleles(), genotypes, master.getNegLog10PError(), master.getFilters(), masterAttributes);
    }

    private static VariantContext findMaster(Collection<VariantContext> unsortedVCs, String masterName) {
        for (VariantContext vc : unsortedVCs) {
            if (vc.getSource().equals(masterName)) {
                return vc;
            }
        }

        throw new ReviewedStingException(String.format("Couldn't find master VCF %s at %s", masterName, unsortedVCs.iterator().next()));
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
                VCs.add(VariantContext.createVariantContextWithPaddedAlleles(vc, false));
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
        double negLog10PError = -1;
        VariantContext vcWithMaxAC = null;
        Map<String, Genotype> genotypes = new TreeMap<String, Genotype>();

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

            negLog10PError = Math.max(negLog10PError, vc.isVariant() ? vc.getNegLog10PError() : -1);

            filters.addAll(vc.getFilters());

            //
            // add attributes
            //
            // special case DP (add it up) and ID (just preserve it)
            //
            if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            if ( vc.hasID() && ! vc.getID().equals(VCFConstants.EMPTY_ID_FIELD) ) rsIDs.add(vc.getID());
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

        if ( ! rsIDs.isEmpty() ) {
            attributes.put(VariantContext.ID_KEY, Utils.join(",", rsIDs));
        }

        VariantContext merged = new VariantContext(name, loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes, negLog10PError, filters, (mergeInfoWithMaxAC ? attributesWithMaxAC : attributes) );
        // Trim the padded bases of all alleles if necessary
        merged = createVariantContextWithTrimmedAlleles(merged);

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
            Map<String, Genotype> genotypes = new TreeMap<String, Genotype>();

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
            for ( Map.Entry<String, Genotype> sample : inputVC.getGenotypes().entrySet() ) {

                List<Allele> originalAlleles = sample.getValue().getAlleles();
                List<Allele> trimmedAlleles = new ArrayList<Allele>();
                for ( Allele a : originalAlleles ) {
                    if ( a.isCalled() )
                        trimmedAlleles.add(originalToTrimmedAlleleMap.get(a));
                    else
                        trimmedAlleles.add(Allele.NO_CALL);
                }
                genotypes.put(sample.getKey(), Genotype.modifyAlleles(sample.getValue(), trimmedAlleles));

            }
            return new VariantContext(inputVC.getSource(), inputVC.getChr(), inputVC.getStart(), inputVC.getEnd(), alleles, genotypes, inputVC.getNegLog10PError(), inputVC.filtersWereApplied() ? inputVC.getFilters() : null, attributes, new Byte(inputVC.getReference().getBases()[0]));

        }

        return inputVC;
    }

    public static Map<String, Genotype> stripPLs(Map<String, Genotype> genotypes) {
        Map<String, Genotype> newGs = new HashMap<String, Genotype>(genotypes.size());

        for ( Map.Entry<String, Genotype> g : genotypes.entrySet() ) {
            newGs.put(g.getKey(), g.getValue().hasLikelihoods() ? removePLs(g.getValue()) : g.getValue());
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

    private static void mergeGenotypes(Map<String, Genotype> mergedGenotypes, VariantContext oneVC, AlleleMapper alleleMapping, boolean uniqifySamples) {
        for ( Genotype g : oneVC.getGenotypes().values() ) {
            String name = mergedSampleName(oneVC.getSource(), g.getSampleName(), uniqifySamples);
            if ( ! mergedGenotypes.containsKey(name) ) {
                // only add if the name is new
                Genotype newG = g;

                if ( uniqifySamples || alleleMapping.needsRemapping() ) {
                    MutableGenotype mutG = new MutableGenotype(name, g);
                    if ( alleleMapping.needsRemapping() ) mutG.setAlleles(alleleMapping.remap(g.getAlleles()));
                    newG = mutG;
                }

                mergedGenotypes.put(name, newG);
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
        Map<String, Genotype> newGenotypes = new HashMap<String, Genotype>(vc.getNSamples());
        for ( Map.Entry<String, Genotype> genotype : vc.getGenotypes().entrySet() ) {
            List<Allele> newAlleles = new ArrayList<Allele>();
            for ( Allele allele : genotype.getValue().getAlleles() ) {
                Allele newAllele = alleleMap.get(allele);
                if ( newAllele == null )
                    newAllele = Allele.NO_CALL;
                newAlleles.add(newAllele);
            }
            newGenotypes.put(genotype.getKey(), Genotype.modifyAlleles(genotype.getValue(), newAlleles));
        }

        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), alleleMap.values(), newGenotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, vc.getAttributes());

    }

    public static VariantContext purgeUnallowedGenotypeAttributes(VariantContext vc, Set<String> allowedAttributes) {
        if ( allowedAttributes == null )
            return vc;

        Map<String, Genotype> newGenotypes = new HashMap<String, Genotype>(vc.getNSamples());
        for ( Map.Entry<String, Genotype> genotype : vc.getGenotypes().entrySet() ) {
            Map<String, Object> attrs = new HashMap<String, Object>();
            for ( Map.Entry<String, Object> attr : genotype.getValue().getAttributes().entrySet() ) {
                if ( allowedAttributes.contains(attr.getKey()) )
                    attrs.put(attr.getKey(), attr.getValue());
            }
            newGenotypes.put(genotype.getKey(), Genotype.modifyAttributes(genotype.getValue(), attrs));
        }

        return VariantContext.modifyGenotypes(vc, newGenotypes);
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

    public abstract static class AlleleMergeRule {
        // vc1, vc2 are ONLY passed to allelesShouldBeMerged() if mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2) AND allSamplesAreMergeable(vc1, vc2):
        abstract public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2);

        public String toString() {
            return "all samples are mergeable";
        }
    }

    // NOTE: returns null if vc1 and vc2 are not merged into a single MNP record

    public static VariantContext mergeIntoMNP(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile, AlleleMergeRule alleleMergeRule) {
        if (!mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2))
            return null;

        // Check that it's logically possible to merge the VCs:
        if (!allSamplesAreMergeable(vc1, vc2))
            return null;

        // Check if there's a "point" in merging the VCs (e.g., annotations could be changed)
        if (!alleleMergeRule.allelesShouldBeMerged(vc1, vc2))
            return null;

        return reallyMergeIntoMNP(vc1, vc2, referenceFile);
    }

    private static VariantContext reallyMergeIntoMNP(VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile) {
        int startInter = vc1.getEnd() + 1;
        int endInter = vc2.getStart() - 1;
        byte[] intermediateBases = null;
        if (startInter <= endInter) {
            intermediateBases = referenceFile.getSubsequenceAt(vc1.getChr(), startInter, endInter).getBases();
            StringUtil.toUpperCase(intermediateBases);
        }
        MergedAllelesData mergeData = new MergedAllelesData(intermediateBases, vc1, vc2); // ensures that the reference allele is added

        Map<String, Genotype> mergedGenotypes = new HashMap<String, Genotype>();
        for (Map.Entry<String, Genotype> gt1Entry : vc1.getGenotypes().entrySet()) {
            String sample = gt1Entry.getKey();
            Genotype gt1 = gt1Entry.getValue();
            Genotype gt2 = vc2.getGenotype(sample);

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            List<Allele> mergedAllelesForSample = new LinkedList<Allele>();

            /* NOTE: Since merged alleles are added to mergedAllelesForSample in the SAME order as in the input VC records,
               we preserve phase information (if any) relative to whatever precedes vc1:
             */
            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next(); // this is OK, since allSamplesAreMergeable()

                Allele mergedAllele = mergeData.ensureMergedAllele(all1, all2);
                mergedAllelesForSample.add(mergedAllele);
            }

            double mergedGQ = Math.max(gt1.getNegLog10PError(), gt2.getNegLog10PError());
            Set<String> mergedGtFilters = new HashSet<String>(); // Since gt1 and gt2 were unfiltered, the Genotype remains unfiltered

            Map<String, Object> mergedGtAttribs = new HashMap<String, Object>();
            PhaseAndQuality phaseQual = calcPhaseForMergedGenotypes(gt1, gt2);
            if (phaseQual.PQ != null)
                mergedGtAttribs.put(ReadBackedPhasingWalker.PQ_KEY, phaseQual.PQ);

            Genotype mergedGt = new Genotype(sample, mergedAllelesForSample, mergedGQ, mergedGtFilters, mergedGtAttribs, phaseQual.isPhased);
            mergedGenotypes.put(sample, mergedGt);
        }

        String mergedName = VariantContextUtils.mergeVariantContextNames(vc1.getSource(), vc2.getSource());
        double mergedNegLog10PError = Math.max(vc1.getNegLog10PError(), vc2.getNegLog10PError());
        Set<String> mergedFilters = new HashSet<String>(); // Since vc1 and vc2 were unfiltered, the merged record remains unfiltered
        Map<String, Object> mergedAttribs = VariantContextUtils.mergeVariantContextAttributes(vc1, vc2);

        VariantContext mergedVc = new VariantContext(mergedName, vc1.getChr(), vc1.getStart(), vc2.getEnd(), mergeData.getAllMergedAlleles(), mergedGenotypes, mergedNegLog10PError, mergedFilters, mergedAttribs);

        mergedAttribs = new HashMap<String, Object>(mergedVc.getAttributes());
        VariantContextUtils.calculateChromosomeCounts(mergedVc, mergedAttribs, true);
        mergedVc = VariantContext.modifyAttributes(mergedVc, mergedAttribs);

        return mergedVc;
    }

    private static class AlleleOneAndTwo {
        private Allele all1;
        private Allele all2;

        public AlleleOneAndTwo(Allele all1, Allele all2) {
            this.all1 = all1;
            this.all2 = all2;
        }

        public int hashCode() {
            return all1.hashCode() + all2.hashCode();
        }

        public boolean equals(Object other) {
            if (!(other instanceof AlleleOneAndTwo))
                return false;

            AlleleOneAndTwo otherAot = (AlleleOneAndTwo) other;
            return (this.all1.equals(otherAot.all1) && this.all2.equals(otherAot.all2));
        }
    }

    private static class MergedAllelesData {
        private Map<AlleleOneAndTwo, Allele> mergedAlleles;
        private byte[] intermediateBases;
        private int intermediateLength;

        public MergedAllelesData(byte[] intermediateBases, VariantContext vc1, VariantContext vc2) {
            this.mergedAlleles = new HashMap<AlleleOneAndTwo, Allele>(); // implemented equals() and hashCode() for AlleleOneAndTwo
            this.intermediateBases = intermediateBases;
            this.intermediateLength = this.intermediateBases != null ? this.intermediateBases.length : 0;

            this.ensureMergedAllele(vc1.getReference(), vc2.getReference(), true);
        }

        public Allele ensureMergedAllele(Allele all1, Allele all2) {
            return ensureMergedAllele(all1, all2, false); // false <-> since even if all1+all2 = reference, it was already created in the constructor
        }

        private Allele ensureMergedAllele(Allele all1, Allele all2, boolean creatingReferenceForFirstTime) {
            AlleleOneAndTwo all12 = new AlleleOneAndTwo(all1, all2);
            Allele mergedAllele = mergedAlleles.get(all12);

            if (mergedAllele == null) {
                byte[] bases1 = all1.getBases();
                byte[] bases2 = all2.getBases();

                byte[] mergedBases = new byte[bases1.length + intermediateLength + bases2.length];
                System.arraycopy(bases1, 0, mergedBases, 0, bases1.length);
                if (intermediateBases != null)
                    System.arraycopy(intermediateBases, 0, mergedBases, bases1.length, intermediateLength);
                System.arraycopy(bases2, 0, mergedBases, bases1.length + intermediateLength, bases2.length);

                mergedAllele = Allele.create(mergedBases, creatingReferenceForFirstTime);
                mergedAlleles.put(all12, mergedAllele);
            }

            return mergedAllele;
        }

        public Set<Allele> getAllMergedAlleles() {
            return new HashSet<Allele>(mergedAlleles.values());
        }
    }

    private static String mergeVariantContextNames(String name1, String name2) {
        return name1 + "_" + name2;
    }

    private static Map<String, Object> mergeVariantContextAttributes(VariantContext vc1, VariantContext vc2) {
        Map<String, Object> mergedAttribs = new HashMap<String, Object>();

        List<VariantContext> vcList = new LinkedList<VariantContext>();
        vcList.add(vc1);
        vcList.add(vc2);

        String[] MERGE_OR_ATTRIBS = {VCFConstants.DBSNP_KEY};
        for (String orAttrib : MERGE_OR_ATTRIBS) {
            boolean attribVal = false;
            for (VariantContext vc : vcList) {
                attribVal = vc.getAttributeAsBoolean(orAttrib, false);
                if (attribVal) // already true, so no reason to continue:
                    break;
            }
            mergedAttribs.put(orAttrib, attribVal);
        }

        // Merge ID fields:
        String iDVal = null;
        for (VariantContext vc : vcList) {
            String val = vc.getAttributeAsString(VariantContext.ID_KEY, null);
            if (val != null && !val.equals(VCFConstants.EMPTY_ID_FIELD)) {
                if (iDVal == null)
                    iDVal = val;
                else
                    iDVal += VCFConstants.ID_FIELD_SEPARATOR + val;
            }
        }
        if (iDVal != null)
            mergedAttribs.put(VariantContext.ID_KEY, iDVal);

        return mergedAttribs;
    }

    private static boolean mergeIntoMNPvalidationCheck(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = VariantContextUtils.getLocation(genomeLocParser, vc1);
        GenomeLoc loc2 = VariantContextUtils.getLocation(genomeLocParser, vc2);

        if (!loc1.onSameContig(loc2))
            throw new ReviewedStingException("Can only merge vc1, vc2 if on the same chromosome");

        if (!loc1.isBefore(loc2))
            throw new ReviewedStingException("Can only merge if vc1 is BEFORE vc2");

        if (vc1.isFiltered() || vc2.isFiltered())
            return false;

        if (!vc1.getSampleNames().equals(vc2.getSampleNames())) // vc1, vc2 refer to different sample sets
            return false;

        if (!allGenotypesAreUnfilteredAndCalled(vc1) || !allGenotypesAreUnfilteredAndCalled(vc2))
            return false;

        return true;
    }

    private static boolean allGenotypesAreUnfilteredAndCalled(VariantContext vc) {
        for (Map.Entry<String, Genotype> gtEntry : vc.getGenotypes().entrySet()) {
            Genotype gt = gtEntry.getValue();
            if (gt.isNoCall() || gt.isFiltered())
                return false;
        }

        return true;
    }

    // Assumes that vc1 and vc2 were already checked to have the same sample names:

    private static boolean allSamplesAreMergeable(VariantContext vc1, VariantContext vc2) {
        // Check that each sample's genotype in vc2 is uniquely appendable onto its genotype in vc1:
        for (Map.Entry<String, Genotype> gt1Entry : vc1.getGenotypes().entrySet()) {
            String sample = gt1Entry.getKey();
            Genotype gt1 = gt1Entry.getValue();
            Genotype gt2 = vc2.getGenotype(sample);

            if (!alleleSegregationIsKnown(gt1, gt2)) // can merge if: phased, or if either is a hom
                return false;
        }

        return true;
    }

    public static boolean alleleSegregationIsKnown(Genotype gt1, Genotype gt2) {
        if (gt1.getPloidy() != gt2.getPloidy())
            return false;

        /* If gt2 is phased or hom, then could even be MERGED with gt1 [This is standard].

           HOWEVER, EVEN if this is not the case, but gt1.isHom(),
           it is trivially known that each of gt2's alleles segregate with the single allele type present in gt1.
         */
        return (gt2.isPhased() || gt2.isHom() || gt1.isHom());
    }

    private static class PhaseAndQuality {
        public boolean isPhased;
        public Double PQ = null;

        public PhaseAndQuality(Genotype gt) {
            this.isPhased = gt.isPhased();
            if (this.isPhased) {
                this.PQ = gt.getAttributeAsDouble(ReadBackedPhasingWalker.PQ_KEY, -1);
                if ( this.PQ == -1 ) this.PQ = null;
            }
        }
    }

    // Assumes that alleleSegregationIsKnown(gt1, gt2):

    private static PhaseAndQuality calcPhaseForMergedGenotypes(Genotype gt1, Genotype gt2) {
        if (gt2.isPhased() || gt2.isHom())
            return new PhaseAndQuality(gt1); // maintain the phase of gt1

        if (!gt1.isHom())
            throw new ReviewedStingException("alleleSegregationIsKnown(gt1, gt2) implies: gt2.genotypesArePhased() || gt2.isHom() || gt1.isHom()");

        /* We're dealing with: gt1.isHom(), gt2.isHet(), !gt2.genotypesArePhased(); so, the merged (het) Genotype is not phased relative to the previous Genotype

           For example, if we're merging the third Genotype with the second one:
           0/1
           1|1
           0/1

           Then, we want to output:
           0/1
           1/2
         */
        return new PhaseAndQuality(gt2); // maintain the phase of gt2 [since !gt2.genotypesArePhased()]
    }

    /* Checks if any sample has a MNP of ALT alleles (segregating together):
     [Assumes that vc1 and vc2 were already checked to have the same sample names && allSamplesAreMergeable(vc1, vc2)]
     */

    public static boolean someSampleHasDoubleNonReferenceAllele(VariantContext vc1, VariantContext vc2) {
        for (Map.Entry<String, Genotype> gt1Entry : vc1.getGenotypes().entrySet()) {
            String sample = gt1Entry.getKey();
            Genotype gt1 = gt1Entry.getValue();
            Genotype gt2 = vc2.getGenotype(sample);

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next(); // this is OK, since allSamplesAreMergeable()

                if (all1.isNonReference() && all2.isNonReference()) // corresponding alleles are alternate
                    return true;
            }
        }

        return false;
    }

    /* Checks if all samples are consistent in their haplotypes:
     [Assumes that vc1 and vc2 were already checked to have the same sample names && allSamplesAreMergeable(vc1, vc2)]
     */

    public static boolean doubleAllelesSegregatePerfectlyAmongSamples(VariantContext vc1, VariantContext vc2) {
        // Check that Alleles at vc1 and at vc2 always segregate together in all samples (including reference):
        Map<Allele, Allele> allele1ToAllele2 = new HashMap<Allele, Allele>();
        Map<Allele, Allele> allele2ToAllele1 = new HashMap<Allele, Allele>();

        // Note the segregation of the alleles for the reference genome:
        allele1ToAllele2.put(vc1.getReference(), vc2.getReference());
        allele2ToAllele1.put(vc2.getReference(), vc1.getReference());

        // Note the segregation of the alleles for each sample (and check that it is consistent with the reference and all previous samples).
        for (Map.Entry<String, Genotype> gt1Entry : vc1.getGenotypes().entrySet()) {
            String sample = gt1Entry.getKey();
            Genotype gt1 = gt1Entry.getValue();
            Genotype gt2 = vc2.getGenotype(sample);

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next();

                Allele all1To2 = allele1ToAllele2.get(all1);
                if (all1To2 == null)
                    allele1ToAllele2.put(all1, all2);
                else if (!all1To2.equals(all2)) // all1 segregates with two different alleles at site 2
                    return false;

                Allele all2To1 = allele2ToAllele1.get(all2);
                if (all2To1 == null)
                    allele2ToAllele1.put(all2, all1);
                else if (!all2To1.equals(all1)) // all2 segregates with two different alleles at site 1
                    return false;
            }
        }

        return true;
    }
}