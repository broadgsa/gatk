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

package org.broadinstitute.sting.gatk.contexts.variantcontext;

import java.io.Serializable;
import java.util.*;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.StringUtil;
import org.apache.commons.jexl2.*;
import org.broad.tribble.util.popgen.HardyWeinbergCalculation;
import org.broad.tribble.util.variantcontext.*;
import org.broadinstitute.sting.gatk.walkers.phasing.ReadBackedPhasingWalker;
import org.broadinstitute.sting.utils.*;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

public class VariantContextUtils {
    final public static JexlEngine engine = new JexlEngine();

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

            attributes.put(VCFConstants.ALLELE_COUNT_KEY, alleleCounts);
            attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs);
        }
        else {
            attributes.put(VCFConstants.ALLELE_COUNT_KEY, 0);
            attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0);
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
    public static boolean match(GenomeLocParser genomeLocParser,VariantContext vc, JexlVCMatchExp exp) {
        return match(genomeLocParser,vc,Arrays.asList(exp)).get(exp);
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
    public static Map<JexlVCMatchExp, Boolean> match(GenomeLocParser genomeLocParser,VariantContext vc, Collection<JexlVCMatchExp> exps) {
        return new JEXLMap(genomeLocParser,exps,vc);

    }

    /**
     * Returns true if exp match VC/g.  See collection<> version for full docs.
     * @param vc   variant context
     * @param g    genotype
     * @param exp   expression
     * @return true if there is a match
     */
    public static boolean match(GenomeLocParser genomeLocParser,VariantContext vc, Genotype g, JexlVCMatchExp exp) {
        return match(genomeLocParser,vc,g,Arrays.asList(exp)).get(exp);
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
    public static Map<JexlVCMatchExp, Boolean> match(GenomeLocParser genomeLocParser,VariantContext vc, Genotype g, Collection<JexlVCMatchExp> exps) {
        return new JEXLMap(genomeLocParser,exps,vc,g);

    }

    public static double computeHardyWeinbergPvalue(VariantContext vc) {
        if ( vc.getChromosomeCount() == 0 )
            return 0.0;
        return HardyWeinbergCalculation.hwCalculate(vc.getHomRefCount(), vc.getHetCount(), vc.getHomVarCount());
    }

    public static VariantContext pruneVariantContext(VariantContext vc) {
        return pruneVariantContext(vc, null);
    }

    public static VariantContext pruneVariantContext(VariantContext vc, Set<String> keysToPreserve ) {
        MutableVariantContext mvc = new MutableVariantContext(vc);

        if ( keysToPreserve == null || keysToPreserve.size() == 0 )
            mvc.clearAttributes();
        else {
            Map<String, Object> d = mvc.getAttributes();
            mvc.clearAttributes();
            for ( String key : keysToPreserve )
                mvc.putAttribute(key, d.get(key));
        }

        Collection<Genotype> gs = mvc.getGenotypes().values();
        mvc.clearGenotypes();
        for ( Genotype g : gs ) {
            MutableGenotype mg = new MutableGenotype(g);
            mg.clearAttributes();
            mvc.addGenotype(mg);
        }

        return mvc;
    }

    public enum GenotypeMergeType {
        UNIQUIFY, PRIORITIZE, UNSORTED, REQUIRE_UNIQUE
    }

    public enum VariantMergeType {
        UNION, INTERSECT
    }

    public static VariantContext simpleMerge(GenomeLocParser genomeLocParser, Collection<VariantContext> unsortedVCs, byte refBase) {
        return simpleMerge(genomeLocParser, unsortedVCs, null, VariantMergeType.INTERSECT, GenotypeMergeType.UNSORTED, false, false, refBase);
    }


    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniqifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name
     *
     * @param unsortedVCs
     * @param priorityListOfVCs
     * @param variantMergeOptions
     * @param genotypeMergeOptions
     * @return
     */
    public static VariantContext simpleMerge(GenomeLocParser genomeLocParser, Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs,
                                             VariantMergeType variantMergeOptions, GenotypeMergeType genotypeMergeOptions,
                                             boolean annotateOrigin, boolean printMessages, byte inputRefBase ) {

        return simpleMerge(genomeLocParser, unsortedVCs, priorityListOfVCs, variantMergeOptions, genotypeMergeOptions, annotateOrigin, printMessages, inputRefBase, "set", false);
    }

    public static VariantContext simpleMerge(GenomeLocParser genomeLocParser, Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs,
                                             VariantMergeType variantMergeOptions, GenotypeMergeType genotypeMergeOptions,
                                             boolean annotateOrigin, boolean printMessages, byte inputRefBase, String setKey,
                                             boolean filteredAreUncalled ) {
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
                VCs.add(VariantContext.createVariantContextWithPaddedAlleles(vc,inputRefBase,false));            
        }
        if ( VCs.size() == 0 ) // everything is filtered out and we're filteredareUncalled
            return null;

        // establish the baseline info from the first VC
        VariantContext first = VCs.get(0);
        String name = first.getSource();
        GenomeLoc loc = getLocation(genomeLocParser,first);

        Set<Allele> alleles = new TreeSet<Allele>();
        Map<String, Genotype> genotypes = new TreeMap<String, Genotype>();
        double negLog10PError = -1;
        Set<String> filters = new TreeSet<String>();
        Map<String, Object> attributes = new TreeMap<String, Object>();
        Set<String> inconsistentAttributes = new HashSet<String>();
        String rsID = null;
        int depth = 0;

        // counting the number of filtered and variant VCs
        int nFiltered = 0, nVariant = 0;

        Allele refAllele = determineReferenceAllele(VCs);
        boolean remapped = false;

        // cycle through and add info from the other VCs, making sure the loc/reference matches

        for ( VariantContext vc : VCs ) {
            if ( loc.getStart() != vc.getStart() ) // || !first.getReference().equals(vc.getReference()) )
                throw new ReviewedStingException("BUG: attempting to merge VariantContexts with different start sites: first="+ first.toString() + " second=" + vc.toString());

            if ( getLocation(genomeLocParser,vc).size() > loc.size() )
                loc = getLocation(genomeLocParser,vc); // get the longest location

            nFiltered += vc.isFiltered() ? 1 : 0;
            nVariant += vc.isVariant() ? 1 : 0;

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
            if ( vc.hasAttribute(VCFConstants.DEPTH_KEY) )
                depth += Integer.valueOf(vc.getAttributeAsString(VCFConstants.DEPTH_KEY));
            if ( rsID == null && vc.hasID() )
                rsID = vc.getID();

            for ( Map.Entry<String, Object> p : vc.getAttributes().entrySet() ) {
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

        // if at least one record was unfiltered and we want a union, clear all of the filters
        if ( variantMergeOptions == VariantMergeType.UNION && nFiltered != VCs.size() )
            filters.clear();

        // we care about where the call came from
        if ( annotateOrigin ) {
            String setValue;
            if ( nFiltered == 0 && nVariant == priorityListOfVCs.size() )                   // nothing was unfiltered
                setValue = "Intersection";
            else if ( nFiltered == VCs.size() )     // everything was filtered out
                setValue = "FilteredInAll";
            else if ( nVariant == 0 )               // everyone was reference
                setValue = "ReferenceInAll";
            else {                                  // we are filtered in some subset
                List<String> s = new ArrayList<String>();
                for ( VariantContext vc : VCs )
                    if ( vc.isVariant() )
                        s.add( vc.isFiltered() ? "filterIn" + vc.getSource() : vc.getSource() );
                setValue = Utils.join("-", s);
            }
            
            if ( setKey != null ) attributes.put(setKey, setValue);
        }

        if ( depth > 0 )
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
        if ( rsID != null )
            attributes.put(VariantContext.ID_KEY, rsID);

        VariantContext merged = new VariantContext(name, loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes, negLog10PError, filters, attributes);
        if ( printMessages && remapped ) System.out.printf("Remapped => %s%n", merged);
        return merged;
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


    public static VariantContext createVariantContextWithTrimmedAlleles(VariantContext inputVC) {
        // see if we need to trim common reference base from all alleles
        boolean trimVC = true;

        // We need to trim common reference base from all alleles if a ref base is common to all alleles
        Allele refAllele = inputVC.getReference();
        if (!inputVC.isVariant())
            trimVC = false;
        else if (refAllele.isNull())
            trimVC = false;
        else {
            for (Allele a : inputVC.getAlternateAlleles()) {
                if (a.length() < 1 || (a.getBases()[0] != refAllele.getBases()[0]))
                    trimVC = false;
            }
        }

        // nothing to do if we don't need to trim bases
        if (trimVC) {
            List<Allele> alleles = new ArrayList<Allele>();
            Map<String, Genotype> genotypes = new TreeMap<String, Genotype>();

            Map<String, Genotype> inputGenotypes = inputVC.getGenotypes();
            // set the reference base for indels in the attributes
            Map<String,Object> attributes = new TreeMap<String,Object>();

            for ( Map.Entry<String, Object> p : inputVC.getAttributes().entrySet() ) {
                attributes.put(p.getKey(), p.getValue());
            }

            attributes.put(VariantContext.REFERENCE_BASE_FOR_INDEL_KEY, new Byte(inputVC.getReference().getBases()[0]));


            for (Allele a : inputVC.getAlleles()) {
                // get bases for current allele and create a new one with trimmed bases
                byte[] newBases = Arrays.copyOfRange(a.getBases(),1,a.length());
                alleles.add(Allele.create(newBases,a.isReference()));
            }

            // now we can recreate new genotypes with trimmed alleles
            for (String sample : inputVC.getSampleNames()) {
                Genotype g = inputGenotypes.get(sample);

                List<Allele> inAlleles = g.getAlleles();
                List<Allele> newGenotypeAlleles = new ArrayList<Allele>();
                for (Allele a : inAlleles) {
                    byte[] newBases = Arrays.copyOfRange(a.getBases(),1,a.length());
                    newGenotypeAlleles.add(Allele.create(newBases, a.isReference()));
                }
                genotypes.put(sample, new Genotype(sample, newGenotypeAlleles, g.getNegLog10PError(),
                        g.getFilters(),g.getAttributes(),g.genotypesArePhased()));

            }
            return new VariantContext(inputVC.getSource(), inputVC.getChr(), inputVC.getStart(), inputVC.getEnd(), alleles, genotypes, inputVC.getNegLog10PError(),
                    inputVC.getFilters(), attributes);

        }
        else
            return inputVC;

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

    public static VariantContext modifyLocation(VariantContext vc, GenomeLoc loc) {
        return new VariantContext(vc.getSource(), loc.getContig(), loc.getStart(), loc.getStop(), vc.getAlleles(), vc.getGenotypes(), vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, vc.getAttributes());
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
     * @param vc the variant context
     * @return the genomeLoc
     */
    public static final GenomeLoc getLocation(GenomeLocParser genomeLocParser,VariantContext vc) {
        return genomeLocParser.createGenomeLoc(vc.getChr(),(int)vc.getStart(),(int)vc.getEnd());
    }

    // NOTE: returns null if vc1 and vc2 are not mergeable into a single MNP record
    public static VariantContext mergeIntoMNP(GenomeLocParser genomeLocParser,VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile) {
        if (!mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2))
            return null;

        // Check that it's logically possible to merge the VCs, and that there's a point in doing so (e.g., annotations could be changed):
        if (!allSamplesAreMergeable(vc1, vc2) || !someSampleHasDoubleNonReferenceAllele(vc1, vc2))
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
        Map<String, Object> mergedAttribs = VariantContextUtils.mergeVariantContextAttributes(vc1.getAttributes(), vc2.getAttributes());
        if (mergedAttribs == null)
            return null;

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

            AlleleOneAndTwo otherAot = (AlleleOneAndTwo)other;
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

    private static Map<String, Object> mergeVariantContextAttributes(Map<String, Object> attribs1, Map<String, Object> attribs2) {
        Map<String, Object> mergedAttribs = new HashMap<String, Object>();

        List<Map<String, Object>> attribsList = new LinkedList<Map<String, Object>>();
        attribsList.add(attribs1);
        attribsList.add(attribs2);

        String[] MERGE_OR_ATTRIBS = {VCFConstants.DBSNP_KEY};
        for (String orAttrib : MERGE_OR_ATTRIBS) {
            boolean attribVal = false;
            for (Map<String, Object> attribs : attribsList) {
                Boolean val = getBooleanAttribute(attribs, orAttrib);
                if (val != null)
                    attribVal = (attribVal || val);
                if (attribVal) // already true, so no reason to continue:
                    break;
            }
            mergedAttribs.put(orAttrib, attribVal);
        }

        // Merge ID fields:
        String iDVal = null;
        for (Map<String, Object> attribs : attribsList) {
            String val = getStringAttribute(attribs, VariantContext.ID_KEY);
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

    public static Boolean getBooleanAttribute(Map<String, Object> attribs, String attribName) {
        Object val = attribs.get(attribName);
        if (val == null || val.equals(VCFConstants.MISSING_VALUE_v4))
            return null;

        try {
            return new Boolean(val.toString());
        }
        catch (Exception e) {// IGNORE unparseable data
            return null;
        }
    }

    public static String getStringAttribute(Map<String, Object> attribs, String attribName) {
        Object val = attribs.get(attribName);
        if (val == null || val.equals(VCFConstants.MISSING_VALUE_v4))
            return null;

        try {
            return val.toString();
        }
        catch (Exception e) {// IGNORE unparseable data
            return null;
        }
    }

    public static Double getDoubleAttribute(Map<String, Object> attribs, String attribName) {
        Object val = attribs.get(attribName);
        if (val == null || val.equals(VCFConstants.MISSING_VALUE_v4))
            return null;

        try {
            return new Double(val.toString());
        }
        catch (Exception e) {// IGNORE unparseable data
            return null;
        }
    }

    private static boolean mergeIntoMNPvalidationCheck(GenomeLocParser genomeLocParser,VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = VariantContextUtils.getLocation(genomeLocParser,vc1);
        GenomeLoc loc2 = VariantContextUtils.getLocation(genomeLocParser,vc2);

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
        return (gt2.genotypesArePhased() || gt2.isHom() || gt1.isHom());
    }

    private static class PhaseAndQuality {
        public boolean isPhased;
        public Double PQ = null;

        public PhaseAndQuality(Genotype gt) {
            this.isPhased = gt.genotypesArePhased();
            if (this.isPhased)
                this.PQ = VariantContextUtils.getDoubleAttribute(gt.getAttributes(), ReadBackedPhasingWalker.PQ_KEY);
        }
    }

    // Assumes that alleleSegregationIsKnown(gt1, gt2):
    private static PhaseAndQuality calcPhaseForMergedGenotypes(Genotype gt1, Genotype gt2) {
        if (gt2.genotypesArePhased() || gt2.isHom())
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
    private static boolean someSampleHasDoubleNonReferenceAllele(VariantContext vc1, VariantContext vc2) {
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
}