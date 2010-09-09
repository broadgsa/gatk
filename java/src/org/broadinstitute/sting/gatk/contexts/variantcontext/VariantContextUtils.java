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
import org.apache.commons.jexl2.*;
import org.broad.tribble.util.popgen.HardyWeinbergCalculation;
import org.broad.tribble.util.variantcontext.*;
import org.broadinstitute.sting.utils.*;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.UserError;

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
     */
    public static VariantContext toVC(String name, GenomeLoc loc, Collection<Allele> alleles, Collection<Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        return new VariantContext(name, loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes != null ? VariantContext.genotypeCollectionToMap(new TreeMap<String, Genotype>(), genotypes) : null, negLog10PError, filters, attributes);
    }

    /**
     * Create a new variant context without genotypes and no Perror, no filters, and no attributes
     * @param name            name
     * @param loc             location
     * @param alleles         alleles
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
     */
    public static VariantContext toVC(String name, GenomeLoc loc, Collection<Allele> alleles, Collection<Genotype> genotypes) {
        return new VariantContext(name, loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    /**
     * Copy constructor
     *
     * @param other the VariantContext to copy
     */
    public static VariantContext toVC(VariantContext other) {
        return new VariantContext(other.getName(), other.getChr(), other.getStart(), other.getEnd(), other.getAlleles(), other.getGenotypes(), other.getNegLog10PError(), other.getFilters(), other.getAttributes());
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
            throw new StingException("BUG: neither names nor exps can be null: names " + Arrays.toString(names) + " exps=" + Arrays.toString(exps) );

        if ( names.length != exps.length )
            throw new UserError("Inconsistent number of provided filter names and expressions: names=" + Arrays.toString(names) + " exps=" + Arrays.toString(exps));

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
		throw new UserError.BadArgumentValue(name, "Invalid expression used (" + expStr + "). Please see the JEXL docs for correct syntax.") ;
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

    public static VariantContext simpleMerge(Collection<VariantContext> unsortedVCs, byte refBase) {
        return simpleMerge(unsortedVCs, null, VariantMergeType.INTERSECT, GenotypeMergeType.UNSORTED, false, false, refBase);
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
    public static VariantContext simpleMerge(Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs,
                                             VariantMergeType variantMergeOptions, GenotypeMergeType genotypeMergeOptions,
                                             boolean annotateOrigin, boolean printMessages, byte inputRefBase ) {

        return simpleMerge(unsortedVCs, priorityListOfVCs, variantMergeOptions, genotypeMergeOptions, annotateOrigin, printMessages, inputRefBase, "set", false);
    }

    public static VariantContext simpleMerge(Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs,
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
                VCs.add(createVariantContextWithPaddedAlleles(vc,inputRefBase));
        }

        if ( VCs.size() == 0 ) // everything is filtered out and we're filteredareUncalled
            return null;

        // establish the baseline info from the first VC
        VariantContext first = VCs.get(0);
        String name = first.getName();
        GenomeLoc loc = getLocation(first);

        Set<Allele> alleles = new TreeSet<Allele>();
        Map<String, Genotype> genotypes = new TreeMap<String, Genotype>();
        double negLog10PError = -1;
        Set<String> filters = new TreeSet<String>();
        Map<String, Object> attributes = new TreeMap<String, Object>();
        Set<String> inconsistentAttributes = new HashSet<String>();
        String rsID = null;
        int depth = 0;

        // counting the number of filtered and varint VCs
        int nFiltered = 0, nVariant = 0;

        Allele refAllele = determineReferenceAllele(VCs);
        boolean remapped = false;

        // cycle through and add info from the other VCs, making sure the loc/reference matches

        for ( VariantContext vc : VCs ) {
            if ( loc.getStart() != vc.getStart() ) // || !first.getReference().equals(vc.getReference()) )
                throw new StingException("BUG: attempting to merge VariantContexts with different start sites: first="+ first.toString() + " second=" + vc.toString());

            if ( getLocation(vc).size() > loc.size() )
                loc = getLocation(vc); // get the longest location

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
            if ( rsID == null && vc.hasAttribute(VariantContext.ID_KEY) )
                rsID = vc.getAttributeAsString(VariantContext.ID_KEY);

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
                        s.add( vc.isFiltered() ? "filterIn" + vc.getName() : vc.getName() );
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
                    throw new UserError("REQUIRE_UNIQUE sample names is true but duplicate names were discovered " + name);
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
                throw new StingException("BUG: equal length references with difference bases: "+ ref + " " + myRef);
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
            if ( refAllele.length() <= myRef.length() ) throw new StingException("BUG: myRef="+myRef+" is longer than refAllele="+refAllele);
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
            return new VariantContext(inputVC.getName(), inputVC.getChr(), inputVC.getStart(), inputVC.getEnd(), alleles, genotypes, inputVC.getNegLog10PError(),
                    inputVC.getFilters(), attributes);

        }
        else
            return inputVC;

    }

    public static VariantContext createVariantContextWithPaddedAlleles(VariantContext inputVC, byte inputRefBase) {
        Allele refAllele = inputVC.getReference();


        // see if we need to pad common reference base from all alleles
        boolean padVC;

        // We need to pad a VC with a common base if the reference allele length is less than the vc location span.
        long locLength = getLocation(inputVC).size();
        if (refAllele.length() == locLength)
            padVC = false;
        else if (refAllele.length() == locLength-1)
            padVC = true;
        else throw new StingException("Badly formed variant context, reference length must be at most one base shorter than location size");


        // nothing to do if we don't need to pad bases
        if (padVC) {
            Byte refByte;

            Map<String,Object> attributes = inputVC.getAttributes();

            if (BaseUtils.isRegularBase(inputRefBase))
                refByte = inputRefBase;
            else if (attributes.containsKey(VariantContext.REFERENCE_BASE_FOR_INDEL_KEY))
                refByte = (Byte)attributes.get(VariantContext.REFERENCE_BASE_FOR_INDEL_KEY);
            else
                throw new StingException("Error when trying to pad Variant Context: either input reference base must be a regular base, or input VC must contain reference base key");

            List<Allele> alleles = new ArrayList<Allele>();
            Map<String, Genotype> genotypes = new TreeMap<String, Genotype>();

            Map<String, Genotype> inputGenotypes = inputVC.getGenotypes();

            for (Allele a : inputVC.getAlleles()) {
                // get bases for current allele and create a new one with trimmed bases
                String newBases = new String(new byte[]{refByte}) + a.getBaseString();
                alleles.add(Allele.create(newBases,a.isReference()));
            }

            // now we can recreate new genotypes with trimmed alleles
            for (String sample : inputVC.getSampleNames()) {
                Genotype g = inputGenotypes.get(sample);

                List<Allele> inAlleles = g.getAlleles();
                List<Allele> newGenotypeAlleles = new ArrayList<Allele>();
                for (Allele a : inAlleles) {
                    if (a.isCalled()) {
                        String newBases = new String(new byte[]{refByte}) + a.getBaseString();
                        newGenotypeAlleles.add(Allele.create(newBases,a.isReference()));
                    }
                    else {
                        // add no-call allele
                        newGenotypeAlleles.add(Allele.NO_CALL);
                    }
                }
                genotypes.put(sample, new Genotype(sample, newGenotypeAlleles, g.getNegLog10PError(),
                        g.getFilters(),g.getAttributes(),g.genotypesArePhased()));

            }
            return new VariantContext(inputVC.getName(), inputVC.getChr(), inputVC.getStart(), inputVC.getEnd(), alleles, genotypes, inputVC.getNegLog10PError(),
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
            int i = priorityListOfVCs.indexOf(vc.getName());
            if ( i == -1 ) throw new UserError.BadArgumentValue(Utils.join(",", priorityListOfVCs), "Priority list " + priorityListOfVCs + " doesn't contain variant context " + vc.getName());
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
            String name = mergedSampleName(oneVC.getName(), g.getSampleName(), uniqifySamples);
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

    public static VariantContext modifyGenotypes(VariantContext vc, Map<String, Genotype> genotypes) {
        return new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, vc.getAttributes());
    }

    public static VariantContext modifyLocation(VariantContext vc, GenomeLoc loc) {
        return new VariantContext(vc.getName(), loc.getContig(), loc.getStart(), loc.getStop(), vc.getAlleles(), vc.getGenotypes(), vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, vc.getAttributes());
    }

    public static VariantContext modifyFilters(VariantContext vc, Set<String> filters) {
        return new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd() , vc.getAlleles(), vc.getGenotypes(), vc.getNegLog10PError(), filters, vc.getAttributes());
    }

    public static VariantContext modifyAttributes(VariantContext vc, Map<String, Object> attributes) {
        return new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), vc.getGenotypes(), vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, attributes);
    }

    public static Genotype modifyName(Genotype g, String name) {
        return new Genotype(name, g.getAlleles(), g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, g.getAttributes(), g.genotypesArePhased());
    }

    public static Genotype modifyAttributes(Genotype g, Map<String, Object> attributes) {
        return new Genotype(g.getSampleName(), g.getAlleles(), g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, attributes, g.genotypesArePhased());
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
            newGenotypes.put(genotype.getKey(), VariantContextUtils.modifyAttributes(genotype.getValue(), attrs));
        }

        return VariantContextUtils.modifyGenotypes(vc, newGenotypes);
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
    public static final GenomeLoc getLocation(VariantContext vc) {
        return GenomeLocParser.createGenomeLoc(vc.getChr(),(int)vc.getStart(),(int)vc.getEnd());
    }

}
