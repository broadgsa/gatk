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

package org.broadinstitute.gatk.utils.genotyper;


import com.google.java.contract.Ensures;
import org.broadinstitute.gatk.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import htsjdk.variant.variantcontext.Allele;

import java.util.*;

/**
 *   Wrapper class that holds a set of maps of the form (Read -> Map(Allele->Double))
 *   For each read, this holds underlying alleles represented by an aligned read, and corresponding relative likelihood.
 */
public class PerReadAlleleLikelihoodMap {
    /** A set of all of the allele, so we can efficiently determine if an allele is already present */
    private final Map<Allele,Integer> allelesSet = new HashMap<>();
    /** A list of the unique allele, as an ArrayList so we can call get(i) efficiently */
    protected final List<Allele> alleles = new ArrayList<>();



    protected final Map<GATKSAMRecord, Map<Allele, Double>> likelihoodReadMap = new LinkedHashMap<>();

    public PerReadAlleleLikelihoodMap() { }

    /**
     * Add a new entry into the Read -> ( Allele -> Likelihood ) map of maps.
     * @param read - the GATKSAMRecord that was evaluated
     * @param a - the Allele against which the GATKSAMRecord was evaluated
     * @param likelihood - the likelihood score resulting from the evaluation of "read" against "a"
     */
    public void add(final GATKSAMRecord read, final Allele a, final Double likelihood) {
        if ( read == null ) throw new IllegalArgumentException("Cannot add a null read to the allele likelihood map");
        if ( a == null ) throw new IllegalArgumentException("Cannot add a null allele to the allele likelihood map");
        if ( likelihood == null ) throw new IllegalArgumentException("Likelihood cannot be null");
        if ( likelihood > 0.0 ) throw new IllegalArgumentException("Likelihood must be negative (L = log(p))");

        if (!allelesSet.containsKey(a)) {
            allelesSet.put(a,alleles.size());
            alleles.add(a);
        }
        Map<Allele,Double> likelihoodMap = likelihoodReadMap.get(read);
        if (likelihoodMap == null){
            // LinkedHashMap will ensure iterating through alleles will be in consistent order
            likelihoodMap = new LinkedHashMap<>();
            likelihoodReadMap.put(read,likelihoodMap);
        }

        likelihoodMap.put(a,likelihood);


    }

    public ReadBackedPileup createPerAlleleDownsampledBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction) {
        return AlleleBiasedDownsamplingUtils.createAlleleBiasedBasePileup(pileup, downsamplingFraction);
    }

    /**
     * For each allele "a" , identify those reads whose most likely allele is "a", and remove a "downsamplingFraction" proportion
     * of those reads from the "likelihoodReadMap". This is used for e.g. sample contamination
     * @param downsamplingFraction - the fraction of supporting reads to remove from each allele. If <=0 all reads kept, if >=1 all reads tossed.
     */
    public void performPerAlleleDownsampling(final double downsamplingFraction) {
        // special case removal of all or no reads
        if ( downsamplingFraction <= 0.0 )
            return;
        if ( downsamplingFraction >= 1.0 ) {
            likelihoodReadMap.clear();
            return;
        }

        // start by stratifying the reads by the alleles they represent at this position
        final Map<Allele, List<GATKSAMRecord>> alleleReadMap = getAlleleStratifiedReadMap();

        // compute the reads to remove and actually remove them
        final List<GATKSAMRecord> readsToRemove = AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(alleleReadMap, downsamplingFraction);
        for ( final GATKSAMRecord read : readsToRemove )
            likelihoodReadMap.remove(read);
    }

    /**
     * Convert the @likelihoodReadMap to a map of alleles to reads, where each read is mapped uniquely to the allele
     * for which it has the greatest associated likelihood
     * @return a map from each allele to a list of reads that 'support' the allele
     */
    protected Map<Allele,List<GATKSAMRecord>> getAlleleStratifiedReadMap() {
        final Map<Allele, List<GATKSAMRecord>> alleleReadMap = new HashMap<>(alleles.size());
        for ( final Allele allele : alleles )
            alleleReadMap.put(allele, new ArrayList<GATKSAMRecord>());

        for ( final Map.Entry<GATKSAMRecord, Map<Allele, Double>> entry : likelihoodReadMap.entrySet() ) {
            final MostLikelyAllele bestAllele = getMostLikelyAllele(entry.getValue());
            if ( bestAllele.isInformative() )
                alleleReadMap.get(bestAllele.getMostLikelyAllele()).add(entry.getKey());
        }

        return alleleReadMap;
    }

    @Ensures("result >=0")
    public int size() {
        return likelihoodReadMap.size();
    }

    /**
     * Helper function to add the read underneath a pileup element to the map
     * @param p                              Pileup element
     * @param a                              Corresponding allele
     * @param likelihood                     Allele likelihood
     */
    public void add(PileupElement p, Allele a, Double likelihood) {
        if (p==null)
            throw new IllegalArgumentException("Pileup element cannot be null");
        if ( p.getRead()==null )
           throw new IllegalArgumentException("Read underlying pileup element cannot be null");
        if ( a == null )
           throw new IllegalArgumentException("Allele for add() cannot be null");

        add(p.getRead(), a, likelihood);
    }

     /**
     * Does the current map contain the key associated with a particular SAM record in pileup?
     * @param p                 Pileup element
     * @return true if the map contains pileup element, else false
     */
    public boolean containsPileupElement(final PileupElement p) {
        return likelihoodReadMap.containsKey(p.getRead());
    }

    public boolean isEmpty() {
        return likelihoodReadMap.isEmpty();
    }

    public Map<GATKSAMRecord,Map<Allele,Double>> getLikelihoodReadMap() {
        return likelihoodReadMap;
    }

    public void clear() {
        allelesSet.clear();
        alleles.clear();
        likelihoodReadMap.clear();
    }

    public Set<GATKSAMRecord> getStoredElements() {
        return likelihoodReadMap.keySet();
    }

//    public Collection<Map<Allele,Double>> getLikelihoodMapValues() {
//        return likelihoodReadMap.values();
//    }

    public int getNumberOfStoredElements() {
        return likelihoodReadMap.size();
    }

    public Map<Allele,Double> getLikelihoodsAssociatedWithPileupElement(final PileupElement p) {
        if (!likelihoodReadMap.containsKey(p.getRead()))
            return null;

        return likelihoodReadMap.get(p.getRead());
    }


    /**
     * Get the log10 likelihood associated with an individual read/allele
     *
     * @param read the read whose likelihood we want
     * @param allele the allele whose likelihood we want
     * @return the log10 likelihood that this read matches this allele
     */
    public double getLikelihoodAssociatedWithReadAndAllele(final GATKSAMRecord read, final Allele allele){
        if (!allelesSet.containsKey(allele) || !likelihoodReadMap.containsKey(read))
            return 0.0;

        return likelihoodReadMap.get(read).get(allele);
    }

    /**
     * Get the most likely alleles estimated across all reads in this object
     *
     * Takes the most likely two alleles according to their diploid genotype likelihoods.  That is, for
     * each allele i and j we compute p(D | i,j) where D is the read likelihoods.  We track the maximum
     * i,j likelihood and return an object that contains the alleles i and j as well as the max likelihood.
     *
     * Note that the second most likely diploid genotype is not tracked so the resulting MostLikelyAllele
     * doesn't have a meaningful get best likelihood.
     *
     * @return a MostLikelyAllele object, or null if this map is empty
     */
    public MostLikelyAllele getMostLikelyDiploidAlleles() {
        if ( isEmpty() ) return null;

        int hap1 = 0;
        int hap2 = 0;
        double maxElement = Double.NEGATIVE_INFINITY;
        for( int iii = 0; iii < alleles.size(); iii++ ) {
            final Allele iii_allele = alleles.get(iii);
            for( int jjj = 0; jjj <= iii; jjj++ ) {
                final Allele jjj_allele = alleles.get(jjj);

                double haplotypeLikelihood = 0.0;
                for( final Map.Entry<GATKSAMRecord, Map<Allele,Double>> entry : likelihoodReadMap.entrySet() ) {
                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    final double likelihood_iii = entry.getValue().get(iii_allele);
                    final double likelihood_jjj = entry.getValue().get(jjj_allele);
                    haplotypeLikelihood += MathUtils.approximateLog10SumLog10(likelihood_iii, likelihood_jjj) + MathUtils.LOG_ONE_HALF;

                    // fast exit.  If this diploid pair is already worse than the max, just stop and look at the next pair
                    if ( haplotypeLikelihood < maxElement ) break;
                }

                // keep track of the max element and associated indices
                if ( haplotypeLikelihood > maxElement ) {
                    hap1 = iii;
                    hap2 = jjj;
                    maxElement = haplotypeLikelihood;
                }
            }
        }

        if ( maxElement == Double.NEGATIVE_INFINITY )
            throw new IllegalStateException("max likelihood is " + maxElement + " indicating something has gone wrong");

        return new MostLikelyAllele(alleles.get(hap1), alleles.get(hap2), maxElement, maxElement);
    }

    /**
     * Given a map from alleles to likelihoods, find the allele with the largest likelihood.
     *
     * @param alleleMap - a map from alleles to likelihoods
     * @return - a MostLikelyAllele object
     */
    @Ensures("result != null")
    public static MostLikelyAllele getMostLikelyAllele( final Map<Allele,Double> alleleMap ) {
        return getMostLikelyAllele(alleleMap, null);
    }

    /**
     * Given a map from alleles to likelihoods, find the allele with the largest likelihood.
     *
     * @param alleleMap - a map from alleles to likelihoods
     * @param onlyConsiderTheseAlleles if not null, we will only consider alleles in this set for being one of the best.
     *                                 this is useful for the case where you've selected a subset of the alleles that
     *                                 the reads have been computed for further analysis.  If null totally ignored
     * @return - a MostLikelyAllele object
     */
    public static MostLikelyAllele getMostLikelyAllele( final Map<Allele,Double> alleleMap, final Set<Allele> onlyConsiderTheseAlleles ) {
        if ( alleleMap == null ) throw new IllegalArgumentException("The allele to likelihood map cannot be null");
        double maxLike = Double.NEGATIVE_INFINITY;
        double prevMaxLike = Double.NEGATIVE_INFINITY;
        Allele mostLikelyAllele = Allele.NO_CALL;
        Allele secondMostLikely = null;

        for (final Map.Entry<Allele,Double> el : alleleMap.entrySet()) {
            if ( onlyConsiderTheseAlleles != null && ! onlyConsiderTheseAlleles.contains(el.getKey()) )
                continue;

            if (el.getValue() > maxLike) {
                prevMaxLike = maxLike;
                maxLike = el.getValue();
                secondMostLikely = mostLikelyAllele;
                mostLikelyAllele = el.getKey();
            } else if( el.getValue() > prevMaxLike ) {
                secondMostLikely = el.getKey();
                prevMaxLike = el.getValue();
            }
        }

        return new MostLikelyAllele(mostLikelyAllele, secondMostLikely, maxLike, prevMaxLike);
    }

    /**
     * Debug method to dump contents of object into string for display
     */
    public String toString() {
        final StringBuilder sb = new StringBuilder();

        sb.append("Alelles in map:");
        for (final Allele a:alleles) {
            sb.append(a.getDisplayString()+",");
        }
        sb.append("\n");
        for (final Map.Entry <GATKSAMRecord, Map<Allele, Double>> el : getLikelihoodReadMap().entrySet() ) {
            for (final Map.Entry<Allele,Double> eli : el.getValue().entrySet()) {
                sb.append("Read "+el.getKey().getReadName()+". Allele:"+eli.getKey().getDisplayString()+" has likelihood="+Double.toString(eli.getValue())+"\n");
            }

        }
        return sb.toString();
    }

    /**
     * Remove reads from this map that are poorly modelled w.r.t. their per allele likelihoods
     *
     * Goes through each read in this map, and if it is poorly modelled removes it from the map.
     *
     * @see #readIsPoorlyModelled(org.broadinstitute.gatk.utils.sam.GATKSAMRecord, java.util.Collection, double)
     * for more information about the poorly modelled test.
     *
     * @param maxErrorRatePerBase see equivalent parameter in #readIsPoorlyModelled
     * @return the list of reads removed from this map because they are poorly modelled
     */
    public List<GATKSAMRecord> filterPoorlyModelledReads(final double maxErrorRatePerBase) {
        final List<GATKSAMRecord> removedReads = new LinkedList<>();
        final Iterator<Map.Entry<GATKSAMRecord, Map<Allele, Double>>> it = likelihoodReadMap.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<GATKSAMRecord, Map<Allele, Double>> record = it.next();
            if ( readIsPoorlyModelled(record.getKey(), record.getValue().values(), maxErrorRatePerBase) ) {
                it.remove();
                removedReads.add(record.getKey());
            }
        }

        return removedReads;
    }

    /**
     * Is this read poorly modelled by all of the alleles in this map?
     *
     * A read is poorly modeled when it's likelihood is below what would be expected for a read
     * originating from one of the alleles given the maxErrorRatePerBase of the reads in general.
     *
     * This function makes a number of key assumptions.  First, that the likelihoods reflect the total likelihood
     * of the read.  In other words, that the read would be fully explained by one of the alleles.  This means
     * that the allele should be something like the full haplotype from which the read might originate.
     *
     * It further assumes that each error in the read occurs with likelihood of -3 (Q30 confidence per base).  So
     * a read with a 10% error rate with Q30 bases that's 100 bp long we'd expect to see 10 real Q30 errors
     * even against the true haplotype.  So for this read to be well modelled by at least one allele we'd expect
     * a likelihood to be >= 10 * -3.
     *
     * @param read the read we want to evaluate
     * @param log10Likelihoods a list of the log10 likelihoods of the read against a set of haplotypes.
     * @param maxErrorRatePerBase the maximum error rate we'd expect for this read per base, in real space.  So
     *                            0.01 means a 1% error rate
     * @return true if none of the log10 likelihoods imply that the read truly originated from one of the haplotypes
     */
    protected boolean readIsPoorlyModelled(final GATKSAMRecord read, final Collection<Double> log10Likelihoods, final double maxErrorRatePerBase) {
        final double maxErrorsForRead = Math.min(2.0, Math.ceil(read.getReadLength() * maxErrorRatePerBase));
        final double log10QualPerBase = -4.0;
        final double log10MaxLikelihoodForTrueAllele = maxErrorsForRead * log10QualPerBase;

        for ( final double log10Likelihood : log10Likelihoods )
            if ( log10Likelihood >= log10MaxLikelihoodForTrueAllele )
                return false;

        return true;
    }

    /**
     * Get an unmodifiable set of the unique alleles in this PerReadAlleleLikelihoodMap
     * @return a non-null unmodifiable map
     */
    public Set<Allele> getAllelesSet() {
        return Collections.unmodifiableSet(allelesSet.keySet());
    }

    /**
     * Loop over all of the reads in this likelihood map and realign them to its most likely haplotype
     * @param haplotypes            the collection of haplotypes
     * @param paddedReferenceLoc    the active region
     */
    public void realignReadsToMostLikelyHaplotype(final Collection<Haplotype> haplotypes, final GenomeLoc paddedReferenceLoc) {

        // we need to remap the Alleles back to the Haplotypes; inefficient but unfortunately this is a requirement currently
        final Map<Allele, Haplotype> alleleToHaplotypeMap = new HashMap<>(haplotypes.size());
        Haplotype refHaplotype = null;
        for ( final Haplotype haplotype : haplotypes ) {
            alleleToHaplotypeMap.put(Allele.create(haplotype.getBases()), haplotype);
            if (refHaplotype == null && haplotype.isReference())
                refHaplotype = haplotype;
        }

        final Map<GATKSAMRecord, Map<Allele, Double>> newLikelihoodReadMap = new LinkedHashMap<>(likelihoodReadMap.size());
        for( final Map.Entry<GATKSAMRecord, Map<Allele, Double>> entry : likelihoodReadMap.entrySet() ) {
            final MostLikelyAllele bestAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(entry.getValue());
            final GATKSAMRecord alignedToRef = AlignmentUtils.createReadAlignedToRef(entry.getKey(), alleleToHaplotypeMap.get(bestAllele.getMostLikelyAllele()), refHaplotype,  paddedReferenceLoc.getStart(), bestAllele.isInformative());
            newLikelihoodReadMap.put(alignedToRef, entry.getValue());
        }

        likelihoodReadMap.clear();
        likelihoodReadMap.putAll(newLikelihoodReadMap);
    }
}
