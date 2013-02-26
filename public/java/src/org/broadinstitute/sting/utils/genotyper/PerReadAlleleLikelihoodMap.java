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

package org.broadinstitute.sting.utils.genotyper;


import com.google.java.contract.Ensures;
import org.broadinstitute.sting.gatk.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.variantcontext.Allele;

import java.io.PrintStream;
import java.util.*;

/**
 *   Wrapper class that holds a set of maps of the form (Read -> Map(Allele->Double))
 *   For each read, this holds underlying alleles represented by an aligned read, and corresponding relative likelihood.
 */
public class PerReadAlleleLikelihoodMap {


    public static final double INFORMATIVE_LIKELIHOOD_THRESHOLD = 0.2;

    protected List<Allele> alleles;
    protected Map<GATKSAMRecord, Map<Allele, Double>> likelihoodReadMap;

    public PerReadAlleleLikelihoodMap() {
        likelihoodReadMap = new LinkedHashMap<GATKSAMRecord,Map<Allele,Double>>();
        alleles = new ArrayList<Allele>();
    }

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
        Map<Allele,Double> likelihoodMap = likelihoodReadMap.get(read);
        if (likelihoodMap == null){
            // LinkedHashMap will ensure iterating through alleles will be in consistent order
            likelihoodMap = new LinkedHashMap<Allele, Double>();
        }
        likelihoodReadMap.put(read,likelihoodMap);

        likelihoodMap.put(a,likelihood);

        if (!alleles.contains(a))
            alleles.add(a);

    }

    public ReadBackedPileup createPerAlleleDownsampledBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction, final PrintStream log) {
        return AlleleBiasedDownsamplingUtils.createAlleleBiasedBasePileup(pileup, downsamplingFraction, log);
    }

    /**
     * For each allele "a" , identify those reads whose most likely allele is "a", and remove a "downsamplingFraction" proportion
     * of those reads from the "likelihoodReadMap". This is used for e.g. sample contamination
     * @param downsamplingFraction - the fraction of supporting reads to remove from each allele. If <=0 all reads kept, if >=1 all reads tossed.
     * @param log - a PrintStream to log the removed reads to (passed through to the utility function)
     */
    public void performPerAlleleDownsampling(final double downsamplingFraction, final PrintStream log) {
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
        final List<GATKSAMRecord> readsToRemove = AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(alleleReadMap, downsamplingFraction, log);
        for ( final GATKSAMRecord read : readsToRemove )
            likelihoodReadMap.remove(read);
    }

    /**
     * Convert the @likelihoodReadMap to a map of alleles to reads, where each read is mapped uniquely to the allele
     * for which it has the greatest associated likelihood
     * @return a map from each allele to a list of reads that 'support' the allele
     */
    protected Map<Allele,List<GATKSAMRecord>> getAlleleStratifiedReadMap() {
        final Map<Allele, List<GATKSAMRecord>> alleleReadMap = new HashMap<Allele, List<GATKSAMRecord>>(alleles.size());
        for ( final Allele allele : alleles )
            alleleReadMap.put(allele, new ArrayList<GATKSAMRecord>());

        for ( final Map.Entry<GATKSAMRecord, Map<Allele, Double>> entry : likelihoodReadMap.entrySet() ) {
            // do not remove reduced reads!
            if ( !entry.getKey().isReducedRead() ) {
                final Allele bestAllele = getMostLikelyAllele(entry.getValue());
                if ( bestAllele != Allele.NO_CALL )
                    alleleReadMap.get(bestAllele).add(entry.getKey());
            }
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
     * @return
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
        alleles.clear();
        likelihoodReadMap.clear();
    }

    public Set<GATKSAMRecord> getStoredElements() {
        return likelihoodReadMap.keySet();
    }

    public Collection<Map<Allele,Double>> getLikelihoodMapValues() {
        return likelihoodReadMap.values();
    }

    public int getNumberOfStoredElements() {
        return likelihoodReadMap.size();
    }

    public Map<Allele,Double> getLikelihoodsAssociatedWithPileupElement(final PileupElement p) {
        if (!likelihoodReadMap.containsKey(p.getRead()))
            return null;

        return likelihoodReadMap.get(p.getRead());
    }


    /**
     * Given a map from alleles to likelihoods, find the allele with the largest likelihood.
     * If the difference between the most-likely allele and the next-most-likely allele is < INFORMATIVE_LIKELIHOOD_THRESHOLD
     * then the most likely allele is set to "no call"
     * @param alleleMap - a map from alleles to likelihoods
     * @return - the most likely allele, or NO_CALL if two or more alleles have likelihoods within INFORMATIVE_LIKELIHOOD_THRESHOLD
     * of one another. By default empty allele maps will return NO_CALL, and allele maps with a single entry will return the
     * corresponding key
     */
    @Ensures("result != null")
    public static Allele getMostLikelyAllele( final Map<Allele,Double> alleleMap ) {
        if ( alleleMap == null ) throw new IllegalArgumentException("The allele to likelihood map cannot be null");
        double maxLike = Double.NEGATIVE_INFINITY;
        double prevMaxLike = Double.NEGATIVE_INFINITY;
        Allele mostLikelyAllele = Allele.NO_CALL;

        for (final Map.Entry<Allele,Double> el : alleleMap.entrySet()) {
            if (el.getValue() > maxLike) {
                prevMaxLike = maxLike;
                maxLike = el.getValue();
                mostLikelyAllele = el.getKey();
            } else if( el.getValue() > prevMaxLike ) {
                prevMaxLike = el.getValue();
            }
        }
        return (maxLike - prevMaxLike > INFORMATIVE_LIKELIHOOD_THRESHOLD ? mostLikelyAllele : Allele.NO_CALL );
    }


    /**
     * Debug method to dump contents of object into string for display
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append("Alelles in map:");
        for (Allele a:alleles) {
            sb.append(a.getDisplayString()+",");

        }
        sb.append("\n");
        for (Map.Entry <GATKSAMRecord, Map<Allele, Double>> el : getLikelihoodReadMap().entrySet() ) {
            for (Map.Entry<Allele,Double> eli : el.getValue().entrySet()) {
                sb.append("Read "+el.getKey().getReadName()+". Allele:"+eli.getKey().getDisplayString()+" has likelihood="+Double.toString(eli.getValue())+"\n");
            }

        }
        return sb.toString();
    }
}
