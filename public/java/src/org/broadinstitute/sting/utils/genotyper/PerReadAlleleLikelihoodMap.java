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
     * Adds a read, allele and corresponding likelihood to map
     * @param read                       SAM record to add
     * @param a                          corresponding allele
     * @param likelihood                 corresponding likelihood
     */
    public void add(GATKSAMRecord read, Allele a, Double likelihood) {
        Map<Allele,Double> likelihoodMap;
        if (likelihoodReadMap.containsKey(read)){
            // seen pileup element before
            likelihoodMap = likelihoodReadMap.get(read);
        }
        else {
            likelihoodMap = new HashMap<Allele, Double>();
            likelihoodReadMap.put(read,likelihoodMap);
        }
        likelihoodMap.put(a,likelihood);

        if (!alleles.contains(a))
            alleles.add(a);

    }

    public ReadBackedPileup createPerAlleleDownsampledBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction, final PrintStream log) {
        return AlleleBiasedDownsamplingUtils.createAlleleBiasedBasePileup(pileup, downsamplingFraction, log);
    }

    public void performPerAlleleDownsampling(final double downsamplingFraction, final PrintStream log) {
        // special case removal of all or no reads
        if ( downsamplingFraction <= 0.0 )
            return;
        if ( downsamplingFraction >= 1.0 ) {
            likelihoodReadMap.clear();
            return;
        }

        // start by stratifying the reads by the alleles they represent at this position
        final Map<Allele, List<GATKSAMRecord>> alleleReadMap = new HashMap<Allele, List<GATKSAMRecord>>(alleles.size());
        for ( Allele allele : alleles )
            alleleReadMap.put(allele, new ArrayList<GATKSAMRecord>());

        for ( Map.Entry<GATKSAMRecord, Map<Allele, Double>> entry : likelihoodReadMap.entrySet() ) {
            // do not remove reduced reads!
            if ( !entry.getKey().isReducedRead() ) {
                final Allele bestAllele = getMostLikelyAllele(entry.getValue());
                if ( bestAllele != Allele.NO_CALL )
                    alleleReadMap.get(bestAllele).add(entry.getKey());
            }
        }

        // compute the reads to remove and actually remove them
        final List<GATKSAMRecord> readsToRemove = AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(alleleReadMap, downsamplingFraction, log);
        for ( final GATKSAMRecord read : readsToRemove )
            likelihoodReadMap.remove(read);
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
        if (p==null || p.getRead()==null || a == null )
            throw new IllegalArgumentException("Invalid parameters passed to PerReadAlleleLikelihoodMap.add");
        add(p.getRead(), a, likelihood);
    }

    /**
     * Does the current map contain the key associated with a particular SAM record in pileup?
     * @param p                 Pileup element
     * @return
     */
    public boolean containsPileupElement(PileupElement p) {
        if (p==null )
           throw new IllegalArgumentException("Invalid pileup element");

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

    public Map<Allele,Double> getLikelihoodsAssociatedWithPileupElement(PileupElement p) {
        if (!likelihoodReadMap.containsKey(p.getRead()))
            return null;

        return likelihoodReadMap.get(p.getRead());
    }


    /**
     * For a given alleleMap, return most likely allele, i.e. the one with highest associated likelihood
     * @param alleleMap          Underlying allele map
     * @return                   Most likely allele. If all alleles are equally likely, returns a no-call allele.
     */
    @Ensures("result != null")
    public static Allele getMostLikelyAllele( final Map<Allele,Double> alleleMap ) {
        double maxLike = Double.NEGATIVE_INFINITY;
        double prevMaxLike = Double.NEGATIVE_INFINITY;
        Allele mostLikelyAllele = Allele.NO_CALL;

        if (alleleMap==null)
            throw new IllegalArgumentException("alleleMap in getMostLikelyAllele() method can't be null");

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
}
