/*
 * Copyright (c) 2011 The Broad Institute
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
package org.broadinstitute.sting.gatk.walkers.genotyper;


//import org.broadinstitute.sting.gatk.walkers.Requires;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;

public class PerReadAlleleLikelihoodMap {
    public static final double INDEL_LIKELIHOOD_THRESH = 0.1;

    private List<Allele> alleles;
    private Map<PileupElement,Map<Allele,Double>> likelihoodReadMap;
    public PerReadAlleleLikelihoodMap() {
        likelihoodReadMap = new LinkedHashMap<PileupElement,Map<Allele,Double>>();
        alleles = new ArrayList<Allele>();
    }

    public void add(PileupElement p, Allele a, Double likelihood) {
        Map<Allele,Double> likelihoodMap;
        if (likelihoodReadMap.containsKey(p)){
            // seen pileup element before
            likelihoodMap = likelihoodReadMap.get(p);
        }
        else {
            likelihoodMap = new HashMap<Allele, Double>();
            likelihoodReadMap.put(p,likelihoodMap);
        }
        likelihoodMap.put(a,likelihood);

        if (!alleles.contains(a))
            alleles.add(a);

    }

    public int size() {
        return likelihoodReadMap.size();
    }

    public void add(GATKSAMRecord read, Allele a, Double likelihood) {
        PileupElement p = new PileupElement(read,-1,false,false,false,false,false,false);
        add(p,a,likelihood);
    }

    public boolean containsPileupElement(PileupElement p) {
        return likelihoodReadMap.containsKey(p);
    }

    public boolean isEmpty() {
        return likelihoodReadMap.isEmpty();
    }

    public Map<PileupElement,Map<Allele,Double>> getLikelihoodReadMap() {
        return likelihoodReadMap;
    }
    public void clear() {
        alleles.clear();
        likelihoodReadMap.clear();
    }

    public Set<PileupElement> getStoredPileupElements() {
        return likelihoodReadMap.keySet();
    }
    /**
     * Returns list of reads greedily associated with a particular allele.
     * Needs to loop for each read, and assign to each allele
     * @param a                          Desired allele
     * @return
     */
    @Requires("a!=null")
    public List<GATKSAMRecord> getReadsAssociatedWithAllele(Allele a) {
        return null;
    }

    public Map<Allele,Double> getLikelihoodsAssociatedWithPileupElement(PileupElement p) {
        if (!likelihoodReadMap.containsKey(p))
            return null;

        return likelihoodReadMap.get(p);
    }

    public static Allele getMostLikelyAllele(Map<Allele,Double> alleleMap) {
        double minLike = Double.POSITIVE_INFINITY, maxLike = Double.NEGATIVE_INFINITY;
        Allele mostLikelyAllele = Allele.NO_CALL;

        for (Map.Entry<Allele,Double> el : alleleMap.entrySet()) {
            if (el.getValue() > maxLike) {
                maxLike = el.getValue();
                mostLikelyAllele = el.getKey();
            }

            if (el.getValue() < minLike)
                minLike = el.getValue();

        }
        if (maxLike-minLike > INDEL_LIKELIHOOD_THRESH)
            return mostLikelyAllele;
        else
            return Allele.NO_CALL;
    }
 }
