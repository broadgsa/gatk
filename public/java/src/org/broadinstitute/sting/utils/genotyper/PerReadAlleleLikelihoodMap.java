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
package org.broadinstitute.sting.utils.genotyper;


import org.broadinstitute.sting.utils.classloader.GATKLiteUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.*;

public abstract class PerReadAlleleLikelihoodMap {

    public static final double INFORMATIVE_LIKELIHOOD_THRESHOLD = 0.1;

    protected List<Allele> alleles;
    protected Map<GATKSAMRecord, Map<Allele, Double>> likelihoodReadMap;

    public abstract void performPerAlleleDownsampling(final double downsamplingFraction, final PrintStream log);
    public abstract ReadBackedPileup createPerAlleleDownsampledBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction, final PrintStream log);

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

    public int size() {
        return likelihoodReadMap.size();
    }

    public void add(PileupElement p, Allele a, Double likelihood) {
        add(p.getRead(), a, likelihood);
    }

    public boolean containsPileupElement(PileupElement p) {
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

    public static Allele getMostLikelyAllele( final Map<Allele,Double> alleleMap ) {
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

    public static PerReadAlleleLikelihoodMap getBestAvailablePerReadAlleleLikelihoodMap() {
        final Class PerReadAlleleLikelihoodMapClass = GATKLiteUtils.getProtectedClassIfAvailable(PerReadAlleleLikelihoodMap.class);
        try {
            Constructor constructor = PerReadAlleleLikelihoodMapClass.getDeclaredConstructor((Class[])null);
            constructor.setAccessible(true);
            return (PerReadAlleleLikelihoodMap)constructor.newInstance();
        }
        catch (Exception e) {
            throw new ReviewedStingException("Unable to create RecalibrationEngine class instance " + PerReadAlleleLikelihoodMapClass.getSimpleName());
        }
    }
}
