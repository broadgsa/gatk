package org.broadinstitute.sting.oneoffprojects.walkers.varianteval.multisample;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 27, 2010
 * Time: 5:47:27 PM
 * To change this template use File | Settings | File Templates.
 */
class MultiSampleConcordanceSet {
    private boolean treatTruthOnlyAsFalseNegative;
    private int minimumDepthForTest;
    private HashSet<VCFConcordanceCalculator> concordanceSet;
    private Set<String> cachedSampleNames;
    private long truthOnlySites;
    private long truthOnlyVariantSites;
    private long variantOnlySites;
    private long overlappingSites;
    private int genotypeQuality;

    public MultiSampleConcordanceSet(int minDepth, boolean assumeRef, int genotypeQuality) {
        concordanceSet = new HashSet<VCFConcordanceCalculator>();
        truthOnlySites = 0l;
        truthOnlyVariantSites = 0l;
        variantOnlySites = 0l;
        overlappingSites = 0l;
        minimumDepthForTest = minDepth;
        treatTruthOnlyAsFalseNegative = assumeRef;
        this.genotypeQuality = genotypeQuality;
    }

    public boolean hasBeenInstantiated() {
        return cachedSampleNames != null;
    }

    public void instantiate(Set<String> samples) {
        cachedSampleNames = samples;
        for ( String s : samples ) {
            concordanceSet.add(new VCFConcordanceCalculator(s,minimumDepthForTest,genotypeQuality));
        }
    }

    public void update(LocusConcordanceInfo info) {
        if ( info.concordanceIsCheckable() ) {
            overlappingSites++;
            for ( VCFConcordanceCalculator concordance : concordanceSet ) {
                concordance.update(info);
            }
        } else if ( info.isTruthOnly() ) {
            truthOnlySites++;
            if ( info.isVariantSite() ) {
                truthOnlyVariantSites++;
                if ( treatTruthOnlyAsFalseNegative ) {
                    for ( VCFConcordanceCalculator concordance : concordanceSet ) {
                        concordance.updateTruthOnly(info);
                    }
                }
            }
        } else {
            variantOnlySites++;
        }
    }

    public Set<VCFConcordanceCalculator> getConcordanceSet() {
        return concordanceSet;
    }

    public long numberOfTruthOnlySites() {
        return truthOnlySites;
    }

    public long numberOfTruthOnlyVariantSites() {
        return truthOnlyVariantSites;
    }

    public long numberOfVariantOnlySites() {
        return variantOnlySites;
    }

    public long numberOfOverlappingSites() {
        return overlappingSites;
    }
}
