package org.broadinstitute.sting.oneoffprojects.walkers.varianteval.multisample;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeEncoding;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 27, 2010
 * Time: 5:48:08 PM
 * To change this template use File | Settings | File Templates.
 */
class VCFConcordanceCalculator {

    private int minimumDepthForUpdate;
    private int minimumGenotypeQuality;
    private String name;
    private Set<GenomeLoc> falsePositiveLoci;
    private Set<GenomeLoc> falseNegativeLoci;
    private Set<GenomeLoc> falseNegativeLociDueToNoCall;
    private Set<GenomeLoc> falseNegativeLociDueToFilters;
    private Set<GenomeLoc> hetsCalledHoms;
    private Set<GenomeLoc> homsCalledHets;
    private Set<GenomeLoc> nonConfidentGenotypeCalls;
    private Set<GenomeLoc> concordantHomCalls;
    private Set<GenomeLoc> concordantHetCalls;
    private Set<GenomeLoc> concordantGenotypeReferenceCalls;
    private Set<GenomeLoc> chipNoCalls;
    private Set<GenomeLoc> ignoredDueToDepth;

    public VCFConcordanceCalculator(String sampleName, int minimumDepth, int minGenQual) {
        name = sampleName;
        falseNegativeLoci = new HashSet<GenomeLoc>();
        falseNegativeLociDueToNoCall = new HashSet<GenomeLoc>();
        falsePositiveLoci = new HashSet<GenomeLoc>();
        falseNegativeLociDueToFilters = new HashSet<GenomeLoc>();
        hetsCalledHoms = new HashSet<GenomeLoc>();
        homsCalledHets = new HashSet<GenomeLoc>();
        nonConfidentGenotypeCalls = new HashSet<GenomeLoc>();
        concordantHomCalls = new HashSet<GenomeLoc>();
        concordantHetCalls = new HashSet<GenomeLoc>();
        concordantGenotypeReferenceCalls = new HashSet<GenomeLoc>();
        chipNoCalls = new HashSet<GenomeLoc>();
        ignoredDueToDepth = new HashSet<GenomeLoc>();
        minimumDepthForUpdate = minimumDepth;
        minimumGenotypeQuality = minGenQual;
    }

    public void update(LocusConcordanceInfo info) {
        compareGenotypes(info.getTruthGenotype(name), info.getVariantGenotype(name), info.getLoc(), info.getReferenceBase() );
    }

    public void updateTruthOnly(LocusConcordanceInfo info) {
        if ( info.getTruthGenotype(name).isVariant( (char) info.getReferenceBase() ) ) {
            falseNegativeLoci.add(info.getLoc());
        } else {
            concordantGenotypeReferenceCalls.add(info.getLoc());
        }
    }

    public void updateFilteredLocus(LocusConcordanceInfo info) {

        if ( info.getTruthGenotype(name).isVariant( (char) info.getReferenceBase()) ) {
            falseNegativeLociDueToFilters.add(info.getLoc());
        }
    }


    public String toString() {
        return String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",name,ignoredDueToDepth.size(),
                concordantGenotypeReferenceCalls.size(),concordantHomCalls.size(),concordantHetCalls.size(),nonConfidentGenotypeCalls.size(),
                homsCalledHets.size(),hetsCalledHoms.size(),falsePositiveLoci.size(),falseNegativeLoci.size(),
                falseNegativeLociDueToNoCall.size(),falseNegativeLociDueToFilters.size());
    }

    private void compareGenotypes(VCFGenotypeRecord truth, VCFGenotypeRecord call, GenomeLoc loc, byte ref) {
        if ( minimumDepthForUpdate > 0 && call.getReadCount() < minimumDepthForUpdate ) {
            ignoredDueToDepth.add(loc);
        } else if ( truth.isNoCall() ) {
            chipNoCalls.add(loc);
        } else if ( truth.isVariant(( char) ref) ) {
            if ( call.isNoCall() ) {
                falseNegativeLociDueToNoCall.add(loc);
            } else if ( ! call.isVariant( (char) ref ) ) {
                falseNegativeLoci.add(loc);
            } else if ( call.isVariant((char) ref) ) {
                // check het vs hom
                checkGenotypeCall(truth,call, loc);
            }

        } else if ( ! truth.isVariant( (char) ref ) ) {

            if ( call.isVariant((char) ref) ) {
                falsePositiveLoci.add(loc);
            } else {
                concordantGenotypeReferenceCalls.add(loc);
            }
        }
    }

    private void checkGenotypeCall( VCFGenotypeRecord truth, VCFGenotypeRecord call, GenomeLoc loc ) {
        if ( ! call.isFiltered() && 10*call.getNegLog10PError() > minimumGenotypeQuality) {

            if ( truth.isHet() && call.isHom() ) {
                hetsCalledHoms.add(loc);
            } else if ( truth.isHom() && call.isHet() ) {
                homsCalledHets.add(loc);
            } else if (  ( truth.isHet() && call.isHet()  ) ) {
                concordantHetCalls.add(loc);
            } else if ( truth.isHom() && call.isHom() ) { // be extra careful
                concordantHomCalls.add(loc);
            }

        } else {
            
            nonConfidentGenotypeCalls.add(loc);
        }

    }
}
