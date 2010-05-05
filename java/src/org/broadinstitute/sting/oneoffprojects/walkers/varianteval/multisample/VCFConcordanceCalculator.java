package org.broadinstitute.sting.oneoffprojects.walkers.varianteval.multisample;

import org.broad.tribble.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.GenomeLoc;

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
    private int falsePositiveLoci;
    private int falseNegativeLoci;
    private int falseNegativeLociDueToNoCall;
    private int falseNegativeLociDueToFilters;
    private int hetsCalledHoms;
    private int homsCalledHets;
    private int nonConfidentGenotypeCalls;
    private int concordantHomCalls;
    private int concordantHetCalls;
    private int concordantGenotypeReferenceCalls;
    private int chipNoCalls;
    private int ignoredDueToDepth;

    public VCFConcordanceCalculator(String sampleName, int minimumDepth, int minGenQual) {
        name = sampleName;
        falseNegativeLoci = 0;
        falseNegativeLociDueToNoCall = 0;
        falsePositiveLoci = 0;
        falseNegativeLociDueToFilters = 0;
        hetsCalledHoms = 0;
        homsCalledHets = 0;
        nonConfidentGenotypeCalls = 0;
        concordantHomCalls = 0;
        concordantHetCalls = 0;
        concordantGenotypeReferenceCalls = 0;
        chipNoCalls = 0;
        ignoredDueToDepth = 0;
        minimumDepthForUpdate = minimumDepth;
        minimumGenotypeQuality = minGenQual;
    }

    public void update(LocusConcordanceInfo info) {
        compareGenotypes(info.getTruthGenotype(name), info.getVariantGenotype(name), info.getLoc(), info.getReferenceBase() );
    }

    public void updateTruthOnly(LocusConcordanceInfo info) {
        if ( info.getTruthGenotype(name).isVariant( (char) info.getReferenceBase() ) ) {
            falseNegativeLoci++;
        } else {
            concordantGenotypeReferenceCalls++;
        }
    }

    public void updateFilteredLocus(LocusConcordanceInfo info) {

        if ( info.getTruthGenotype(name).isVariant( (char) info.getReferenceBase()) ) {
            falseNegativeLociDueToFilters++;
        }
    }


    public String toString() {
        return String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",name,ignoredDueToDepth,
                concordantGenotypeReferenceCalls,concordantHomCalls,concordantHetCalls,nonConfidentGenotypeCalls,
                homsCalledHets,hetsCalledHoms,falsePositiveLoci,falseNegativeLoci,
                falseNegativeLociDueToNoCall,falseNegativeLociDueToFilters);
    }

    private void compareGenotypes(VCFGenotypeRecord truth, VCFGenotypeRecord call, GenomeLoc loc, byte ref) {
        if ( minimumDepthForUpdate > 0 && call.getReadCount() < minimumDepthForUpdate ) {
            ignoredDueToDepth++;
        } else if ( truth.isNoCall() ) {
            chipNoCalls++;
        } else if ( truth.isVariant(( char) ref) ) {
            if ( call.isNoCall() ) {
                falseNegativeLociDueToNoCall++;
            } else if ( ! call.isVariant( (char) ref ) ) {
                falseNegativeLoci++;
            } else if ( call.isVariant((char) ref) ) {
                // check het vs hom
                checkGenotypeCall(truth,call, loc);
            }

        } else if ( ! truth.isVariant( (char) ref ) ) {

            if ( call.isVariant((char) ref) ) {
                falsePositiveLoci++;
            } else {
                concordantGenotypeReferenceCalls++;
            }
        }
    }

    private void checkGenotypeCall( VCFGenotypeRecord truth, VCFGenotypeRecord call, GenomeLoc loc ) {
        if ( ! call.isFiltered() && 10*call.getNegLog10PError() > minimumGenotypeQuality) {

            if ( truth.isHet() && call.isHom() ) {
                hetsCalledHoms++;
            } else if ( truth.isHom() && call.isHet() ) {
                homsCalledHets++;
            } else if (  ( truth.isHet() && call.isHet()  ) ) {
                concordantHetCalls++;
            } else if ( truth.isHom() && call.isHom() ) { // be extra careful
                concordantHomCalls++;
            }

        } else {
            
            nonConfidentGenotypeCalls++;
        }

    }
}
