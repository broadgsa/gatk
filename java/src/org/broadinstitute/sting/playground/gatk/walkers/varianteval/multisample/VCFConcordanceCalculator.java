package org.broadinstitute.sting.playground.gatk.walkers.varianteval.multisample;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;

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
    private String name;
    private Set<GenomeLoc> falsePositiveLoci;
    private Set<GenomeLoc> falseNegativeLoci;
    private Set<GenomeLoc> falseNegativeLociDueToNoCall;
    private Set<GenomeLoc> hetsCalledHoms;
    private Set<GenomeLoc> homsCalledHets;
    private Set<GenomeLoc> concordantCalls;
    private Set<GenomeLoc> concordantGenotypeReferenceCalls;
    private Set<GenomeLoc> chipNoCalls;

    public VCFConcordanceCalculator(String sampleName) {
        name = sampleName;
        falseNegativeLoci = new HashSet<GenomeLoc>();
        falseNegativeLociDueToNoCall = new HashSet<GenomeLoc>();
        falsePositiveLoci = new HashSet<GenomeLoc>();
        hetsCalledHoms = new HashSet<GenomeLoc>();
        homsCalledHets = new HashSet<GenomeLoc>();
        concordantCalls = new HashSet<GenomeLoc>();
        concordantGenotypeReferenceCalls = new HashSet<GenomeLoc>();
        chipNoCalls = new HashSet<GenomeLoc>();
    }

    public void update(LocusConcordanceInfo info) {
        compareGenotypes(info.getTruthGenotype(name), info.getVariantGenotype(name), info.getLoc(), info.getReferenceBase() );
    }

    public String toString() {
        return String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d",name,concordantGenotypeReferenceCalls.size(),concordantCalls.size(),homsCalledHets.size(),hetsCalledHoms.size(),falsePositiveLoci.size(),falseNegativeLoci.size(),falseNegativeLociDueToNoCall.size());
    }

    private void compareGenotypes(VCFGenotypeRecord truth, VCFGenotypeRecord call, GenomeLoc loc, byte ref) {
        if ( truth.isNoCall() ) {
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
        if ( truth.isHet() && call.isHom() ) {
            hetsCalledHoms.add(loc);
        } else if ( truth.isHom() && call.isHet() ) {
            homsCalledHets.add(loc);
        } else if (  ( truth.isHet() && call.isHet() ) || ( truth.isHom() && call.isHom() ) ) { // be extra careful
            concordantCalls.add(loc);
        }
    }
}
