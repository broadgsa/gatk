package org.broadinstitute.sting.oneoffprojects.walkers.varianteval.multisample;

import org.broad.tribble.vcf.VCFGenotypeRecord;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 27, 2010
 * Time: 5:48:36 PM
 * To change this template use File | Settings | File Templates.
 */
class LocusConcordanceInfo {

    public enum ConcordanceType {
        TRUTH_SET,TRUTH_SET_VARIANT_FILTERED,VARIANT_SET,BOTH_SETS
    }

    private ConcordanceType concordanceType;
    private VCFRecord variantVCFRecord;
    private VCFRecord truthVCFRecord;
    private ReferenceContext reference;

    public LocusConcordanceInfo(ConcordanceType type, VCFRecord truthRecord, VCFRecord variantRecord, ReferenceContext ref) {
        concordanceType = type;
        variantVCFRecord = variantRecord;
        truthVCFRecord = truthRecord;
        reference = ref;
    }

    public boolean concordanceIsCheckable() {
        return concordanceType == ConcordanceType.BOTH_SETS;
    }

    public VCFGenotypeRecord getTruthGenotype(String sample) {
        return truthVCFRecord.getGenotype(sample);
    }

    public VCFGenotypeRecord getVariantGenotype(String sample) {
        return variantVCFRecord.getGenotype(sample);
    }

    public Set<String> getOverlappingSamples() {
        Set<String> variantSamples = new HashSet<String>( Arrays.asList(variantVCFRecord.getSampleNames()) );
        variantSamples.retainAll(Arrays.asList(truthVCFRecord.getSampleNames()));
        return variantSamples;
    }

    public byte getReferenceBase() {
        return reference.getBase();
    }

    public boolean isTruthOnly () {
        return concordanceType == ConcordanceType.TRUTH_SET;
    }

    public boolean isVariantSite() {
        for ( VCFGenotypeRecord g : truthVCFRecord.getVCFGenotypeRecords() ) {
            if ( g.isVariant(reference.getBaseAsChar()) ) {
                return true;
            }
        }

        return false;
    }

    public boolean isVariantFiltered() {
        return this.concordanceType == ConcordanceType.TRUTH_SET_VARIANT_FILTERED;
    }

    public GenomeLoc getLoc() {
        if ( concordanceType == ConcordanceType.TRUTH_SET || concordanceType == ConcordanceType.BOTH_SETS || concordanceType == ConcordanceType.TRUTH_SET_VARIANT_FILTERED) {
            return GenomeLocParser.createGenomeLoc(truthVCFRecord.getChr(),truthVCFRecord.getStart());
        } else {
            return GenomeLocParser.createGenomeLoc( variantVCFRecord.getChr(),variantVCFRecord.getStart());
        }
    }

}
