package org.broadinstitute.sting.oneoffprojects.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotation;
import org.broadinstitute.sting.oneoffprojects.refdata.HapmapVCFROD;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @author chartl
 * @date Feb 1, 2010
 */
public class ThousandGenomesAnnotator implements VariantAnnotation {

    public String getKeyName() {
        return "1KG";
    }

    public VCFInfoHeaderLine getDescription() {
        return new VCFInfoHeaderLine(getKeyName(),
                1,VCFInfoHeaderLine.INFO_TYPE.String,"Is this site seen in Pilot1 or Pilot2 of 1KG?");
    }

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, VariantContext vc) {
        if ( tracker == null ) {
            return null;
        }

        RODRecordList pilot1 = tracker.getTrackData("pilot1",null);
        RODRecordList pilot2 = tracker.getTrackData("pilot2",null);

        if ( pilot1 == null && pilot2 == null) {
            return "0";
        } else {
            if ( pilot1 != null && ! ( (HapmapVCFROD) pilot1.get(0)).getRecord().isFiltered() ) {
                return "1";
            } else if ( pilot2 != null && ! ( (HapmapVCFROD) pilot2.get(0)).getRecord().isFiltered() ) {
                return "1";
            } else {
                return "0";
            }
        }
    }
}

