package org.broadinstitute.sting.oneoffprojects.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.oneoffprojects.refdata.HapmapVCFROD;

import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Arrays;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @author chartl
 * @date Feb 1, 2010
 */
public class ThousandGenomesAnnotator implements InfoFieldAnnotation {

    public List<String> getKeyNames() {
        return Arrays.asList("1KG");
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0),
                1,VCFInfoHeaderLine.INFO_TYPE.String,"Is this site seen in Pilot1 or Pilot2 of 1KG?"));
    }

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, VariantContext vc) {
        if ( tracker == null )
            return null;

        HapmapVCFROD pilot1 = tracker.lookup("pilot1",HapmapVCFROD.class);
        HapmapVCFROD pilot2 = tracker.lookup("pilot2",HapmapVCFROD.class);

        String result;

        if ( pilot1 == null && pilot2 == null) {
            result = "0";
        } else {
            if ( pilot1 != null && ! pilot1.getRecord().isFiltered() ) {
                result = "1";
            } else if ( pilot2 != null && !pilot2.getRecord().isFiltered() ) {
                result = "1";
            } else {
                result = "0";
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), result);
        return map;
    }
}

