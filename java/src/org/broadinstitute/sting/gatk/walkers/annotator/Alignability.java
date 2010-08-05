package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.TabularROD;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;

import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Arrays;


public class Alignability implements InfoFieldAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker,
						   ReferenceContext ref, 
						   Map<String, StratifiedAlignmentContext> stratifiedContexts, 
						   VariantContext vc)
	{
		TabularROD record = tracker.lookup("alignability",TabularROD.class);
        if (record == null)
            return null;

        if (record.get("alignability") == null)
            throw new RuntimeException("ERROR: alignability column not defined in alignability input.\n");

        int value = Integer.parseInt(record.get("alignability"));

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", value));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("Alignability"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Integer, "Alignability according to a mask file (3 is best, 0 is worst)")); }
}
