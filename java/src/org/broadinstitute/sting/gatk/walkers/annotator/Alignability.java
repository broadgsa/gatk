package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.TabularROD;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.HashMap;
import java.util.Map;


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
        map.put(getKeyName(), String.format("%d", value));
        return map;
    }

    public String getKeyName() { return "Alignability"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(getKeyName(), 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Alignability according to a mask file (3 is best, 0 is worst)"); }
}
