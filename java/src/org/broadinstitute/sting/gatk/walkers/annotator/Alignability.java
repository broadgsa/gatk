package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.TabularROD;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;


public class Alignability implements VariantAnnotation {

    public String annotate(RefMetaDataTracker tracker, 
						   ReferenceContext ref, 
						   Map<String, StratifiedAlignmentContext> stratifiedContexts, 
						   Variation variation) 
	{
		TabularROD record = (TabularROD)(tracker.lookup("alignability", null));
		int value;
		if (record == null) { value = 3; }
		else 
		{ 
			if (record.get("alignability") == null) 
			{ 
				throw new RuntimeException("ERROR: alignability column not defined in alignability input.\n");
			}
			value = Integer.parseInt(record.get("alignability")); 
		}
		return String.format("%d", value);
    }

    public String getKeyName() { return "Alignability"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(getKeyName(), 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Alignability according to a mask file (3 is best, 0 is worst)"); }
}
