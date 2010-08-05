package org.broadinstitute.sting.oneoffprojects.walkers.annotator;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Arrays;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Mar 29, 2010
 */
public class InsertSizeDistribution implements InfoFieldAnnotation {
    private final long INSERT_SIZE_LOWER_BOUND = 500;
    public List<String> getKeyNames() { return Arrays.asList("INSIZE"); }
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0),1, VCFHeaderLineType.Integer,"Do not use this if your name is not Chris")); }

    public Map<String,Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, VariantContext variant) {
        int weirdInsertSizeReads = 0;
        for ( String sample : context.keySet() ) {
            ReadBackedPileup pileup = context.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for (PileupElement e : pileup ) {
                if ( Math.abs(e.getRead().getInferredInsertSize()) > INSERT_SIZE_LOWER_BOUND ) {
                    weirdInsertSizeReads++;
                }
            }
        }

        Map<String,Object> toReturn = new HashMap<String,Object>();
        toReturn.put(getKeyNames().get(0),String.format("%d",weirdInsertSizeReads));
        return toReturn;
    }
}
