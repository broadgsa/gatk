package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Total (unfiltered) depth over all samples.
 *
 * This and AD are complementary fields that are two important ways of thinking about the depth of the data for this sample
 * at this site.  The DP field describe the total depth of reads that passed the Unified Genotypers internal
 * quality control metrics (like MAPQ > 17, for example), whatever base was present in the read at this site.
 * The AD values (one for each of REF and ALT fields) is the count of all reads that carried with them the
 * REF and ALT alleles. The reason for this distinction is that the DP is in some sense reflective of the
 * power I have to determine the genotype of the sample at this site, while the AD tells me how many times
 * I saw each of the REF and ALT alleles in the reads, free of any bias potentially introduced by filtering
 * the reads. If, for example, I believe there really is a an A/T polymorphism at a site, then I would like
 * to know the counts of A and T bases in this sample, even for reads with poor mapping quality that would
 * normally be excluded from the statistical calculations going into GQ and QUAL.
 *
 * Note that the DP is affected by downsampling (-dcov) though, so the max value one can obtain for N samples with
 * -dcov D is N * D
 */
public class DepthOfCoverage extends InfoFieldAnnotation implements StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int depth = 0;
        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() )
            depth += sample.getValue().hasBasePileup() ? sample.getValue().getBasePileup().depthOfCoverage() : sample.getValue().getExtendedEventPileup().depthOfCoverage();
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", depth));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.DEPTH_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Integer, "Filtered Depth")); }
}
