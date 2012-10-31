package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.IndelUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Rough category of indel type (insertion, deletion, multi-allelic, other)
 */
public class IndelType extends InfoFieldAnnotation implements ExperimentalAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {

        int run;
        if (vc.isMixed()) {
            Map<String, Object> map = new HashMap<String, Object>();
            map.put(getKeyNames().get(0), String.format("%s", "MIXED"));
            return map;

        }
        else if ( vc.isIndel() ) {
            String type="";
            if (!vc.isBiallelic())
                type = "MULTIALLELIC_INDEL";
            else {
                if (vc.isSimpleInsertion())
                    type = "INS.";
                else if (vc.isSimpleDeletion())
                    type = "DEL.";
                else
                    type = "OTHER.";
                ArrayList<Integer> inds = IndelUtils.findEventClassificationIndex(vc, ref);
                for (int k : inds) {
                    type = type+ IndelUtils.getIndelClassificationName(k)+".";
                }
            }
            Map<String, Object> map = new HashMap<String, Object>();
            map.put(getKeyNames().get(0), String.format("%s", type));
            return map;

        } else {
            return null;
        }

    }

    public List<String> getKeyNames() { return Arrays.asList("IndelType"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("IndelType", 1, VCFHeaderLineType.String, "Indel type description")); }

}
