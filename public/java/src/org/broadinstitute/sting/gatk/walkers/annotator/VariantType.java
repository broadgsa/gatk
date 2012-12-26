package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.IndelUtils;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Assigns a roughly correct category of the variant type (SNP, MNP, insertion, deletion, etc.)
 */
public class VariantType extends InfoFieldAnnotation implements ExperimentalAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {

        StringBuffer type = new StringBuffer("");
        if ( vc.isVariant() && !vc.isBiallelic() )
            type.append("MULTIALLELIC_");

        if ( !vc.isIndel() ) {
            type.append(vc.getType().toString());
        } else {
            if (vc.isSimpleInsertion())
                type.append("INSERTION.");
            else if (vc.isSimpleDeletion())
                type.append("DELETION.");
            else
                type.append("COMPLEX.");
            ArrayList<Integer> inds = IndelUtils.findEventClassificationIndex(vc, ref);
            type.append(IndelUtils.getIndelClassificationName(inds.get(0)));

            for (int i = 1; i < inds.size(); i++ ) {
                type.append(".");
                type.append(IndelUtils.getIndelClassificationName(inds.get(i)));
            }
        }

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%s", type));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("VariantType"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("VariantType", 1, VCFHeaderLineType.String, "Variant type description")); }

}
