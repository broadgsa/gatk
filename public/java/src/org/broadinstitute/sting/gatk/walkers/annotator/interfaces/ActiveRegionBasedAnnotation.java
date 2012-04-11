package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.List;
import java.util.Map;

// TODO -- make this an abstract class when we move away from InfoFieldAnnotation
public interface ActiveRegionBasedAnnotation extends AnnotationType {
    // return annotations for the given contexts split by sample and then allele
    public abstract Map<String, Object> annotate(final Map<String, Map<Allele, List<GATKSAMRecord>>> stratifiedContexts, final VariantContext vc);

    // return the descriptions used for the VCF INFO meta field
    public abstract List<VCFInfoHeaderLine> getDescriptions();
}