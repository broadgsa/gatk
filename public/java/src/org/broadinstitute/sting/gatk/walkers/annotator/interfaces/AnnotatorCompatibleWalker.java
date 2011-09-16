package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.List;

public interface AnnotatorCompatibleWalker {

    // getter methods for various used bindings
    public abstract RodBinding<VariantContext> getVariantRodBinding();
    public abstract RodBinding<VariantContext> getSnpEffRodBinding();
    public abstract RodBinding<VariantContext> getDbsnpRodBinding();
    public abstract List<RodBinding<VariantContext>> getCompRodBindings();
    public abstract List<RodBinding<VariantContext>> getResourceRodBindings();
}
