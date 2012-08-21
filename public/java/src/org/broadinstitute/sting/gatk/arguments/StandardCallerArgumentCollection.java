package org.broadinstitute.sting.gatk.arguments;

import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 8/20/12
 * A collection of arguments that are common to the various callers.
 * This is pulled out so that every caller isn't exposed to the arguments from every other caller.
 */

public class StandardCallerArgumentCollection {
    /**
     * The expected heterozygosity value used to compute prior likelihoods for any locus. The default priors are:
     * het = 1e-3, P(hom-ref genotype) = 1 - 3 * het / 2, P(het genotype) = het, P(hom-var genotype) = het / 2
     */
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = UnifiedGenotyperEngine.HUMAN_SNP_HETEROZYGOSITY;

    @Argument(fullName = "genotyping_mode", shortName = "gt_mode", doc = "Specifies how to determine the alternate alleles to use for genotyping", required = false)
    public GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY;

    @Argument(fullName = "output_mode", shortName = "out_mode", doc = "Specifies which type of calls we should output", required = false)
    public UnifiedGenotyperEngine.OUTPUT_MODE OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY;

    /**
     * The minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls. Only genotypes with
     * confidence >= this threshold are emitted as called sites. A reasonable threshold is 30 for high-pass calling (this
     * is the default).
     */
    @Argument(fullName = "standard_min_confidence_threshold_for_calling", shortName = "stand_call_conf", doc = "The minimum phred-scaled confidence threshold at which variants should be called", required = false)
    public double STANDARD_CONFIDENCE_FOR_CALLING = 30.0;

    /**
     * This argument allows you to emit low quality calls as filtered records.
     */
    @Argument(fullName = "standard_min_confidence_threshold_for_emitting", shortName = "stand_emit_conf", doc = "The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)", required = false)
    public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;

    /**
     * When the UnifiedGenotyper is put into GENOTYPE_GIVEN_ALLELES mode it will genotype the samples using only the alleles provide in this rod binding
     */
    @Input(fullName="alleles", shortName = "alleles", doc="The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES", required=false)
    public RodBinding<VariantContext> alleles;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or GENOTYPE_GIVEN ALLELES),
     * then only this many alleles will be used.  Note that genotyping sites with many alternate alleles is both CPU and memory intensive and it
     * scales exponentially based on the number of alternate alleles.  Unless there is a good reason to change the default value, we highly recommend
     * that you not play around with this parameter.
     */
    @Advanced
    @Argument(fullName = "max_alternate_alleles", shortName = "maxAltAlleles", doc = "Maximum number of alternate alleles to genotype", required = false)
    public int MAX_ALTERNATE_ALLELES = 3;
}
