package org.broadinstitute.sting.gatk.arguments;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalcFactory;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.PrintStream;

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
     *
     * As of GATK 2.2 the genotyper can handle a very large number of events, so the default maximum has been increased to 6.
     */
    @Advanced
    @Argument(fullName = "max_alternate_alleles", shortName = "maxAltAlleles", doc = "Maximum number of alternate alleles to genotype", required = false)
    public int MAX_ALTERNATE_ALLELES = 6;

    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
    @Advanced
    @Argument(fullName = "p_nonref_model", shortName = "pnrm", doc = "Non-reference probability calculation model to employ", required = false)
    public AFCalcFactory.Calculation AFmodel = AFCalcFactory.Calculation.getDefaultModel();

    /**
     * If this fraction is greater is than zero, the caller will aggressively attempt to remove contamination through biased down-sampling of reads.
     * Basically, it will ignore the contamination fraction of reads for each alternate allele.  So if the pileup contains N total bases, then we
     * will try to remove (N * contamination fraction) bases for each alternate allele.
     */
    @Argument(fullName = "contamination_fraction_to_filter", shortName = "contamination", doc = "Fraction of contamination in sequencing data (for all samples) to aggressively remove", required = false)
    public double CONTAMINATION_FRACTION = DEFAULT_CONTAMINATION_FRACTION;
    public static final double DEFAULT_CONTAMINATION_FRACTION = 0.05;

    @Hidden
    @Argument(fullName = "logRemovedReadsFromContaminationFiltering", shortName="contaminationLog", required=false)
    public PrintStream contaminationLog = null;

    @Hidden
    @Argument(shortName = "logExactCalls", doc="x", required=false)
    public File exactCallsLog = null;

    public StandardCallerArgumentCollection() { }

    // Developers must remember to add any newly added arguments to the list here as well otherwise they won't get changed from their default value!
    public StandardCallerArgumentCollection(final StandardCallerArgumentCollection SCAC) {
        this.alleles = SCAC.alleles;
        this.GenotypingMode = SCAC.GenotypingMode;
        this.heterozygosity = SCAC.heterozygosity;
        this.MAX_ALTERNATE_ALLELES = SCAC.MAX_ALTERNATE_ALLELES;
        this.OutputMode = SCAC.OutputMode;
        this.STANDARD_CONFIDENCE_FOR_CALLING = SCAC.STANDARD_CONFIDENCE_FOR_CALLING;
        this.STANDARD_CONFIDENCE_FOR_EMITTING = SCAC.STANDARD_CONFIDENCE_FOR_EMITTING;
        this.CONTAMINATION_FRACTION = SCAC.CONTAMINATION_FRACTION;
        this.contaminationLog = SCAC.contaminationLog;
        this.exactCallsLog = SCAC.exactCallsLog;
        this.AFmodel = SCAC.AFmodel;
    }
}
