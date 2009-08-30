package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.BasicGenotype;
import org.broadinstitute.sting.utils.genotype.confidence.BayesianConfidenceScore;

import java.io.File;
import java.util.List;
import java.util.ArrayList;

@ReadFilters(ZeroMappingQualityReadFilter.class)
public class SingleSampleGenotyper extends LocusWalker<SSGGenotypeCall, GenotypeWriter> {
    // Control output settings
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = false)
    public File VARIANTS_FILE = null;

    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File format to be used", required = false)
    public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold on which variants should be filtered", required = false)
    public Double LOD_THRESHOLD = Double.MIN_VALUE;

    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confident genotypes or just the variants?", required = false)
    public boolean GENOTYPE = false;

    // Control how the genotype hypotheses are weighed
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = DiploidGenotypePriors.HUMAN_HETEROZYGOSITY;

    @Argument(fullName = "baseModel", shortName = "m", doc = "Base substitution model to employ -- EMPIRICAL is the recommended default, but it's possible to select the ONE_STATE and THREE_STATE models for comparison purposes", required = false)
    public BaseMismatchModel baseModel = BaseMismatchModel.EMPIRICAL;

    @Argument(fullName = "verbose", shortName = "v", doc = "EXPERIMENTAL", required = false)
    public boolean VERBOSE = false;

    @Argument(fullName = "platform", shortName = "pl", doc = "Causes the genotyper to assume that reads without PL header TAG are this platform.  Defaults to null, indicating that the system will throw a runtime exception when such reads are detected", required = false)
    public EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform defaultPlatform = null;

     /**
     * Filter out loci to ignore (at an ambiguous base in the reference or a locus with zero coverage).
     *
     * @param tracker the meta data tracker
     * @param ref     the reference base
     * @param context contextual information around the locus
     *
     * @return true if we should look at this locus, false otherwise
     */
    public boolean filter(RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref) != -1 && context.getReads().size() != 0);
    }

    /** Initialize the walker with some sensible defaults */
    public void initialize() {
        // nothing to do
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param context contextual information around the locus
     */
    public SSGGenotypeCall map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext context) {
        char ref = refContext.getBase();
        DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, heterozygosity, DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);

        // setup GenotypeLike object
        GenotypeLikelihoods gl = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);

        gl.setVerbose(VERBOSE);
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        gl.add(pileup, true);
        gl.validate();

        // lets setup the locus
        List<Genotype> lst = new ArrayList<Genotype>();
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            lst.add(new BasicGenotype(pileup.getLocation(), g.toString(),new BayesianConfidenceScore(gl.getLikelihood(g))));
        }

        //System.out.printf("At %s%n", new SSGGenotypeCall(! GENOTYPE, ref, 2, lst, g.getPosteriors(), pileup));
        return new SSGGenotypeCall(! GENOTYPE, ref, 2, lst, gl.getPosteriors(), pileup);
    }

    /**
     * Determine whether we're at a Hapmap site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a Hapmap site, false if otherwise
     */
    private boolean isHapmapSite(RefMetaDataTracker tracker) {
        return tracker.lookup("hapmap", null) != null;
    }

    /**
     * Determine whether we're at a dbSNP site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a dbSNP site, false if otherwise
     */
    private boolean isDbSNPSite(RefMetaDataTracker tracker) {
        return tracker.lookup("dbsnp", null) != null;
    }

    /**
     * Initialize values appropriately for the reduce step.
     *
     * @return an empty string
     */
    public GenotypeWriter reduceInit() {
        if ( VARIANTS_FILE != null )
            return GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), VARIANTS_FILE);
        else
            return GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), out);
    }

    /**
     * If we've found a LOD >= 5 variant, output it to disk.
     *
     * @param call an GenotypeCall object for the variant.
     * @param sum  accumulator for the reduce.
     *
     * @return an empty string
     */
    public GenotypeWriter reduce(SSGGenotypeCall call, GenotypeWriter sum) {
        if (call != null && (GENOTYPE || call.isVariation())) {
            //System.out.printf("Score %s%n", call.getConfidenceScore());
            if (call.getConfidenceScore().getScore() > LOD_THRESHOLD) {
                //System.out.printf("Call %s%n", call);                
                sum.addGenotypeCall(call);
            }
        }
        return sum;
    }

    /** Close the variant file. */
    public void onTraversalDone(GenotypeWriter sum) {
        sum.close();
    }
}


