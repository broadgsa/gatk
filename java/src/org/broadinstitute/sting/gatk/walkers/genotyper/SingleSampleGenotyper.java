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

    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File to which metrics should be written", required = false)
    public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold on which variants should be filtered", required = false)
    public Double LOD_THRESHOLD = Double.MIN_VALUE;

    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confidient genotypes or just the variants?", required = false)
    public boolean GENOTYPE = false;

    // Control how the genotype hypotheses are weighed
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = GenotypeLikelihoods.HUMAN_HETEROZYGOSITY;

    // todo -- remove dbSNP awareness until we understand how it should be used
    //@Argument(fullName = "priors_dbsnp", shortName = "pdbsnp", doc = "Comma-separated prior likelihoods for dbSNP loci (homref,het,homvar)", required = false)
    //public String PRIORS_DBSNP = "0.999,1e-3,1e-5";

    @Argument(fullName = "keepQ0Bases", shortName = "keepQ0Bases", doc = "If true, then Q0 bases will be included in the genotyping calculation, and treated as Q1 -- this is really not a good idea", required = false)
    public boolean keepQ0Bases = false;


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
        NewHotnessGenotypeLikelihoods G = new NewHotnessGenotypeLikelihoods(ref, GenotypeLikelihoods.HUMAN_HETEROZYGOSITY);
        G.filterQ0Bases(! keepQ0Bases);
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        G.add(pileup, true);
        G.validate();

        // lets setup the locus
        List<Genotype> lst = new ArrayList<Genotype>();
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            lst.add(new BasicGenotype(pileup.getLocation(), g.toString(),new BayesianConfidenceScore(G.getLikelihood(g))));
        }

        return new SSGGenotypeCall(! GENOTYPE, ref, 2, lst, G.getPosteriors(), pileup);
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


