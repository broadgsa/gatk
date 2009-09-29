package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.io.File;
import java.util.*;

import net.sf.samtools.SAMReadGroupRecord;


@ReadFilters({ZeroMappingQualityReadFilter.class, MissingReadGroupFilter.class})
public class UnifiedGenotyper extends LocusWalker<Integer, Integer> {

    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = false)
    public File VARIANTS_FILE = null;

    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File format to be used", required = false)
    public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;


    // the model used for calculating genotypes
    private GenotypeCalculationModel gcm;

    // output writer
    private GenotypeWriter writer;


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

    /** Enable deletions in the pileup */
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /** Initialize the walker with some sensible defaults */
    public void initialize() {

        // get all of the unique sample names
        HashSet<String> samples = new HashSet<String>();
        List<SAMReadGroupRecord> readGroups = getToolkit().getSAMFileHeader().getReadGroups();
        for ( SAMReadGroupRecord readGroup : readGroups )
            samples.add(readGroup.getSample());

        // print them out for debugging (need separate loop to ensure uniqueness)
        for ( String sample : samples )
            logger.debug("SAMPLE: " + sample);

        // create the output writer stream
        if ( VARIANTS_FILE != null )
            writer = GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), VARIANTS_FILE);
        else
            writer = GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), out);

        gcm = GenotypeCalculationModelFactory.makeGenotypeCalculation(UAC, samples, writer);
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param context contextual information around the locus
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext context) {
        char ref = Character.toUpperCase(refContext.getBase());
        if ( !BaseUtils.isRegularBase(ref) )
            return 0;

        // an optimization to speed things up when overly covered
        if ( UAC.MAX_READS_IN_PILEUP > 0 && context.getReads().size() > UAC.MAX_READS_IN_PILEUP )
            return 0;

        DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, UAC.heterozygosity, DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);
        return ( gcm.calculateGenotype(tracker, ref, context, priors) ?  1 : 0 );
    }

    /**
     * Determine whether we're at a Hapmap site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a Hapmap site, false if otherwise
     */
    private static boolean isHapmapSite(RefMetaDataTracker tracker) {
        return tracker.getTrackData("hapmap", null) != null;
    }

    /**
     * Determine whether we're at a dbSNP site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a dbSNP site, false if otherwise
     */
    private static boolean isDbSNPSite(RefMetaDataTracker tracker) {
        return rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null)) != null;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /** Close the variant file. */
    public void onTraversalDone(Integer sum) {
        logger.info("Processed " + sum + " loci that are callable for SNPs");
        writer.close();
    }
}
