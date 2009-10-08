/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MissingReadGroupFilter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.ArgumentCollection;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;

import java.io.File;
import java.util.HashSet;
import java.util.List;


@ReadFilters({ZeroMappingQualityReadFilter.class, MissingReadGroupFilter.class})
public class UnifiedGenotyper extends LocusWalker<List<GenotypeCall>, Integer> {

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
            writer = GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), VARIANTS_FILE,
                                                  "UnifiedGenotyper",
                                                  this.getToolkit().getArguments().referenceFile.getName());
        else
            writer = GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), out);

        gcm = GenotypeCalculationModelFactory.makeGenotypeCalculation(UAC, samples, writer, logger);
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param context contextual information around the locus
     */
    public List<GenotypeCall> map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext context) {
        char ref = Character.toUpperCase(refContext.getBase());
        if ( !BaseUtils.isRegularBase(ref) )
            return null;

        // an optimization to speed things up when overly covered
        if ( UAC.MAX_READS_IN_PILEUP > 0 && context.getReads().size() > UAC.MAX_READS_IN_PILEUP )
            return null;

        DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, UAC.heterozygosity, DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);
        return gcm.calculateGenotype(tracker, ref, context, priors);
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

    public Integer reduce(List<GenotypeCall> value, Integer sum) {
        return sum + 1;
    }

    /** Close the variant file. */
    public void onTraversalDone(Integer sum) {
        logger.info("Processed " + sum + " loci that are callable for SNPs");
        writer.close();
    }
}
