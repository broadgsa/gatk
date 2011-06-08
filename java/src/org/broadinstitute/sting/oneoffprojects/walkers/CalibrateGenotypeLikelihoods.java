/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers;

import net.sf.samtools.SAMReadGroupRecord;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.*;

import static org.broadinstitute.sting.utils.IndelUtils.isInsideExtendedIndel;

/**
 * Computes raw GL calibration data for read groups in BAMs against a comp VCF track of genotypes
 *
 * @author depristo
 * @since May, 2011
 * @help.summary Computes raw GL calibration data for read groups in BAMs against a comp VCF track of genotypes
 */

@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="alleles",type=VariantContext.class))
@Allows(value={DataSource.READS, DataSource.REFERENCE})

// Ugly fix because RodWalkers don't have access to reads
@By(DataSource.REFERENCE)
@Reference(window=@Window(start=-200,stop=200))
public class CalibrateGenotypeLikelihoods extends RodWalker<CalibrateGenotypeLikelihoods.Data, CalibrateGenotypeLikelihoods.Data> implements TreeReducible<CalibrateGenotypeLikelihoods.Data> {
    public static final String COMP_NAME = "alleles";

    @Argument(fullName="minimum_base_quality_score", shortName="mbq", doc="Minimum base quality score for calling a genotype", required=false)
    private int mbq = -1;

    @Argument(fullName="maximum_deletion_fraction", shortName="deletions", doc="Maximum deletion fraction for calling a genotype", required=false)
    private double deletions = -1;

    //@Argument(fullName="standard_min_confidence_threshold_for_calling", shortName="stand_call_conf", doc="the minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls", required=false)
    private double callConf = 0;

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    Set<String> samples;

    /**
     * Trivial wrapper class.  Data is a collection of Datum.
     */
    public static class Data {
        Collection<Datum> values;

        public Data() { this(new LinkedList<Datum>()); }
        public Data(Collection<Datum> data) { this.values = data; }

        final public static Data EMPTY_DATA = new Data(Collections.<Datum>emptyList());
    }

    /**
     * The raw datapoints we are tracking for a specific site for a specific sample.
     * read group id and sample name.  The PL object.
     * the ref and alt alleles. The type of the variant context.  And the genotype of the
     * comp. track at this site.
     */
    public static class Datum implements Comparable<Datum> {
        final String rgID, sample;
        final GenotypeLikelihoods pl;
        final String ref, alt;
        final VariantContext.Type siteType;
        final Genotype.Type genotypeType;

        @Override
        public int compareTo(Datum o) {
            int bySample = sample.compareTo(o.sample);
            int byRG = rgID.compareTo(o.rgID);
            return bySample != 0 ? bySample : byRG;
        }

        public Datum(String ref, String alt, String sample, String rgID, GenotypeLikelihoods pl, VariantContext.Type siteType, Genotype.Type genotypeType) {
            this.ref = ref;
            this.alt = alt;
            this.sample = sample;
            this.rgID = rgID;
            this.pl = pl;
            this.siteType = siteType;
            this.genotypeType = genotypeType;
        }
    }

    private UnifiedGenotyperEngine snpEngine;
    private UnifiedGenotyperEngine indelEngine;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        // We only operate over the samples in the BAM file
        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        logger.info("Samples: " + samples);
        if ( samples.size() > 1 ) // todo -- remove me when we support multiple samples
            throw new UserException.BadInput("CalibrateGenotypeLikelihoods does not currently support comparison of multiple samples simulatenously.  To enable, see TODO in code");

        List<ReferenceOrderedDataSource> rodList = this.getToolkit().getRodDataSources();
        if ( rodList.size() != 1 )
            throw new UserException.BadInput("You should provide exactly one genotype VCF");
        if ( !rodList.get(0).getName().equals(COMP_NAME))
            throw new UserException.BadInput("The ROD track has to be named \""+ COMP_NAME +"\". Not " + rodList.get(0).getName());

        // Filling in SNP calling arguments for UG
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES;
        uac.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        uac.NO_SLOD = true;
        if (mbq >= 0) uac.MIN_BASE_QUALTY_SCORE = mbq;
        if (deletions >= 0) uac.MAX_DELETION_FRACTION = deletions;
        uac.STANDARD_CONFIDENCE_FOR_CALLING = callConf;
        uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;
        snpEngine = new UnifiedGenotyperEngine(getToolkit(), uac);

        // Adding the INDEL calling arguments for UG
        uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.INDEL;
        indelEngine = new UnifiedGenotyperEngine(getToolkit(), uac);
    }

    @Override
    // todo -- remove me when the new indel genotyping is done
    public boolean generateExtendedEvents() { return true; }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------
    public Data map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        if ( tracker == null || tracker.getNBoundRodTracks() == 0 )
            return Data.EMPTY_DATA;

        // Grabs a usable VariantContext from the Alleles ROD
        VariantContext vcComp = SNPGenotypeLikelihoodsCalculationModel.getSNPVCFromAllelesRod(tracker, ref, false, logger);
        if( vcComp == null )
            return Data.EMPTY_DATA;

        Data data = new Data();
        for ( String sample : samples ) {
            // What's the genotype of our sample at this record?
            Genotype compGT = getGenotype(tracker, ref, sample, COMP_NAME);
            if ( compGT == null || compGT.isNoCall() )
                continue;

            // For each read group
            // todo -- this only works with a single sample right now.  For multi-sample BAMs
            // todo -- this loop needs to be refactored so that the spliting by read group only happens once
            // todo -- and the read groups appropriate to each sample is used.
            Map<SAMReadGroupRecord,AlignmentContext> byRG = AlignmentContextUtils.splitContextByReadGroup(context, getToolkit().getSAMFileHeader().getReadGroups());
            //byRG.put(new SAMReadGroupRecord("ALL"), context);     // uncomment to include a synthetic RG for all RG for the sample
            for ( Map.Entry<SAMReadGroupRecord, AlignmentContext> rgAC : byRG.entrySet() ) {
                VariantCallContext call;
                if ( vcComp.isIndel() ) {
                    throw new UserException.BadInput("CalibrateGenotypeLikelihoods does not currently support indel GL calibration.  This capability needs to be tested and verified to be working with the new genotyping code for indels in UG");
                    //call = indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, rgAC.getValue());
                } else {
                    call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, rgAC.getValue());
                }

                if ( call == null )
                    throw new ReviewedStingException("Unexpected genotyping failure " + sample + " at " + ref.getLocus() + " call " + call);

                Genotype rgGT = call.getGenotype(sample);

                if ( rgGT != null && ! rgGT.isNoCall() && rgGT.getLikelihoods().getAsVector() != null ) {
                    Datum d = new Datum(vcComp.getReference().getBaseString(), vcComp.getAlternateAllele(0).getBaseString(),
                            sample, rgAC.getKey().getReadGroupId(), rgGT.getLikelihoods(), vcComp.getType(), compGT.getType());
                    data.values.add(d);
                }
            }
        }

        return data;
    }

    /**
     * Convenience function that determines the genotype in the comp VC for sample
     *
     * @param tracker
     * @param ref
     * @param sample
     * @param rod
     * @return
     */
    private Genotype getGenotype(RefMetaDataTracker tracker, ReferenceContext ref, String sample, String rod) {
        for ( VariantContext vc : tracker.getVariantContexts(ref, rod, null, ref.getLocus(), true, false) ) {
            if ( vc.isNotFiltered() && vc.hasGenotype(sample) )
                return vc.getGenotype(sample);
            else
                return null;
        }

        return null;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Data reduceInit() {
        return new Data();
    }

    @Override
    public Data treeReduce( final Data sum1, final Data sum2) {
        sum2.values.addAll(sum1.values);
        return sum2;
    }

    @Override
    public Data reduce( final Data mapValue, final Data reduceSum ) {
        return treeReduce(mapValue, reduceSum);
    }

    @Override
    public void onTraversalDone(Data data) {
        // print the header
        List<String> pGNames = Arrays.asList("QofAAGivenD", "QofABGivenD", "QofBBGivenD");
        List<String> fields = Arrays.asList("sample", "rg", "ref", "alt", "siteType", "pls", "comp", "pGGivenDType", "pGGivenD");
        out.println(Utils.join("\t", fields));

        // determine the priors by counting all of the events we've seen in comp
        double[] counts = new double[]{1, 1, 1};
        for ( Datum d : data.values ) { counts[d.genotypeType.ordinal()-1]++; }
        double sum = MathUtils.sum(counts);
        logger.info(String.format("Types %s %s %s", Genotype.Type.values()[1], Genotype.Type.values()[2], Genotype.Type.values()[3]));
        logger.info(String.format("Counts %.0f %.0f %.0f %.0f", counts[0], counts[1], counts[2], sum));
        double[] log10priors = new double[]{Math.log10(counts[0] / sum), Math.log10(counts[1] / sum), Math.log10(counts[2] / sum)};
        logger.info(String.format("Priors %.2f %.2f %.2f", log10priors[0], log10priors[1], log10priors[2]));

        // emit the molten data set
        for ( Datum d : data.values ) {
            double[] log10pGGivenD = d.pl.getAsVector().clone();
            for ( int i = 0; i < log10priors.length; i++ ) log10pGGivenD[i] += log10priors[i];
            double[] pOfGGivenD = MathUtils.normalizeFromLog10(log10pGGivenD, false);
            for ( int i = 0; i < pGNames.size(); i++ ) {
                int q = QualityUtils.probToQual(pOfGGivenD[i], Math.pow(10.0, -9.9));
                if ( q > 1 ) { // tons of 1s, and not interesting
                    out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d%n",
                            d.sample, d.rgID, d.ref, d.alt, d.siteType, d.pl.getAsString(), d.genotypeType.toString(),
                            pGNames.get(i), q);
                }
            }
        }
    }
}
