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
import org.broad.tribble.util.variantcontext.MutableVariantContext;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

import static org.broadinstitute.sting.utils.IndelUtils.isInsideExtendedIndel;

/**
 * Validates the calls on a ROD track using a BAM dataset.
 *
 * @author carneiro
 * @since Mar 3, 2011
 * @help.summary Validates the calls on a ROD track using a BAM dataset.
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

    public static class Data {
        Collection<Datum> values;

        public Data() { this(new LinkedList<Datum>()); }
        public Data(Collection<Datum> data) { this.values = data; }

        final public static Data EMPTY_DATA = new Data(Collections.<Datum>emptyList());
    }

    public static class Datum implements Comparable<Datum> {
        String rgID, sample;
        GenotypeLikelihoods pl;
        VariantContext.Type siteType;
        Genotype.Type genotypeType;

        @Override
        public int compareTo(Datum o) {
            int bySample = sample.compareTo(o.sample);
            int byRG = rgID.compareTo(o.rgID);
            return bySample != 0 ? bySample : byRG;
        }

        public Datum(String sample, String rgID, GenotypeLikelihoods pl, VariantContext.Type siteType, Genotype.Type genotypeType) {
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
        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        logger.info("Samples: " + samples);

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
    public boolean generateExtendedEvents() { return true; }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------
    public Data map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        if ( tracker == null || tracker.getNBoundRodTracks() == 0 )
            return Data.EMPTY_DATA;

        VariantContext vcComp = SNPGenotypeLikelihoodsCalculationModel.getSNPVCFromAllelesRod(tracker, ref, false, logger);
        if( vcComp == null )
            return Data.EMPTY_DATA;

        //todo - not sure I want this, may be misleading to filter extended indel events.
        if (isInsideExtendedIndel(vcComp,  ref))
            return Data.EMPTY_DATA;

        Data data = new Data();
        for ( String sample : samples ) {
            Genotype compGT = getGenotype(tracker, ref, sample, COMP_NAME);
            if ( compGT == null || compGT.isNoCall() )
                continue;

            Map<SAMReadGroupRecord,AlignmentContext> byRG = AlignmentContextUtils.splitContextByReadGroup(context, getToolkit().getSAMFileHeader().getReadGroups());
            //byRG.put(new SAMReadGroupRecord("ALL"), context);
            for ( Map.Entry<SAMReadGroupRecord, AlignmentContext> rgAC : byRG.entrySet() ) {
                VariantCallContext call;
                if ( vcComp.isIndel() ) {
                    call = indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, rgAC.getValue());
                } else {
                    call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, rgAC.getValue());
                }

                if ( call == null )
                    throw new ReviewedStingException("Unexpected genotyping failure " + sample + " at " + ref.getLocus() + " call " + call);

                Genotype rgGT = call.getGenotype(sample);

                if ( rgGT != null && ! rgGT.isNoCall() && rgGT.getLikelihoods().getAsVector() != null ) {
                    Datum d = new Datum(sample, rgAC.getKey().getReadGroupId(), rgGT.getLikelihoods(), vcComp.getType(), compGT.getType());
                    data.values.add(d);
                }
            }
        }

        return data;
    }

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
        List<String> fields = Arrays.asList("sample", "rg", "siteType", "pls", "comp", "pGGivenDType", "pGGivenD");
        out.println(Utils.join("\t", fields));

        double[] counts = new double[]{1, 1, 1};
        for ( Datum d : data.values ) { counts[d.genotypeType.ordinal()-1]++; }
        double sum = MathUtils.sum(counts);
        logger.info(String.format("Types %s %s %s", Genotype.Type.values()[1], Genotype.Type.values()[2], Genotype.Type.values()[3]));
        logger.info(String.format("Counts %.0f %.0f %.0f %.0f", counts[0], counts[1], counts[2], sum));
        double[] log10priors = new double[]{Math.log10(counts[0] / sum), Math.log10(counts[1] / sum), Math.log10(counts[2] / sum)};
        logger.info(String.format("Priors %.2f %.2f %.2f", log10priors[0], log10priors[1], log10priors[2]));

        for ( Datum d : data.values ) {
            double[] log10pGGivenD = d.pl.getAsVector().clone();
            for ( int i = 0; i < log10priors.length; i++ ) log10pGGivenD[i] += log10priors[i];
            double[] pOfGGivenD = MathUtils.normalizeFromLog10(log10pGGivenD, false);
            for ( int i = 0; i < pGNames.size(); i++ ) {
                int q = QualityUtils.probToQual(pOfGGivenD[i], Math.pow(10.0, -9.9));
                if ( q > 1 ) { // tons of 1s, and not interesting
                    out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%d%n",
                            d.sample, d.rgID, d.siteType, d.pl.getAsString(), d.genotypeType.toString(),
                            pGNames.get(i), q);
                }
            }
        }
    }
}
