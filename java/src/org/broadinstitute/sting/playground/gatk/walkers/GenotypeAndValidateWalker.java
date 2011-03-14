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

import com.google.common.collect.ImmutableSet;
import org.broad.tribble.util.variantcontext.MutableVariantContext;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import sun.management.counter.Variability;

import java.util.*;

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


public class GenotypeAndValidateWalker extends RodWalker<GenotypeAndValidateWalker.CountedData, GenotypeAndValidateWalker.CountedData> {

    @Output(doc="File to which validated variants should be written", required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="minimum_base_quality_score", shortName="mbq", doc="Minimum base quality score for calling a genotype", required=false)
    private int mbq = -1;

    @Argument(fullName="maximum_deletion_fraction", shortName="deletions", doc="Maximum deletion fraction for calling a genotype", required=false)
    private double deletions = -1;

    @Argument(fullName="condition_on_depth", shortName="depth", doc="Condition validation on a minimum depth of coverage by the reads", required=false)
    private int minDepth = -1;


    private String compName = "alleles";
    private UnifiedGenotyperEngine snpEngine;
    private UnifiedGenotyperEngine indelEngine;

    public static class CountedData {
        private long numTP = 0L;
        private long numTN = 0L;
        private long numFP = 0L;
        private long numFN = 0L;
        private long numUncovered = 0L;
        private long numConfidentCalls = 0L;
        private long numNotConfidentCalls = 0L;

        /**
         * Adds the values of other to this, returning this
         * @param other the other object
         */
        public void add(CountedData other) {
            numTP += other.numTP;
            numTN += other.numTN;
            numFP += other.numFP;
            numFN += other.numFN;
            numUncovered += other.numUncovered;
            numNotConfidentCalls += other.numNotConfidentCalls;
            numConfidentCalls    += other.numConfidentCalls;
        }
    }



    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

        List<ReferenceOrderedDataSource> rodList = this.getToolkit().getRodDataSources();
        if ( rodList.size() != 1 )
            throw new UserException.BadInput("You should provide exactly one genotype VCF");
        if ( !rodList.get(0).getName().equals(compName))
            throw new UserException.BadInput("The ROD track has to be named \""+ compName +"\". Not " + rodList.get(0).getName());


        // Initialize VCF header
        Map<String, VCFHeader> header = VCFUtils.getVCFHeadersFromRodPrefix(getToolkit(), compName);
        Set<String> samples = SampleUtils.getSampleList(header, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(header.values(), logger);
        headerLines.add(new VCFHeaderLine("source", "GenotypeAndValidate"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));


        // Filling in SNP calling arguments for UG
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES;
        uac.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        if (mbq >= 0) uac.MIN_BASE_QUALTY_SCORE = mbq;
        if (deletions >= 0) uac.MAX_DELETION_FRACTION = deletions;
        snpEngine = new UnifiedGenotyperEngine(getToolkit(), uac);

        // Adding the INDEL calling arguments for UG
        uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.DINDEL;
        indelEngine = new UnifiedGenotyperEngine(getToolkit(), uac);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public CountedData map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        final CountedData counter = new CountedData();

        // For some reason RodWalkers get map calls with null trackers
        if( tracker == null )
            return counter;

        VariantContext vcComp = tracker.getVariantContext(ref, compName, null, context.getLocation(), false);
        if( vcComp == null )
            return counter;

        // Do not operate on variants that are not covered to the optional minimum depth
        if (!context.hasReads() || minDepth > 0 && context.getBasePileup().getBases().length < minDepth) {
            counter.numUncovered = 1L;
            return counter;
        }

        if ((vcComp.getFilters().contains("TP") && vcComp.getFilters().contains("FP")) || (!vcComp.getFilters().contains("TP") && !vcComp.getFilters().contains("FP")))
            throw new UserException.BadInput("Variant has the wrong filter annotation -- either missing FP/TP or has both. " + vcComp.getChr() + ":" + vcComp.getStart());

        VariantCallContext call;
        if ( vcComp.isSNP() )
            call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context);
        else if ( vcComp.isIndel() ) {
            call = indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context);
            if (call.vc == null) // variant context will be null on an extended indel event and I just want to call it one event.
                return counter;
        }
        else {
            logger.info("Not SNP or INDEL " + vcComp.getChr() + ":" + vcComp.getStart() + " " + vcComp.getAlleles());
            return counter;
        }

        if (!call.confidentlyCalled) {
            counter.numNotConfidentCalls = 1L;
            if (vcComp.getFilters().contains("TP"))
                counter.numFN = 1L;
            else
                counter.numTN = 1L;
        }
        else {
            counter.numConfidentCalls = 1L;
            if (vcComp.getFilters().contains("TP"))
                counter.numTP = 1L;
            else
                counter.numFP = 1L;
        }
        if (!vcComp.hasAttribute("callStatus")) {
            MutableVariantContext mvc = new MutableVariantContext(vcComp);
            mvc.putAttribute("callStatus", call.confidentlyCalled ? "confident" : "notConfident" );
            vcfWriter.add(mvc, ref.getBase());
        }
        else
            vcfWriter.add(vcComp, ref.getBase());
        return counter;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public CountedData reduceInit() {
        return new CountedData();
    }

    public CountedData reduce( final CountedData mapValue, final CountedData reduceSum ) {
        reduceSum.add(mapValue);
        return reduceSum;
    }

    public void onTraversalDone( CountedData reduceSum ) {
        logger.info("TP = " + reduceSum.numTP);
        logger.info("TN = " + reduceSum.numTN);
        logger.info("FP = " + reduceSum.numFP);
        logger.info("FN = " + reduceSum.numFN);
        logger.info("PPV = " + ((double) reduceSum.numTP /( reduceSum.numTP + reduceSum.numFP)));
        logger.info("NPV = " + ((double) reduceSum.numTN /( reduceSum.numTN + reduceSum.numFN)));
        logger.info("Uncovered = " + reduceSum.numUncovered);
        logger.info("confidently called = " + reduceSum.numConfidentCalls);
        logger.info("not confidently called = " + reduceSum.numNotConfidentCalls );
    }
}
