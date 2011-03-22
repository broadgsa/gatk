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

package org.broadinstitute.sting.playground.gatk.walkers;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
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


public class GenotypeAndValidateWalker extends RodWalker<GenotypeAndValidateWalker.CountedData, GenotypeAndValidateWalker.CountedData> implements TreeReducible<GenotypeAndValidateWalker.CountedData> {

    @Output(doc="File to which validated variants should be written", required=false)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName ="set_bam_truth", shortName ="bt", doc="Use the calls on the reads (bam file) as the truth dataset and validate the calls on the VCF", required=false)
    private boolean bamIsTruth = false;

    @Argument(fullName="minimum_base_quality_score", shortName="mbq", doc="Minimum base quality score for calling a genotype", required=false)
    private int mbq = -1;

    @Argument(fullName="maximum_deletion_fraction", shortName="deletions", doc="Maximum deletion fraction for calling a genotype", required=false)
    private double deletions = -1;

    @Argument(fullName="standard_min_confidence_threshold_for_calling", shortName="stand_call_conf", doc="he minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls", required=false)
    private double callConf = -1;

    @Argument(fullName="standard_min_confidence_threshold_for_emitting", shortName="stand_emit_conf", doc="the minimum phred-scaled Qscore threshold to emit low confidence calls", required=false)
    private double emitConf = -1;

    @Argument(fullName="condition_on_depth", shortName="depth", doc="Condition validation on a minimum depth of coverage by the reads", required=false)
    private int minDepth = -1;

    @Argument(fullName ="sample", shortName ="sn", doc="Name of the sample to validate (in case your VCF/BAM has more than one sample)", required=false)
    private String sample = "";




    private String compName = "alleles";
    private UnifiedGenotyperEngine snpEngine;
    private UnifiedGenotyperEngine indelEngine;

    public static class CountedData {
        private long nAltCalledAlt = 0L;
        private long nAltCalledRef = 0L;
        private long nRefCalledAlt = 0L;
        private long nRefCalledRef = 0L;
        private long nNotConfidentCalls = 0L;
        private long nUncovered = 0L;

        /**
         * Adds the values of other to this, returning this
         * @param other the other object
         */
        public void add(CountedData other) {
            nAltCalledAlt += other.nAltCalledAlt;
            nAltCalledRef += other.nAltCalledRef;
            nRefCalledAlt += other.nRefCalledAlt;
            nRefCalledRef += other.nRefCalledRef;
            nUncovered += other.nUncovered;
            nNotConfidentCalls += other.nNotConfidentCalls;
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
        if (vcfWriter != null) {
            Map<String, VCFHeader> header = VCFUtils.getVCFHeadersFromRodPrefix(getToolkit(), compName);
            Set<String> samples = SampleUtils.getSampleList(header, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
            Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(header.values(), logger);
            headerLines.add(new VCFHeaderLine("source", "GenotypeAndValidate"));
            vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
        }

        // Filling in SNP calling arguments for UG
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES;
        if (!bamIsTruth) uac.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        if (mbq >= 0) uac.MIN_BASE_QUALTY_SCORE = mbq;
        if (deletions >= 0) uac.MAX_DELETION_FRACTION = deletions;
        if (emitConf >= 0) uac.STANDARD_CONFIDENCE_FOR_EMITTING = emitConf;
        if (callConf >= 0) uac.STANDARD_CONFIDENCE_FOR_CALLING = callConf;

        snpEngine = new UnifiedGenotyperEngine(getToolkit(), uac);

        // Adding the INDEL calling arguments for UG
        uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.DINDEL;
        indelEngine = new UnifiedGenotyperEngine(getToolkit(), uac);

        // make sure we have callConf set to the threshold set by the UAC so we can use it later.
        callConf = uac.STANDARD_CONFIDENCE_FOR_CALLING;
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

        //todo - not sure I want this, may be misleading to filter extended indel events.
        if (isInsideExtendedIndel(vcComp,  ref))
            return counter;

        // Do not operate on variants that are not covered to the optional minimum depth
        if (!context.hasReads() || (minDepth > 0 && context.getBasePileup().getBases().length < minDepth)) {
            counter.nUncovered = 1L;
            return counter;
        }

        VariantCallContext call;
        if ( vcComp.isSNP() )
            call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context);
        else if ( vcComp.isIndel() ) {
            call = indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context);
        }
        else {
            logger.info("Not SNP or INDEL " + vcComp.getChr() + ":" + vcComp.getStart() + " " + vcComp.getAlleles());
            return counter;
        }

        if (bamIsTruth) {
            if (call.confidentlyCalled) {
                // If truth is a confident REF call
                if (call.isVariant()) {
                    if (vcComp.isVariant())
                        counter.nAltCalledAlt = 1L;  // todo -- may wanna check if the alts called are the same?
                    else
                        counter.nAltCalledRef = 1L;
                }
                // If truth is a confident ALT call
                else {
                    if (vcComp.isVariant())
                        counter.nRefCalledAlt = 1L;
                    else
                        counter.nRefCalledRef = 1L;
                }
            }
            else {
                counter.nNotConfidentCalls = 1L;
            }
        }
        else {
            if (!vcComp.hasAttribute("GV"))
                throw new UserException.BadInput("Variant has no GV annotation in the INFO field. " + vcComp.getChr() + ":" + vcComp.getStart());



            if (call.isCalledAlt(callConf)) {
                if (vcComp.getAttribute("GV").equals("T"))
                    counter.nAltCalledAlt = 1L;
                else
                    counter.nRefCalledAlt = 1L;
            }
            else if (call.isCalledRef(callConf)) {
                if (vcComp.getAttribute("GV").equals("T"))
                    counter.nAltCalledRef = 1L;
                else
                    counter.nRefCalledRef = 1L;
            }
            else {
                counter.nNotConfidentCalls = 1L;
            }
        }

        if (vcfWriter != null) {
            if (!vcComp.hasAttribute("callStatus")) {
                MutableVariantContext mvc = new MutableVariantContext(vcComp);
                mvc.putAttribute("callStatus", call.isCalledAlt(callConf) ? "ALT" : "REF" );
                vcfWriter.add(mvc, ref.getBase());
            }
            else
                vcfWriter.add(vcComp, ref.getBase());
        }
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

    public CountedData treeReduce( final CountedData sum1, final CountedData sum2) {
        sum2.add(sum1);
        return sum2;
    }

    public CountedData reduce( final CountedData mapValue, final CountedData reduceSum ) {
        reduceSum.add(mapValue);
        return reduceSum;
    }

    public void onTraversalDone( CountedData reduceSum ) {
        double ppv = 100 * ((double) reduceSum.nAltCalledAlt /( reduceSum.nAltCalledAlt + reduceSum.nRefCalledAlt));
        double npv = 100 * ((double) reduceSum.nRefCalledRef /( reduceSum.nRefCalledRef + reduceSum.nAltCalledRef));
        logger.info(String.format("Resulting Truth Table Output\n\n" +
                                  "---------------------------------------------------\n" +
                                  "\t\t|\tALT\t|\tREF\t\n"  +
                                  "---------------------------------------------------\n" +
                                  "called alt\t|\t%d\t|\t%d\n" +
                                  "called ref\t|\t%d\t|\t%d\n" +
                                  "---------------------------------------------------\n" +
                                  "positive predictive value: %f%%\n" +
                                  "negative predictive value: %f%%\n" +
                                  "---------------------------------------------------\n" +
                                  "not confident: %d\n" +
                                  "not covered: %d\n" +
                                  "---------------------------------------------------\n", reduceSum.nAltCalledAlt, reduceSum.nRefCalledAlt, reduceSum.nAltCalledRef, reduceSum.nRefCalledRef, ppv, npv, reduceSum.nNotConfidentCalls, reduceSum.nUncovered));
    }
}
