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

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.*;

/**
 * Validates the calls on a ROD track using a BAM dataset.
 *
 * @author carneiro
 * @since Mar 3, 2011
 * @help.summary Validates the calls on a ROD track using a BAM dataset.
 */

//@Requires(value={DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_ORDERED_DATA})
@Allows(value={DataSource.READS, DataSource.REFERENCE})

public class GenotypeAndValidateWalker extends RodWalker<GenotypeAndValidateWalker.CountedData, GenotypeAndValidateWalker.CountedData> {

    @Argument(fullName="minimum_base_quality_score", shortName="mbq", doc="Minimum base quality score for calling a genotype", required=false)
    private int mbq = -1;

    @Argument(fullName="maximum_deletion_fraction", shortName="deletions", doc="Maximum deletion fraction for calling a genotype", required=false)
    private double deletions = -1;


    @Output( doc="File to write the two-way truth table", required=false)
    private PrintStream printStream = null;

    private String compName = "alleles";
    private UnifiedGenotyperEngine snpEngine;
    private UnifiedGenotyperEngine indelEngine;


    public static enum VARIANT_TYPE {
        TRUE_POSITIVE,
        FALSE_POSITIVE,
    }

    public static class CountedData {
        private long numTP = 0L;
        private long numTN = 0L;
        private long numFP = 0L;
        private long numFN = 0L;
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

        // Filling in SNP calling arguments for UG
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES;
        uac.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        if (mbq > 0) uac.MIN_BASE_QUALTY_SCORE = mbq;
        if (deletions > 0) uac.MAX_DELETION_FRACTION = deletions;
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
        vcComp = vcComp.subContextFromGenotypes( vcComp.getGenotypes().values() );
        if ((vcComp.getFilters().contains("TP") && vcComp.getFilters().contains("FP")) || (!vcComp.getFilters().contains("TP") && !vcComp.getFilters().contains("FP")))
            throw new UserException.BadInput("Variant has the wrong filter annotation -- either missing FP/TP or has both. " + vcComp.getChr() + ":" + vcComp.getStart());

        VariantCallContext call;
        if ( vcComp.isSNP() )
            call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context);

        else if ( vcComp.isIndel() )
            call = indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context);

        else
            return counter;

        VARIANT_TYPE variantType;
        if (vcComp.getFilters().contains("TP"))
            variantType = VARIANT_TYPE.TRUE_POSITIVE;
        else
            variantType = VARIANT_TYPE.FALSE_POSITIVE;

        if (!call.confidentlyCalled) {
            counter.numNotConfidentCalls = 1L;
            if (variantType == VARIANT_TYPE.TRUE_POSITIVE)
                counter.numFN = 1L;
            else
                counter.numTN = 1L;
        }
        else {
            counter.numConfidentCalls = 1L;
            if (variantType == VARIANT_TYPE.TRUE_POSITIVE)
                counter.numTP = 1L;
            else
                counter.numFP = 1L;
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

    public CountedData reduce( final CountedData mapValue, final CountedData reduceSum ) {
        reduceSum.add(mapValue);
        return reduceSum;
    }

    public void onTraversalDone( CountedData reduceSum ) {
        logger.info("TP = " + reduceSum.numTP);
        logger.info("TN = " + reduceSum.numTN);
        logger.info("FP = " + reduceSum.numFP);
        logger.info("FN = " + (reduceSum.numFN) );
        logger.info("confidently called = " + reduceSum.numConfidentCalls);
        logger.info("not confidently called = " + reduceSum.numNotConfidentCalls );
    }
}
