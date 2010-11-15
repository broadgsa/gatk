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

package org.broadinstitute.sting.playground.gatk.walkers.validationgenotyper;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.TreeSet;

/**
 * This is the validation genotyper. It is not ready to be used by anybody.
 *
 * @author rpoplin
 * @since Oct 20, 2010
 * @help.summary This is the validation genotyper. It is not ready to be used by anybody.
 */

@Allows(value={DataSource.READS, DataSource.REFERENCE})
public class ValidationGenotyper extends LocusWalker<ValidationGenotyper.CountedData, ValidationGenotyper.CountedData> implements TreeReducible<ValidationGenotyper.CountedData> {

        @Output( doc="The output filtered VCF file", required=true)
    private PrintStream printStream = null;

    public static class CountedData {
        private long numTP = 0L;
        private long numTN = 0L;
        private long numFP = 0L;
        private long numFilteringFN = 0L;
        private long numCallingFN = 0L;

        /**
         * Adds the values of other to this, returning this
         * @param other the other object
         */
        public void add(CountedData other) {
            numTP += other.numTP;
            numTN += other.numTN;
            numFP += other.numFP;
            numFilteringFN += other.numFilteringFN;
            numCallingFN += other.numCallingFN;
        }
    }

    public enum VARIANT_STATUS {
        CALLED,
        FILTERED,
        MISSING
    }

    final private ArrayList<String> evalNames = new ArrayList<String>();
    final private ArrayList<String> compNames = new ArrayList<String>();
    final private TreeSet<String> overlappingSamples = new TreeSet<String>();
    private UnifiedGenotyperEngine engine;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

        for( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if( d.getName().toLowerCase().startsWith("eval") ) {
                evalNames.add( d.getName() );
            } else if( d.getName().toLowerCase().startsWith("comp") ) {
                compNames.add( d.getName() );
            } else {
                throw new UserException.BadInput("Don't know what to do with input ROD track named: " + d.getName());
            }
        }

        if( evalNames.size() != 1 || compNames.size() != 1 ) {
            throw new UserException.BadInput("Expecting to see exactly one eval track and exactly one comp track");
        }

        final TreeSet<String> evalSamples = new TreeSet<String>();
        evalSamples.addAll(SampleUtils.getUniqueSamplesFromRods(getToolkit(), evalNames));
        final TreeSet<String> compSamples = new TreeSet<String>();
        compSamples.addAll(SampleUtils.getUniqueSamplesFromRods(getToolkit(), compNames));
        for( final String sample : evalSamples ) {
            if( compSamples.contains( sample ) ) {
                overlappingSamples.add( sample );
            }
        }

        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.ALL_BASES_MODE = true;
        engine = new UnifiedGenotyperEngine(getToolkit(),uac);

        logger.info( "Overlapping samples = " + overlappingSamples );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public CountedData map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        final CountedData counter = new CountedData();

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return counter;
        }

        VariantContext vcEval = tracker.getVariantContext(ref, evalNames.get(0), null, context.getLocation(), false);
        VariantContext vcComp = tracker.getVariantContext(ref, compNames.get(0), null, context.getLocation(), false);

        if( vcEval != null ) { vcEval = vcEval.subContextFromGenotypes( vcEval.getGenotypes(overlappingSamples).values() ); }
        if( vcComp != null ) { vcComp = vcComp.subContextFromGenotypes( vcComp.getGenotypes(overlappingSamples).values() ); }

        VARIANT_STATUS evalStatus;
        VARIANT_STATUS compStatus;

        // First set the variant status variable for both the eval and comp then decide the site's T/F status
        if( vcEval != null && vcEval.isSNP() && vcEval.isPolymorphic() ) {
            if( !vcEval.isFiltered() ) {
                evalStatus = VARIANT_STATUS.CALLED;
            } else {
                evalStatus = VARIANT_STATUS.FILTERED;
            }
        } else {
            evalStatus = VARIANT_STATUS.MISSING;
        }
        if( vcComp != null && vcComp.isSNP() && vcComp.isPolymorphic() ) {
            if( !vcComp.isFiltered() ) {
                compStatus = VARIANT_STATUS.CALLED;
            } else {
                compStatus = VARIANT_STATUS.FILTERED;
            }
        } else {
            compStatus = VARIANT_STATUS.MISSING;
        }

        if( evalStatus == VARIANT_STATUS.CALLED && compStatus == VARIANT_STATUS.CALLED ) { counter.numTP = 1L; }
        else if( evalStatus == VARIANT_STATUS.CALLED && compStatus == VARIANT_STATUS.FILTERED ) {
            counter.numFP = 1L;
            if( printStream!= null ) {
                printStream.println(vcEval.getChr() + ":" + vcEval.getStart() ); // Used to create interval lists of FP variants
            }
        }
        else if( evalStatus == VARIANT_STATUS.CALLED && compStatus == VARIANT_STATUS.MISSING ) {
            VariantCallContext call = engine.runGenotyper(tracker, ref, context);
            if( call != null && call.confidentlyCalled && call.vc != null && call.vc.getType() == VariantContext.Type.NO_VARIATION ) {
                counter.numFP = 1L;
                if( printStream!= null ) {
                    printStream.println(vcEval.getChr() + ":" + vcEval.getStart() ); // Used to create interval lists of FP variants
                }
            }
        }
        else if( evalStatus == VARIANT_STATUS.FILTERED && compStatus == VARIANT_STATUS.CALLED ) { counter.numFilteringFN = 1L; }
        //if( evalStatus == VARIANT_STATUS.FILTERED && compStatus == VARIANT_STATUS.FILTERED ) { counter.numTP = 1L; }
        //if( evalStatus == VARIANT_STATUS.FILTERED && compStatus == VARIANT_STATUS.MISSING ) { counter.numTP = 1L; }
        else if( evalStatus == VARIANT_STATUS.MISSING && compStatus == VARIANT_STATUS.CALLED ) { counter.numCallingFN = 1L; }
        //if( evalStatus == VARIANT_STATUS.MISSING && compStatus == VARIANT_STATUS.FILTERED ) { counter.numTP = 1L; }
        //if( evalStatus == VARIANT_STATUS.MISSING && compStatus == VARIANT_STATUS.MISSING ) { counter.numTP = 1L; }
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

    public CountedData treeReduce( final CountedData sum1, final CountedData sum2) {
        sum2.add(sum1);
        return sum2;
    }

    public void onTraversalDone( CountedData reduceSum ) {
        logger.info("TP = " + reduceSum.numTP);
        logger.info("TN = " + reduceSum.numTN);
        logger.info("FP = " + reduceSum.numFP);
        logger.info("FN = " + (reduceSum.numFilteringFN + reduceSum.numCallingFN) );
        logger.info("  filtering FN = " + reduceSum.numFilteringFN );
        logger.info("  calling   FN = " + reduceSum.numCallingFN );
    }
}
