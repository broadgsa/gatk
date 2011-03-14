/*
 * Copyright (c) 2011 The Broad Institute
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

package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Applies calibrated variant cluster parameters to variant calls to produce an accurate and informative variant quality score.
 *
 * User: rpoplin
 * Date: 3/12/11
 *
 * @help.summary Takes variant calls as .vcf files, learns a Gaussian mixture model over the variant annotations and evaluates the variants
 */

public class ContrastiveRecalibrator extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> implements TreeReducible<ExpandingArrayList<VariantDatum>> {

    public static final String VQS_LOD_KEY = "VQSLOD";

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Output(fullName="recal_file", shortName="recalFile", doc="The output recal file used by ApplyRecalibration", required=true)
    private PrintStream RECAL_FILE;
    @Output(fullName="tranches_file", shortName="tranchesFile", doc="The output tranches file used by ApplyRecalibration", required=true)
    private PrintStream TRANCHES_FILE;

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    //BUGBUG: use VariantRecalibrationArgumentCollection
    @Argument(fullName="use_annotation", shortName="an", doc="The names of the annotations which should used for calculations", required=true)
    private String[] USE_ANNOTATIONS = null;
    @Argument(fullName="FDRtranche", shortName="tranche", doc="The levels of novel false discovery rate (FDR, implied by ti/tv) at which to slice the data. (in percent, that is 1.0 for 1 percent)", required=false)
    private double[] FDR_TRANCHES = new double[]{0.1, 1.0, 10.0, 100.0};

    /////////////////////////////
    // Debug Arguments
    /////////////////////////////
    @Hidden
    @Argument(fullName = "debugFile", shortName = "debugFile", doc = "Print debugging information here", required=false)
    private File DEBUG_FILE = null;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantDataManager dataManager;
    private Set<String> inputNames = new HashSet<String>();
    private VariantRecalibratorEngine engine = new VariantRecalibratorEngine(new VariantRecalibratorArgumentCollection()); //BUGBUG: doesn't do anything with the args at the moment

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        dataManager = new VariantDataManager( new ArrayList<String>(Arrays.asList(USE_ANNOTATIONS)) );

        for( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if( d.getName().toLowerCase().startsWith("input") ) {
                inputNames.add(d.getName());
                logger.info("Found input variant track with name " + d.getName());
            } else if ( d.getName().equals("dbsnp") ) { //BUGBUG: values are hard coded for now until command line tagging up for Rod bindings
                logger.info("Found dbSNP  track: \tKnown = true \tTraining = false \tTruth = false \tPrior = Q10.0");
                dataManager.addTrainingSet(new TrainingSet("dbsnp", true, false, false, 10.0));
            } else if ( d.getName().equals("hapmap") ) {
                logger.info("Found HapMap track: \tKnown = false \tTraining = true \tTruth = true \tPrior = Q12.0");
                dataManager.addTrainingSet(new TrainingSet("hapmap", false, true, true, 12.0));
            } else if ( d.getName().equals("1kg") ) {
                logger.info("Found 1KG    track: \tKnown = false \tTraining = true \tTruth = false \tPrior = Q6.0");
                dataManager.addTrainingSet(new TrainingSet("1kg", false, true, false, 6.0));
            } else if ( d.getName().equals("omni") ) {
                logger.info("Found Omni   track: \tKnown = false \tTraining = true \tTruth = true \tPrior = Q15.0");
                dataManager.addTrainingSet( new TrainingSet("omni", false, true, true, 15.0) );
            } else {
                logger.info("WARNING! Not evaluating ROD binding: " + d.getName());
            }
        }

        if( !dataManager.checkHasTrainingSet() ) {
            throw new UserException.CommandLineException("No training set found! Please provide sets of known polymorphic loci to be used as training data using the dbsnp, hapmap, 1kg, or omni rod bindings.");
        }
        if( !dataManager.checkHasTruthSet() ) {
            throw new UserException.CommandLineException("No truth set found! Please provide sets of known polymorphic loci to be used as training data using the dbsnp, hapmap, 1kg, or omni rod bindings.");
        }
        if( !dataManager.checkHasKnownSet() ) {
            throw new UserException.CommandLineException("No known set found! Please provide sets of known polymorphic loci to be used as training data using the dbsnp, hapmap, 1kg, or omni rod bindings.");
        }

        if( inputNames.size() == 0 ) {
            throw new UserException.BadInput( "No input variant tracks found. Input variant binding names must begin with 'input'." );
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public ExpandingArrayList<VariantDatum> map( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {
        final ExpandingArrayList<VariantDatum> mapList = new ExpandingArrayList<VariantDatum>();

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return mapList;
        }

        VariantContext thisVC = null;
        int maxAC = -1;
        for( final VariantContext vc : tracker.getVariantContexts(ref, inputNames, null, context.getLocation(), true, false) ) {
            if( vc != null ) {
                final int ac = vc.getAttributeAsInt("AC");
                if( ac > maxAC ) {
                    thisVC = vc;
                    maxAC = ac;
                }
            }
        }

        //for( final VariantContext vc : tracker.getVariantContexts(ref, inputNames, null, context.getLocation(), true, false) ) {
            if( thisVC != null ) { // BUGBUG: filtered?
                final VariantDatum datum = new VariantDatum();
                datum.annotations = dataManager.decodeAnnotations( ref.getGenomeLocParser(), thisVC, true ); //BUGBUG: when run with HierarchicalMicroScheduler this is non-deterministic because order of calls depends on load of machine  
                datum.pos = context.getLocation();
                datum.originalQual = thisVC.getPhredScaledQual();
                datum.isTransition = thisVC.isSNP() && ( VariantContextUtils.getSNPSubstitutionType(thisVC).compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0 );
                dataManager.parseTrainingSets( tracker, ref, context, thisVC, datum );
                final double priorFactor = QualityUtils.qualToProb( datum.prior );
                datum.prior = Math.log10( priorFactor ) - Math.log10( 1.0 - priorFactor );
                mapList.add( datum );
            }
        //}

        return mapList;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public ExpandingArrayList<VariantDatum> reduceInit() {
        return new ExpandingArrayList<VariantDatum>();
    }

    public ExpandingArrayList<VariantDatum> reduce( final ExpandingArrayList<VariantDatum> mapValue, final ExpandingArrayList<VariantDatum> reduceSum ) {
        reduceSum.addAll( mapValue );
        return reduceSum;
    }

    public ExpandingArrayList<VariantDatum> treeReduce( final ExpandingArrayList<VariantDatum> lhs, final ExpandingArrayList<VariantDatum> rhs ) {
        rhs.addAll( lhs );
        return rhs;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // on traversal done
    //
    //---------------------------------------------------------------------------------------------------------------

    public void onTraversalDone( final ExpandingArrayList<VariantDatum> reduceSum ) {
        dataManager.setData( reduceSum );
        dataManager.normalizeData();
        engine.evaluateData( dataManager.getData(), engine.generateModel( dataManager.getTrainingData() ), false );
        engine.evaluateData( dataManager.getData(), engine.generateModel( dataManager.selectWorstVariants( 0.07f ) ), true );

        // old tranches stuff
        int nCallsAtTruth = TrancheManager.countCallsAtTruth( dataManager.getData(), Double.NEGATIVE_INFINITY );
        //logger.info(String.format("Truth set size is %d, raw calls at these sites %d, maximum sensitivity of %.2f",
        //        nTruthSites, nCallsAtTruth, (100.0*nCallsAtTruth / Math.max(nTruthSites, nCallsAtTruth))));
        TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric( nCallsAtTruth );
        List<Tranche> tranches = TrancheManager.findTranches( dataManager.getData(), FDR_TRANCHES, metric, DEBUG_FILE ); //BUGBUG: recreated here to match the integration tests
        TRANCHES_FILE.print(Tranche.tranchesString( tranches ));

        logger.info( "Writing out recalibration table..." );
        dataManager.writeOutRecalibrationTable( RECAL_FILE );
    }
}
