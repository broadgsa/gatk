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
    @Argument(fullName="TStranche", shortName="tranche", doc="The levels of novel false discovery rate (FDR, implied by ti/tv) at which to slice the data. (in percent, that is 1.0 for 1 percent)", required=false)
    private double[] TS_TRANCHES = new double[] {100.0, 99.9, 99.0, 90.0};
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified the optimizer will use variants even if the specified filter name is marked in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;

    /////////////////////////////
    // Debug Arguments
    /////////////////////////////
    @Hidden
    @Argument(fullName = "debugFile", shortName = "debugFile", doc = "Print debugging information here", required=false)
    private File DEBUG_FILE = null;
    @Hidden
    @Argument(fullName = "trustAllPolymorphic", shortName = "allPoly", doc = "Trust that all the input training sets' unfiltered records contain only polymorphic sites to drastically speed up the computation.", required = false)
    protected Boolean TRUST_ALL_POLYMORPHIC = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantDataManager dataManager;
    private final Set<String> ignoreInputFilterSet = new TreeSet<String>();
    private final Set<String> inputNames = new HashSet<String>();
    private final VariantRecalibratorEngine engine = new VariantRecalibratorEngine(new VariantRecalibratorArgumentCollection()); //BUGBUG: doesn't do anything with the args at the moment

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        dataManager = new VariantDataManager( new ArrayList<String>(Arrays.asList(USE_ANNOTATIONS)) );

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( Arrays.asList(IGNORE_INPUT_FILTERS) );
        }

        for( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if( d.getName().toLowerCase().startsWith("input") ) {
                inputNames.add(d.getName());
                logger.info( "Found input variant track with name " + d.getName() );
            } else {
                dataManager.addTrainingSet( new TrainingSet(d.getName(), d.getTags()) );
            }
        }

        if( !dataManager.checkHasTrainingSet() ) {
            throw new UserException.CommandLineException( "No training set found! Please provide sets of known polymorphic loci marked with the training=true ROD binding tag. For example, -B:hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }
        if( !dataManager.checkHasTruthSet() ) {
            throw new UserException.CommandLineException( "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true ROD binding tag. For example, -B:hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }
        if( !dataManager.checkHasKnownSet() ) {
            throw new UserException.CommandLineException( "No known set found! Please provide sets of known polymorphic loci marked with the known=true ROD binding tag. For example, -B:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 dbsnpFile.vcf" );
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

        for( final VariantContext vc : tracker.getVariantContexts(ref, inputNames, null, context.getLocation(), true, false) ) {
            if( vc != null && ( vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()) ) ) {
                final VariantDatum datum = new VariantDatum();
                datum.annotations = dataManager.decodeAnnotations( ref.getGenomeLocParser(), vc, true ); //BUGBUG: when run with HierarchicalMicroScheduler this is non-deterministic because order of calls depends on load of machine
                datum.pos = context.getLocation();
                datum.originalQual = vc.getPhredScaledQual();
                datum.isTransition = vc.isSNP() && vc.isBiallelic() && VariantContextUtils.isTransition(vc);
                dataManager.parseTrainingSets( tracker, ref, context, vc, datum, TRUST_ALL_POLYMORPHIC );
                final double priorFactor = QualityUtils.qualToProb( datum.prior );
                datum.prior = Math.log10( priorFactor ) - Math.log10( 1.0 - priorFactor );
                mapList.add( datum );
            }
        }

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
        List<Tranche> tranches = TrancheManager.findTranches( dataManager.getData(), TS_TRANCHES, metric, DEBUG_FILE ); //BUGBUG: recreated here to match the integration tests
        TRANCHES_FILE.print(Tranche.tranchesString( tranches ));

        logger.info( "Writing out recalibration table..." );
        dataManager.writeOutRecalibrationTable( RECAL_FILE );
    }
}
