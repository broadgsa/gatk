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

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.R.RScriptExecutor;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.io.Resource;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Create a Gaussian mixture model by looking at the annotations values over a high quality subset of the input call set and then evaluate all input variants.
 *
 * <p>
 * This walker is the first pass in a two-stage processing step. This walker is designed to be used in conjunction with ApplyRecalibration walker.
 *
 * <p>
 * The purpose of the variant recalibrator is to assign a well-calibrated probability to each variant call in a call set.
 * One can then create highly accurate call sets by filtering based on this single estimate for the accuracy of each call.
 * The approach taken by variant quality score recalibration is to develop a continuous, covarying estimate of the relationship
 * between SNP call annotations (QD, SB, HaplotypeScore, HRun, for example) and the the probability that a SNP is a true genetic
 * variant versus a sequencing or data processing artifact. This model is determined adaptively based on "true sites" provided
 * as input, typically HapMap 3 sites and those sites found to be polymorphic on the Omni 2.5M SNP chip array. This adaptive
 * error model can then be applied to both known and novel variation discovered in the call set of interest to evaluate the
 * probability that each call is real. The score that gets added to the INFO field of each variant is called the VQSLOD. It is
 * the log odds ratio of being a true variant versus being false under the trained Gaussian mixture model.
 *
 * <p>
 * NOTE: Please see our <a href="http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3">best practices wiki page</a> for our recommendations on which annotations to use for specific project designs.
 *
 * <p>
 * NOTE: In order to create the model reporting plots Rscript needs to be in your environment PATH (this is the scripting version of R, not the interactive version).
 * See <a target="r-project" href="http://www.r-project.org">http://www.r-project.org</a> for more info on how to download and install R.
 *
 * <p>
 * See <a href="http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration">the GATK wiki for a tutorial and example recalibration accuracy plots.</a>
 *
 * <h2>Input</h2>
 * <p>
 * The input raw variants to be recalibrated.
 * <p>
 * Known, truth, and training sets to be used by the algorithm. How these various sets are used is described below.
 *
 * <h2>Output</h2>
 * <p>
 * A recalibration table file in CSV format that is used by the ApplyRecalibration walker.
 * <p>
 * A tranches file which shows various metrics of the recalibration callset as a function of making several slices through the data.
 *
 * <h2>Example</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T VariantRecalibrator \
 *   -R reference/human_g1k_v37.fasta \
 *   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf \
 *   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
 *   -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
 *   -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 dbsnp_132.b37.vcf \
 *   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ \
 *   -recalFile path/to/output.recal \
 *   -tranchesFile path/to/output.tranches \
 *   -rscriptFile path/to/output.plots.R
 * </pre>
 *
 */

@PartitionBy(PartitionType.NONE)
public class VariantRecalibrator extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> implements TreeReducible<ExpandingArrayList<VariantDatum>> {

    public static final String VQS_LOD_KEY = "VQSLOD"; // Log odds ratio of being a true variant versus being false under the trained gaussian mixture model
    public static final String CULPRIT_KEY = "culprit"; // The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out
    private static final String PLOT_TRANCHES_RSCRIPT = "plot_Tranches.R";

    @ArgumentCollection private VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();

    /////////////////////////////
    // Inputs
    /////////////////////////////
    /**
     * These calls should be unfiltered and annotated with the error covariates that are intended to use for modeling.
     */
    @Input(fullName="input", shortName = "input", doc="The raw input variants to be recalibrated", required=true)
    public List<RodBinding<VariantContext>> input;

    /**
     * Any set of VCF files to use as lists of training, truth, or known sites.
     * Training - Input variants which are found to overlap with these training sites are used to build the Gaussian mixture model.
     * Truth - When deciding where to set the cutoff in VQSLOD sensitivity to these truth sites is used.
     * Known - The known / novel status of a variant isn't used by the algorithm itself and is only used for reporting / display purposes.
     * Bad - In addition to using the worst 3% of variants as compared to the Gaussian mixture model, we can also supplement the list with a database of known bad variants.
     */
    @Input(fullName="resource", shortName = "resource", doc="A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm", required=false)
    public List<RodBinding<VariantContext>> resource = Collections.emptyList();

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Output(fullName="recal_file", shortName="recalFile", doc="The output recal file used by ApplyRecalibration", required=true)
    private PrintStream RECAL_FILE;
    @Output(fullName="tranches_file", shortName="tranchesFile", doc="The output tranches file used by ApplyRecalibration", required=true)
    private File TRANCHES_FILE;

    /////////////////////////////
    // Additional Command Line Arguments
    /////////////////////////////
    /**
     * The expected transition / tranversion ratio of true novel variants in your targeted region (whole genome, exome, specific
     * genes), which varies greatly by the CpG and GC content of the region. See expected Ti/Tv ratios section of the GATK best
     * practices wiki documentation for more information. Normal whole genome values are 2.15 and for whole exome 3.2. Note
     * that this parameter is used for display purposes only and isn't used anywhere in the algorithm!
     */
    @Argument(fullName="target_titv", shortName="titv", doc="The expected novel Ti/Tv ratio to use when calculating FDR tranches and for display on the optimization curve output figures. (approx 2.15 for whole genome experiments). ONLY USED FOR PLOTTING PURPOSES!", required=false)
    private double TARGET_TITV = 2.15;

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(fullName="use_annotation", shortName="an", doc="The names of the annotations which should used for calculations", required=true)
    private String[] USE_ANNOTATIONS = null;

    /**
     * Add truth sensitivity slices through the call set at the given values. The default values are 100.0, 99.9, 99.0, and 90.0
     * which will result in 4 estimated tranches in the final call set: the full set of calls (100% sensitivity at the accessible
     * sites in the truth set), a 99.9% truth sensitivity tranche, along with progressively smaller tranches at 99% and 90%.
     */
    @Argument(fullName="TStranche", shortName="tranche", doc="The levels of novel false discovery rate (FDR, implied by ti/tv) at which to slice the data. (in percent, that is 1.0 for 1 percent)", required=false)
    private double[] TS_TRANCHES = new double[] {100.0, 99.9, 99.0, 90.0};
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified the variant recalibrator will use variants even if the specified filter name is marked in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;
    @Output(fullName="rscript_file", shortName="rscriptFile", doc="The output rscript file generated by the VQSR to aid in visualization of the input data and learned model", required=false)
    private File RSCRIPT_FILE = null;
    @Argument(fullName="ts_filter_level", shortName="ts_filter_level", doc="The truth sensitivity level at which to start filtering, used here to indicate filtered variants in the model reporting plots", required=false)
    private double TS_FILTER_LEVEL = 99.0;

    /////////////////////////////
    // Debug Arguments
    /////////////////////////////
    @Hidden
    @Argument(fullName = "trustAllPolymorphic", shortName = "allPoly", doc = "Trust that all the input training sets' unfiltered records contain only polymorphic sites to drastically speed up the computation.", required = false)
    protected Boolean TRUST_ALL_POLYMORPHIC = false;
    //@Hidden
    //@Argument(fullName = "projectConsensus", shortName = "projectConsensus", doc = "Perform 1000G project consensus. This implies an extra prior factor based on the individual participant callsets passed in with consensus=true rod binding tags.", required = false)
    //protected Boolean PERFORM_PROJECT_CONSENSUS = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantDataManager dataManager;
    private PrintStream tranchesStream;
    private final Set<String> ignoreInputFilterSet = new TreeSet<String>();
    private final VariantRecalibratorEngine engine = new VariantRecalibratorEngine( VRAC );

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        dataManager = new VariantDataManager( new ArrayList<String>(Arrays.asList(USE_ANNOTATIONS)), VRAC );

        if (RSCRIPT_FILE != null && !RScriptExecutor.RSCRIPT_EXISTS)
            Utils.warnUser(logger, String.format(
                    "Rscript not found in environment path. %s will be generated but PDF plots will not.",
                    RSCRIPT_FILE));

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( Arrays.asList(IGNORE_INPUT_FILTERS) );
        }

        try {
            tranchesStream = new PrintStream(TRANCHES_FILE);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(TRANCHES_FILE, e);
        }

        for( RodBinding<VariantContext> rod : resource ) {
            dataManager.addTrainingSet( new TrainingSet( rod ) );
        }

        if( !dataManager.checkHasTrainingSet() ) {
            throw new UserException.CommandLineException( "No training set found! Please provide sets of known polymorphic loci marked with the training=true ROD binding tag. For example, -B:hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }
        if( !dataManager.checkHasTruthSet() ) {
            throw new UserException.CommandLineException( "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true ROD binding tag. For example, -B:hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
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

        for( final VariantContext vc : tracker.getValues(input, context.getLocation()) ) {
            if( vc != null && ( vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()) ) ) {
                if( checkRecalibrationMode( vc, VRAC.MODE ) ) {
                    final VariantDatum datum = new VariantDatum();

                    // Populate the datum with lots of fields from the VariantContext, unfortunately the VC is too big so we just pull in only the things we absolutely need.
                    dataManager.decodeAnnotations( datum, vc, true ); //BUGBUG: when run with HierarchicalMicroScheduler this is non-deterministic because order of calls depends on load of machine
                    datum.contig = vc.getChr();
                    datum.start = vc.getStart();
                    datum.stop = vc.getEnd();
                    datum.originalQual = vc.getPhredScaledQual();
                    datum.isSNP = vc.isSNP() && vc.isBiallelic();
                    datum.isTransition = datum.isSNP && VariantContextUtils.isTransition(vc);

                    // Loop through the training data sets and if they overlap this loci then update the prior and training status appropriately
                    dataManager.parseTrainingSets( tracker, context.getLocation(), vc, datum, TRUST_ALL_POLYMORPHIC );
                    double priorFactor = QualityUtils.qualToProb( datum.prior );
                    //if( PERFORM_PROJECT_CONSENSUS ) { // BUGBUG: need to resurrect this functionality?
                    //    final double consensusPrior = QualityUtils.qualToProb( 1.0 + 5.0 * datum.consensusCount );
                    //    priorFactor = 1.0 - ((1.0 - priorFactor) * (1.0 - consensusPrior));
                    //}
                    datum.prior = Math.log10( priorFactor ) - Math.log10( 1.0 - priorFactor );

                    mapList.add( datum );
                }
            }
        }

        return mapList;
    }

    public static boolean checkRecalibrationMode( final VariantContext vc, final VariantRecalibratorArgumentCollection.Mode mode ) {
        return mode == VariantRecalibratorArgumentCollection.Mode.BOTH ||
                (mode == VariantRecalibratorArgumentCollection.Mode.SNP && vc.isSNP()) ||
	            (mode == VariantRecalibratorArgumentCollection.Mode.INDEL && (vc.isIndel() || vc.isMixed()));
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
        dataManager.normalizeData(); // Each data point is now (x - mean) / standard deviation

        // Generate the positive model using the training data and evaluate each variant
        final GaussianMixtureModel goodModel = engine.generateModel( dataManager.getTrainingData() );
        engine.evaluateData( dataManager.getData(), goodModel, false );

        // Generate the negative model using the worst performing data and evaluate each variant contrastively
        final GaussianMixtureModel badModel = engine.generateModel( dataManager.selectWorstVariants( VRAC.PERCENT_BAD_VARIANTS, VRAC.MIN_NUM_BAD_VARIANTS ) );
        engine.evaluateData( dataManager.getData(), badModel, true );
        engine.calculateWorstPerformingAnnotation( dataManager.getData(), goodModel, badModel );

        // Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
        final int nCallsAtTruth = TrancheManager.countCallsAtTruth( dataManager.getData(), Double.NEGATIVE_INFINITY );
        final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric( nCallsAtTruth );
        final List<Tranche> tranches = TrancheManager.findTranches( dataManager.getData(), TS_TRANCHES, metric );
        tranchesStream.print(Tranche.tranchesString( tranches ));

        // Find the filtering lodCutoff for display on the model PDFs. Red variants are those which were below the cutoff and filtered out of the final callset.
        double lodCutoff = 0.0;
        for( final Tranche tranche : tranches ) {
            if( MathUtils.compareDoubles(tranche.ts, TS_FILTER_LEVEL, 0.0001)==0 ) {
                lodCutoff = tranche.minVQSLod;
            }
        }

        logger.info( "Writing out recalibration table..." );
        dataManager.writeOutRecalibrationTable( RECAL_FILE );
        if( RSCRIPT_FILE != null ) {
            logger.info( "Writing out visualization Rscript file...");
            createVisualizationScript( dataManager.getRandomDataForPlotting( 6000 ), goodModel, badModel, lodCutoff );
        }

        // Execute the RScript command to plot the table of truth values
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(PLOT_TRANCHES_RSCRIPT, VariantRecalibrator.class));
        executor.addArgs(TRANCHES_FILE.getAbsoluteFile(), TARGET_TITV);
        // Print out the command line to make it clear to the user what is being executed and how one might modify it
        logger.info("Executing: " + executor.getApproximateCommandLine());
        executor.exec();
    }

    private void createVisualizationScript( final ExpandingArrayList<VariantDatum> randomData, final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel, final double lodCutoff ) {
        PrintStream stream;
        try {
            stream = new PrintStream(RSCRIPT_FILE);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(RSCRIPT_FILE, e);
        }

        // We make extensive use of the ggplot2 R library: http://had.co.nz/ggplot2/
        stream.println("library(ggplot2)");
        // For compactPDF in R 2.13+
        stream.println("library(tools)");

        createArrangeFunction( stream );

        stream.println("outputPDF <- \"" + RSCRIPT_FILE + ".pdf\"");
        stream.println("pdf(outputPDF)"); // Unfortunately this is a huge pdf file, BUGBUG: need to work on reducing the file size

        for(int iii = 0; iii < USE_ANNOTATIONS.length; iii++) {
            for( int jjj = iii + 1; jjj < USE_ANNOTATIONS.length; jjj++) {
                logger.info( "Building " + USE_ANNOTATIONS[iii] + " x " + USE_ANNOTATIONS[jjj] + " plot...");

                final ExpandingArrayList<VariantDatum> fakeData = new ExpandingArrayList<VariantDatum>();
                double minAnn1 = 100.0, maxAnn1 = -100.0, minAnn2 = 100.0, maxAnn2 = -100.0;
                for( final VariantDatum datum : randomData ) {
                    minAnn1 = Math.min(minAnn1, datum.annotations[iii]);
                    maxAnn1 = Math.max(maxAnn1, datum.annotations[iii]);
                    minAnn2 = Math.min(minAnn2, datum.annotations[jjj]);
                    maxAnn2 = Math.max(maxAnn2, datum.annotations[jjj]);
                }
                // Create a fake set of data which spans the full extent of these two annotation dimensions in order to calculate the model PDF projected to 2D
                for(double ann1 = minAnn1; ann1 <= maxAnn1; ann1+=0.1) {
                    for(double ann2 = minAnn2; ann2 <= maxAnn2; ann2+=0.1) {
                        final VariantDatum datum = new VariantDatum();
                        datum.prior = 0.0;
                        datum.annotations = new double[randomData.get(0).annotations.length];
                        datum.isNull = new boolean[randomData.get(0).annotations.length];
                        for(int ann=0; ann< datum.annotations.length; ann++) {
                            datum.annotations[ann] = 0.0;
                            datum.isNull[ann] = true;
                        }
                        datum.annotations[iii] = ann1;
                        datum.annotations[jjj] = ann2;
                        datum.isNull[iii] = false;
                        datum.isNull[jjj] = false;
                        fakeData.add(datum);
                    }
                }

                engine.evaluateData( fakeData, goodModel, false );
                engine.evaluateData( fakeData, badModel, true );

                stream.print("surface <- c(");
                for( final VariantDatum datum : fakeData ) {
                    stream.print(String.format("%.3f, %.3f, %.3f, ", datum.annotations[iii], datum.annotations[jjj], Math.min(4.0, Math.max(-4.0, datum.lod))));
                }
                stream.println("NA,NA,NA)");
                stream.println("s <- matrix(surface,ncol=3,byrow=T)");

                stream.print("data <- c(");
                for( final VariantDatum datum : randomData ) {
                    stream.print(String.format("%.3f, %.3f, %.3f, %d, %d,", datum.annotations[iii], datum.annotations[jjj], (datum.lod < lodCutoff ? -1.0 : 1.0),
                            (datum.atAntiTrainingSite ? -1 : (datum.atTrainingSite ? 1 : 0)), (datum.isKnown ? 1 : -1)));
                }
                stream.println("NA,NA,NA,NA,1)");
                stream.println("d <- matrix(data,ncol=5,byrow=T)");

                final String surfaceFrame = "sf." + USE_ANNOTATIONS[iii] + "." + USE_ANNOTATIONS[jjj];
                final String dataFrame = "df." + USE_ANNOTATIONS[iii] + "." + USE_ANNOTATIONS[jjj];

                stream.println(surfaceFrame + " <- data.frame(x=s[,1], y=s[,2], lod=s[,3])");
                stream.println(dataFrame + " <- data.frame(x=d[,1], y=d[,2], retained=d[,3], training=d[,4], novelty=d[,5])");
                stream.println("dummyData <- " + dataFrame + "[1,]");
                stream.println("dummyData$x <- NaN");
                stream.println("dummyData$y <- NaN");
                stream.println("p <- ggplot(data=" + surfaceFrame + ", aes(x=x, y=y)) + opts(panel.background = theme_rect(colour = NA), panel.grid.minor = theme_line(colour = NA), panel.grid.major = theme_line(colour = NA))");
                stream.println("p1 = p + opts(title=\"model PDF\") + labs(x=\""+ USE_ANNOTATIONS[iii] +"\", y=\""+ USE_ANNOTATIONS[jjj] +"\") + geom_tile(aes(fill = lod)) + scale_fill_gradient(high=\"green\", low=\"red\")");
                stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=retained, alpha=I(1/7),legend=FALSE) + opts(panel.background = theme_rect(colour = NA), panel.grid.minor = theme_line(colour = NA), panel.grid.major = theme_line(colour = NA))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=retained),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p2 = p + q + labs(x=\""+ USE_ANNOTATIONS[iii] +"\", y=\""+ USE_ANNOTATIONS[jjj] +"\") + scale_colour_gradient(name=\"outcome\", high=\"black\", low=\"red\",breaks=c(-1,1),labels=c(\"filtered\",\"retained\"))");
                stream.println("p <- qplot(x,y,data="+ dataFrame + "["+dataFrame+"$training != 0,], color=training, alpha=I(1/7)) + opts(panel.background = theme_rect(colour = NA), panel.grid.minor = theme_line(colour = NA), panel.grid.major = theme_line(colour = NA))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=training),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p3 = p + q + labs(x=\""+ USE_ANNOTATIONS[iii] +"\", y=\""+ USE_ANNOTATIONS[jjj] +"\") + scale_colour_gradient(high=\"green\", low=\"purple\",breaks=c(-1,1), labels=c(\"neg\", \"pos\"))");
                stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=novelty, alpha=I(1/7)) + opts(panel.background = theme_rect(colour = NA), panel.grid.minor = theme_line(colour = NA), panel.grid.major = theme_line(colour = NA))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=novelty),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p4 = p + q + labs(x=\""+ USE_ANNOTATIONS[iii] +"\", y=\""+ USE_ANNOTATIONS[jjj] +"\") + scale_colour_gradient(name=\"novelty\", high=\"blue\", low=\"red\",breaks=c(-1,1), labels=c(\"novel\",\"known\"))");
                stream.println("arrange(p1, p2, p3, p4, ncol=2)");
            }
        }
        stream.println("dev.off()");

        stream.println("if (exists(\"compactPDF\")) {");
        stream.println("compactPDF(ouputPDF)");
        stream.println("}");

        stream.close();

        // Execute Rscript command to generate the clustering plots
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(RSCRIPT_FILE);
        logger.info("Executing: " + executor.getApproximateCommandLine());
        executor.exec();
     }

    // The Arrange function is how we place the 4 model plots on one page
    // from http://gettinggeneticsdone.blogspot.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html
    private void createArrangeFunction( final PrintStream stream ) {
        stream.println("vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)");
        stream.println("arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {");
        stream.println("dots <- list(...)");
        stream.println("n <- length(dots)");
        stream.println("if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}");
        stream.println("if(is.null(nrow)) { nrow = ceiling(n/ncol)}");
        stream.println("if(is.null(ncol)) { ncol = ceiling(n/nrow)}");
        stream.println("grid.newpage()");
        stream.println("pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )");
        stream.println("ii.p <- 1");
        stream.println("for(ii.row in seq(1, nrow)){");
        stream.println("ii.table.row <- ii.row ");
        stream.println("if(as.table) {ii.table.row <- nrow - ii.table.row + 1}");
        stream.println("for(ii.col in seq(1, ncol)){");
        stream.println("ii.table <- ii.p");
        stream.println("if(ii.p > n) break");
        stream.println("print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))");
        stream.println("ii.p <- ii.p + 1");
        stream.println("}");
        stream.println("}");
        stream.println("}");
    }
}
