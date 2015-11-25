/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.varianteval;

import com.google.java.contract.Requires;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.SAMSequenceRecord;
import oracle.jrockit.jfr.StringConstantPool;
import org.apache.log4j.Logger;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.engine.samples.Trio;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.varianteval.evaluators.*;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.IntervalStratification;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.VariantEvalUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)
 *
 * <p>
 * Given a variant callset, it is common to calculate various quality control metrics. These metrics include the number of
 * raw or filtered SNP counts; ratio of transition mutations to transversions; concordance of a particular sample's calls
 * to a genotyping chip; number of singletons per sample; etc. Furthermore, it is often useful to stratify these metrics
 * by various criteria like functional class (missense, nonsense, silent), whether the site is CpG site, the amino acid
 * degeneracy of the site, etc. VariantEval facilitates these calculations in two ways: by providing several built-in
 * evaluation and stratification modules, and by providing a framework that permits the easy development of new evaluation
 * and stratification modules.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more variant sets to evaluate plus any number of comparison sets.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Evaluation tables detailing the results of the eval modules which were applied.
 * For example:
 * </p>
 * <pre>
 * output.eval.grp:
 * ##:GATKReport.v0.1 CountVariants : Counts different classes of variants in the sample
 * CountVariants  CompRod   CpG      EvalRod  JexlExpression  Novelty  nProcessedLoci  nCalledLoci  nRefLoci  nVariantLoci  variantRate ...
 * CountVariants  dbsnp     CpG      eval     none            all      65900028        135770       0         135770        0.00206024  ...
 * CountVariants  dbsnp     CpG      eval     none            known    65900028        47068        0         47068         0.00071423  ...
 * CountVariants  dbsnp     CpG      eval     none            novel    65900028        88702        0         88702         0.00134601  ...
 * CountVariants  dbsnp     all      eval     none            all      65900028        330818       0         330818        0.00502000  ...
 * CountVariants  dbsnp     all      eval     none            known    65900028        120685       0         120685        0.00183133  ...
 * CountVariants  dbsnp     all      eval     none            novel    65900028        210133       0         210133        0.00318866  ...
 * CountVariants  dbsnp     non_CpG  eval     none            all      65900028        195048       0         195048        0.00295976  ...
 * CountVariants  dbsnp     non_CpG  eval     none            known    65900028        73617        0         73617         0.00111710  ...
 * CountVariants  dbsnp     non_CpG  eval     none            novel    65900028        121431       0         121431        0.00184265  ...
 * ...
 * </pre>
 * </p>
 *
 * <h3>Usage examples</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T VariantEval \
 *   -R reference.fasta \
 *   -o output.eval.grp \
 *   --eval:set1 set1.vcf \
 *   --eval:set2 set2.vcf \
 *   [--comp comp.vcf]
 * </pre>
 *
 * Count Mendelian violations for each family in a callset with multiple families (and provided pedigree)
 * <pre>
 * Java -jar GenomeAnalysisTK.jar \
 *   -T VariantEval \
 *   -R reference.fasta \
 *   -o output.MVs.byFamily.table \
 *   --eval multiFamilyCallset.vcf \
 *   -noEV -noST \
 *   -ST Family \
 *   -EV MendelianViolationEvaluator
 * </pre>
 *
 * <h3>Caveat</h3>
 *
 * <p>Some stratifications and evaluators are incompatible with each other due to their respective memory requirements,
 * such as AlleleCount and VariantSummary, or Sample and VariantSummary. If you specify such a combination, the program
 * will output an error message and ask you to disable one of these options. We do not currently provide an exhaustive
 * list of incompatible combinations, so we recommend trying out combinations that you are interested in on a dummy
 * command line, to rapidly ascertain whether it will work or not.</p>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-50, stop=50))
@PartitionBy(PartitionType.NONE)
public class VariantEval extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    public static final String IS_SINGLETON_KEY = "ISSINGLETON";

    @Output
    protected PrintStream out;

    /**
     * The variant file(s) to evaluate.
     */
    @Input(fullName="eval", shortName = "eval", doc="Input evaluation file(s)", required=true)
    public List<RodBinding<VariantContext>> evals;

    /**
     * The variant file(s) to compare against.
     */
    @Input(fullName="comp", shortName = "comp", doc="Input comparison file(s)", required=false)
    public List<RodBinding<VariantContext>> compsProvided = Collections.emptyList();
    private List<RodBinding<VariantContext>> comps = new ArrayList<RodBinding<VariantContext>>();

    /**
     * dbSNP comparison VCF.  By default, the dbSNP file is used to specify the set of "known" variants.
     * Other sets can be specified with the -knownName (--known_names) argument.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * Some analyses want to count overlap not with dbSNP (which is in general very open) but
     * actually want to itemize their overlap specifically with a set of gold standard sites
     * such as HapMap, OMNI, or the gold standard indels.  This argument provides a mechanism
     * for communicating which file to use
     */
    @Input(fullName="goldStandard", shortName = "gold", doc="Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison", required=false)
    public RodBinding<VariantContext> goldStandard = null;

    /**
     * Note that the --list argument requires a fully resolved and correct command-line to work.
     */
    @Argument(fullName="list", shortName="ls", doc="List the available eval modules and exit", required=false)
    protected Boolean LIST = false;

    // Partitioning the data arguments
    @Argument(shortName="select", doc="One or more stratifications to use when evaluating the data", required=false)
    protected ArrayList<String> SELECT_EXPS = new ArrayList<String>();

    @Argument(shortName="selectName", doc="Names to use for the list of stratifications (must be a 1-to-1 mapping)", required=false)
    protected ArrayList<String> SELECT_NAMES = new ArrayList<String>();

    @Argument(fullName="sample", shortName="sn", doc="Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context", required=false)
    protected Set<String> SAMPLE_EXPRESSIONS;

    /**
     * List of rod tracks to be used for specifying "known" variants other than dbSNP.
     */
    @Argument(shortName="knownName", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected HashSet<String> KNOWN_NAMES = new HashSet<String>();
    List<RodBinding<VariantContext>> knowns = new ArrayList<RodBinding<VariantContext>>();

    // Stratification arguments
    @Argument(fullName="stratificationModule", shortName="ST", doc="One or more specific stratification modules to apply to the eval track(s) (in addition to the standard stratifications, unless -noS is specified)", required=false)
    protected String[] STRATIFICATIONS_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardStratifications", shortName="noST", doc="Do not use the standard stratification modules by default (instead, only those that are specified with the -S option)", required=false)
    protected Boolean NO_STANDARD_STRATIFICATIONS = false;

    /**
     * See the -list argument to view available modules.
     */
    @Argument(fullName="evalModule", shortName="EV", doc="One or more specific eval modules to apply to the eval track(s) (in addition to the standard modules, unless -noEV is specified)", required=false)
    protected String[] MODULES_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardModules", shortName="noEV", doc="Do not use the standard modules by default (instead, only those that are specified with the -EV option)", required=false)
    protected Boolean NO_STANDARD_MODULES = false;

    @Argument(fullName="minPhaseQuality", shortName="mpq", doc="Minimum phasing quality", required=false)
    protected double MIN_PHASE_QUALITY = 10.0;

    @Argument(shortName="mvq", fullName="mendelianViolationQualThreshold", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation. Default is 50.", required=false)
    protected double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 50;

    @Argument(shortName="ploidy", fullName="samplePloidy", doc="Per-sample ploidy (number of chromosomes per sample)", required=false)
    protected int ploidy = GATKVariantContextUtils.DEFAULT_PLOIDY;

    @Argument(fullName="ancestralAlignments", shortName="aa", doc="Fasta file with ancestral alleles", required=false)
    private File ancestralAlignmentsFile = null;

    @Argument(fullName="requireStrictAlleleMatch", shortName="strict", doc="If provided only comp and eval tracks with exactly matching reference and alternate alleles will be counted as overlapping", required=false)
    private boolean requireStrictAlleleMatch = false;

    @Argument(fullName="keepAC0", shortName="keepAC0", doc="If provided, modules that track polymorphic sites will not require that a site have AC > 0 when the input eval has genotypes", required=false)
    private boolean keepSitesWithAC0 = false;

    @Hidden
    @Argument(fullName="numSamples", shortName="numSamples", doc="If provided, modules that track polymorphic sites will not require that a site have AC > 0 when the input eval has genotypes", required=false)
    private int numSamplesFromArgument = 0;

    /**
     * If true, VariantEval will treat -eval 1 -eval 2 as separate tracks from the same underlying
     * variant set, and evaluate the union of the results.  Useful when you want to do -eval chr1.vcf -eval chr2.vcf etc.
     */
    @Argument(fullName="mergeEvals", shortName="mergeEvals", doc="If provided, all -eval tracks will be merged into a single eval track", required=false)
    public boolean mergeEvals = false;

    /**
     * File containing tribble-readable features for the IntervalStratificiation
     */
    @Input(fullName="stratIntervals", shortName="stratIntervals", doc="File containing tribble-readable features for the IntervalStratificiation", required=false)
    public IntervalBinding<Feature> intervalsFile = null;

    /**
     * File containing tribble-readable features containing known CNVs.  For use with VariantSummary table.
     */
    @Input(fullName="knownCNVs", shortName="knownCNVs", doc="File containing tribble-readable features describing a known list of copy number variants", required=false)
    public IntervalBinding<Feature> knownCNVsFile = null;
    Map<String, IntervalTree<GenomeLoc>> knownCNVsByContig = Collections.emptyMap();

    // Variables
    private Set<SortableJexlVCMatchExp> jexlExpressions = new TreeSet<SortableJexlVCMatchExp>();

    private boolean isSubsettingSamples;
    private Set<String> sampleNamesForEvaluation = new LinkedHashSet<String>();
    private Set<String> familyNamesForEvaluation = new LinkedHashSet<String>();
    private Set<String> sampleNamesForStratification = new LinkedHashSet<String>();
    private Set<String> familyNamesForStratification = new LinkedHashSet<String>();

    // important stratifications
    private boolean byFilterIsEnabled = false;
    private boolean perSampleIsEnabled = false;
    private boolean perFamilyIsEnabled = false;

    // Public constants
    final private static String ALL_SAMPLE_NAME = "all";
    final private static String ALL_FAMILY_NAME = "all";

    // the number of processed bp for this walker
    long nProcessedLoci = 0;

    // Utility class
    private final VariantEvalUtils variantEvalUtils = new VariantEvalUtils(this);

    // Ancestral alignments
    private IndexedFastaSequenceFile ancestralAlignments = null;

    // The set of all possible evaluation contexts
    StratificationManager<VariantStratifier, EvaluationContext> stratManager;
    //Set<DynamicStratification> dynamicStratifications = Collections.emptySet();

    /**
     * Initialize the stratifications, evaluations, evaluation contexts, and reporting object
     */
    public void initialize() {
        // Just list the modules, and exit quickly.
        if (LIST) { variantEvalUtils.listModulesAndExit(); }

        // maintain the full list of comps
        comps.addAll(compsProvided);
        if ( dbsnp.dbsnp.isBound() ) {
            comps.add(dbsnp.dbsnp);
            knowns.add(dbsnp.dbsnp);
        }

        // Add a dummy comp track if none exists
        if ( comps.size() == 0 )
            comps.add(new RodBinding<VariantContext>(VariantContext.class, "none", "UNBOUND", "", new Tags()));

        // Set up set of additional knowns
        for ( RodBinding<VariantContext> compRod : comps ) {
            if ( KNOWN_NAMES.contains(compRod.getName()) )
                knowns.add(compRod);
        }

        // Now that we have all the rods categorized, determine the sample list from the eval rods.
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), evals);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        // Load the sample list, using an intermediate tree set to sort the samples
        final Set<String> allSampleNames = SampleUtils.getSamplesFromCommandLineInput(vcfSamples);
        sampleNamesForEvaluation.addAll(new TreeSet<String>(SampleUtils.getSamplesFromCommandLineInput(vcfSamples, SAMPLE_EXPRESSIONS)));
        isSubsettingSamples = ! sampleNamesForEvaluation.containsAll(allSampleNames);
        familyNamesForEvaluation.addAll(getSampleDB().getFamilyIDs());

        //If stratifying by sample name, assign a stratification for each sample we're evaluating (based on commandline args)...
        if (Arrays.asList(STRATIFICATIONS_TO_USE).contains("Sample") ) {
            sampleNamesForStratification.addAll(sampleNamesForEvaluation);
        }
        //...and also a stratification for the sum over all samples
        sampleNamesForStratification.add(ALL_SAMPLE_NAME);

        //If stratifying by sample name, assign a stratification for each family...
        if ( Arrays.asList(STRATIFICATIONS_TO_USE).contains("Family") ) {
            familyNamesForStratification.addAll(familyNamesForEvaluation);
        }
        //...and also a stratification for the sum over all families
        familyNamesForStratification.add(ALL_FAMILY_NAME);

        // Initialize select expressions
        for (VariantContextUtils.JexlVCMatchExp jexl : VariantContextUtils.initializeMatchExps(SELECT_NAMES, SELECT_EXPS)) {
            SortableJexlVCMatchExp sjexl = new SortableJexlVCMatchExp(jexl.name, jexl.exp);
            jexlExpressions.add(sjexl);
        }

        // Initialize the set of stratifications and evaluations to use
        // The list of stratifiers and evaluators to use
        final List<VariantStratifier> stratificationObjects = variantEvalUtils.initializeStratificationObjects(NO_STANDARD_STRATIFICATIONS, STRATIFICATIONS_TO_USE);
        final Set<Class<? extends VariantEvaluator>> evaluationClasses = variantEvalUtils.initializeEvaluationObjects(NO_STANDARD_MODULES, MODULES_TO_USE);

        checkForIncompatibleEvaluatorsAndStratifiers(stratificationObjects, evaluationClasses);

        for ( VariantStratifier vs : stratificationObjects ) {
            if ( vs.getName().equals("Filter") )
                byFilterIsEnabled = true;
            else if ( vs.getName().equals("Sample") )
                perSampleIsEnabled = true;
            else if ( vs.getName().equals("Family"))
                perFamilyIsEnabled = true;
        }

        if (perSampleIsEnabled && perFamilyIsEnabled)
            throw new UserException.BadArgumentValue("ST", "Variants cannot be stratified by sample and family at the same time");

        if (perFamilyIsEnabled && getSampleDB().getTrios().isEmpty())
            throw new UserException.BadArgumentValue("ST", "Cannot stratify by family without *.ped file");


        if ( intervalsFile != null ) {
            boolean fail = true;
            for ( final VariantStratifier vs : stratificationObjects ) {
                if ( vs.getClass().equals(IntervalStratification.class) )
                    fail = false;
            }
            if ( fail )
                throw new UserException.BadArgumentValue("ST", "stratIntervals argument provided but -ST IntervalStratification not provided");
        }

        // Initialize the evaluation contexts
        createStratificationStates(stratificationObjects, evaluationClasses);

        // Load ancestral alignments
        if (ancestralAlignmentsFile != null) {
            try {
                ancestralAlignments = new IndexedFastaSequenceFile(ancestralAlignmentsFile);
            } catch (FileNotFoundException e) {
                throw new ReviewedGATKException(String.format("The ancestral alignments file, '%s', could not be found", ancestralAlignmentsFile.getAbsolutePath()));
            }
        }

        // initialize CNVs
        if ( knownCNVsFile != null ) {
            knownCNVsByContig = createIntervalTreeByContig(knownCNVsFile);
        }
    }

    final void checkForIncompatibleEvaluatorsAndStratifiers( final List<VariantStratifier> stratificationObjects,
                                                             Set<Class<? extends VariantEvaluator>> evaluationClasses) {
        for ( final VariantStratifier vs : stratificationObjects ) {
            for ( Class<? extends VariantEvaluator> ec : evaluationClasses )
                if ( vs.getIncompatibleEvaluators().contains(ec) )
                    throw new UserException.BadArgumentValue("ST and ET", 
                            "The selected stratification " + vs.getName() + 
                                    " and evaluator " + ec.getSimpleName() +
                                    " are incompatible due to combinatorial memory requirements." +
                                    " Please disable one");
        }
    }
    
    final void createStratificationStates(final List<VariantStratifier> stratificationObjects, final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        final List<VariantStratifier> strats = new ArrayList<VariantStratifier>(stratificationObjects);
        stratManager = new StratificationManager<VariantStratifier, EvaluationContext>(strats);

        logger.info("Creating " + stratManager.size() + " combinatorial stratification states");
        for ( int i = 0; i < stratManager.size(); i++ ) {
            EvaluationContext ec = new EvaluationContext(this, evaluationObjects);
            stratManager.set(i, ec);
        }
    }    
    
    public final Map<String, IntervalTree<GenomeLoc>> createIntervalTreeByContig(final IntervalBinding<Feature> intervals) {
        final Map<String, IntervalTree<GenomeLoc>> byContig = new HashMap<String, IntervalTree<GenomeLoc>>();

        final List<GenomeLoc> locs = intervals.getIntervals(getToolkit().getGenomeLocParser());

        // set up the map from contig -> interval tree
        for ( final String contig : getContigNames() )
            byContig.put(contig, new IntervalTree<GenomeLoc>());

        for ( final GenomeLoc loc : locs ) {
            byContig.get(loc.getContig()).put(loc.getStart(), loc.getStop(), loc);
        }

        return byContig;
    }

    /**
     * Collect relevant information from each variant in the supplied VCFs
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        // we track the processed bp and expose this for modules instead of wasting CPU power on calculating
        // the same thing over and over in evals that want the processed bp
        synchronized (this) {
            nProcessedLoci += context.getSkippedBases() + (ref == null ? 0 : 1);
        }

        if (tracker != null) {
            String aastr = (ancestralAlignments == null) ? null : new String(ancestralAlignments.getSubsequenceAt(ref.getLocus().getContig(), ref.getLocus().getStart(), ref.getLocus().getStop()).getBases());

//            // update the dynamic stratifications
//            for (final VariantContext vc : tracker.getValues(evals, ref.getLocus())) {
//                // don't worry -- DynamicStratification only work with one eval object
//                for ( final DynamicStratification ds :  dynamicStratifications ) {
//                    ds.update(vc);
//                }
//            }

            //      --------- track ---------           sample  - VariantContexts -
            HashMap<RodBinding<VariantContext>, HashMap<String, Collection<VariantContext>>> evalVCs = variantEvalUtils.bindVariantContexts(tracker, ref, evals, byFilterIsEnabled, true, perSampleIsEnabled, perFamilyIsEnabled, mergeEvals);
            HashMap<RodBinding<VariantContext>, HashMap<String, Collection<VariantContext>>> compVCs = variantEvalUtils.bindVariantContexts(tracker, ref, comps, byFilterIsEnabled, false, false, false, false);

            // for each eval track
            for ( final RodBinding<VariantContext> evalRod : evals ) {
                final Map<String, Collection<VariantContext>> emptyEvalMap = Collections.emptyMap();
                final Map<String, Collection<VariantContext>> evalSet = evalVCs.containsKey(evalRod) ? evalVCs.get(evalRod) : emptyEvalMap;

                Set<String> statificationLevels;

                // for each sample stratifier
                if (perFamilyIsEnabled)
                    statificationLevels = familyNamesForStratification;
                else
                    statificationLevels = sampleNamesForStratification;
                for ( final String stratLevelName : statificationLevels ) {
                    Collection<VariantContext> evalSetBySample = evalSet.get(stratLevelName);

                    if ( evalSetBySample == null ) {
                        evalSetBySample = new HashSet<VariantContext>(1);
                        evalSetBySample.add(null);
                    }

                    // for each eval in the track
                    for ( VariantContext eval : evalSetBySample ) {
                        // deal with ancestral alleles if requested
                        if ( eval != null && aastr != null ) {
                            eval = new VariantContextBuilder(eval).attribute("ANCESTRALALLELE", aastr).make();
                        }

                        // for each comp track
                        for ( final RodBinding<VariantContext> compRod : comps ) {
                            // no sample stratification for comps
                            final HashMap<String, Collection<VariantContext>> compSetHash = compVCs.get(compRod);
                            final Collection<VariantContext> compSet = (compSetHash == null || compSetHash.size() == 0) ? Collections.<VariantContext>emptyList() : compVCs.get(compRod).values().iterator().next();

                            // find the comp
                            final VariantContext comp = findMatchingComp(eval, compSet);

                            Collection<EvaluationContext> contextsForStratification;
                            if (perFamilyIsEnabled)
                                contextsForStratification = getEvaluationContexts(tracker, ref, eval, evalRod.getName(), comp, compRod.getName(), null, stratLevelName);
                            else {
                                String familyID;
                                if (stratLevelName.equals("all"))
                                    familyID = "all";
                                else
                                    familyID = getSampleDB().getSample(stratLevelName).getFamilyID();
                                contextsForStratification = getEvaluationContexts(tracker, ref, eval, evalRod.getName(), comp, compRod.getName(), stratLevelName, familyID);
                            }
                            for ( EvaluationContext nec : contextsForStratification ) {

                                // eval against the comp
                                synchronized (nec) {
                                    nec.apply(tracker, ref, context, comp, eval);
                                }

                                // eval=null against all comps of different type that aren't bound to another eval
                                for ( VariantContext otherComp : compSet ) {
                                    if ( otherComp != comp && ! compHasMatchingEval(otherComp, evalSetBySample) ) {
                                        synchronized (nec) {
                                            nec.apply(tracker, ref, context, otherComp, null);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if ( mergeEvals ) break; // stop processing the eval tracks
            }
        }

        return null;
    }

    /**
     * Given specific eval and comp VCs and the sample name, return an iterable
     * over all of the applicable state keys.
     *
     * this code isn't structured yet for efficiency.  Here we currently are
     * doing the following inefficient algorithm:
     *
     * for each strat:
     *   get list of relevant states that eval and comp according to strat
     *   add this list of states to a list of list states
     *
     * then
     *
     * ask the strat manager to look up all of the keys associated with the combinations
     * of these states.  For example, suppose we have a single variant S.  We have active
     * strats EvalRod, CompRod, and Novelty.  We produce a list that looks like:
     *
     *   L = [[Eval], [Comp], [All, Novel]]
     *
     * We then go through the strat manager tree to produce the keys associated with these states:
     *
     *   K = [0, 1] where EVAL x COMP x ALL = 0 and EVAL x COMP x NOVEL = 1
     *
     * It's clear that a better
     *
     * TODO -- create an inline version that doesn't create the intermediate list of list
     *
     * @param tracker
     * @param ref
     * @param eval
     * @param evalName
     * @param comp
     * @param compName
     * @param sampleName
     * @return
     */
    protected Collection<EvaluationContext> getEvaluationContexts(final RefMetaDataTracker tracker,
                                                                  final ReferenceContext ref,
                                                                  final VariantContext eval,
                                                                  final String evalName,
                                                                  final VariantContext comp,
                                                                  final String compName,
                                                                  final String sampleName,
                                                                  final String familyName) {
        final List<List<Object>> states = new LinkedList<List<Object>>();
        for ( final VariantStratifier vs : stratManager.getStratifiers() ) {
            states.add(vs.getRelevantStates(ref, tracker, comp, compName, eval, evalName, sampleName, familyName));
        }
        return stratManager.values(states);
    }


    @Requires({"comp != null", "evals != null"})
    private boolean compHasMatchingEval(final VariantContext comp, final Collection<VariantContext> evals) {
        // find all of the matching comps
        for ( final VariantContext eval : evals ) {
            if ( eval != null && doEvalAndCompMatch(comp, eval, requireStrictAlleleMatch) != EvalCompMatchType.NO_MATCH )
                return true;
        }

        // nothing matched
        return false;
    }

    private enum EvalCompMatchType { NO_MATCH, STRICT, LENIENT }

    @Requires({"eval != null", "comp != null"})
    private EvalCompMatchType doEvalAndCompMatch(final VariantContext eval, final VariantContext comp, boolean requireStrictAlleleMatch) {
        if ( comp.getType() == VariantContext.Type.NO_VARIATION || eval.getType() == VariantContext.Type.NO_VARIATION )
            // if either of these are NO_VARIATION they are LENIENT matches
            return EvalCompMatchType.LENIENT;

        if ( comp.getType() != eval.getType() )
            return EvalCompMatchType.NO_MATCH;

        // find the comp which matches both the reference allele and alternate allele from eval
        final Allele altEval = eval.getAlternateAlleles().size() == 0 ? null : eval.getAlternateAllele(0);
        final Allele altComp = comp.getAlternateAlleles().size() == 0 ? null : comp.getAlternateAllele(0);
        if ((altEval == null && altComp == null) || (altEval != null && altEval.equals(altComp) && eval.getReference().equals(comp.getReference())))
            return EvalCompMatchType.STRICT;
        else
            return requireStrictAlleleMatch ? EvalCompMatchType.NO_MATCH : EvalCompMatchType.LENIENT;
    }

    private VariantContext findMatchingComp(final VariantContext eval, final Collection<VariantContext> comps) {
        // if no comps, return null
        if ( comps == null || comps.isEmpty() )
            return null;

        // if no eval, return any comp
        if ( eval == null )
            return comps.iterator().next();

        // find all of the matching comps
        VariantContext lenientMatch = null;
        for ( final VariantContext comp : comps ) {
            switch ( doEvalAndCompMatch(comp, eval, requireStrictAlleleMatch) ) {
                case STRICT:
                    return comp;
                case LENIENT:
                    if ( lenientMatch == null ) lenientMatch = comp;
                    break;
                case NO_MATCH:
                    // do nothing
            }
        }

        // nothing matched, just return lenientMatch, which might be null
        return lenientMatch;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) { return null; }

    @Override
    public Integer reduceInit() { return null; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return null; }

    /**
     * Output the finalized report
     *
     * @param result  an integer that doesn't get used for anything
     */
    public void onTraversalDone(Integer result) {
        logger.info("Finalizing variant report");
        
        // go through the evaluations and finalize them
        for ( final EvaluationContext nec : stratManager.values() )
            for ( final VariantEvaluator ve : nec.getVariantEvaluators() )
                ve.finalizeEvaluation();

        //send data to MetricsCollection
        CompOverlap compOverlap = null;
        IndelSummary indelSummary = null;
        CountVariants countVariants = null;
        MultiallelicSummary multiallelicSummary = null;
        TiTvVariantEvaluator tiTvVariantEvaluator = null;
        MetricsCollection metricsCollection = null;
        for(final EvaluationContext nec: stratManager.values()) {
            for(final VariantEvaluator ve : nec.getVariantEvaluators()) {
                if (ve instanceof CompOverlap)
                    compOverlap = (CompOverlap) ve;
                else if (ve instanceof IndelSummary)
                    indelSummary = (IndelSummary) ve;
                else if (ve instanceof CountVariants)
                    countVariants = (CountVariants) ve;
                else if (ve instanceof MultiallelicSummary)
                    multiallelicSummary = (MultiallelicSummary) ve;
                else if (ve instanceof TiTvVariantEvaluator)
                    tiTvVariantEvaluator = (TiTvVariantEvaluator) ve;
                else if (ve instanceof MetricsCollection)
                    metricsCollection = (MetricsCollection) ve;
            }

        if(metricsCollection != null)
            metricsCollection.setData(compOverlap.concordantRate, indelSummary.n_SNPs, countVariants.nSNPs, indelSummary.n_indels, multiallelicSummary.nIndels, indelSummary.insertion_to_deletion_ratio, countVariants.insertionDeletionRatio, tiTvVariantEvaluator.tiTvRatio);
        }

        VariantEvalReportWriter.writeReport(out, stratManager, stratManager.getStratifiers(), stratManager.get(0).getVariantEvaluators());
    }

    // Accessors
    public Logger getLogger() { return logger; }

    public double getMinPhaseQuality() { return MIN_PHASE_QUALITY; }

    public int getSamplePloidy() { return ploidy; }
    public double getMendelianViolationQualThreshold() { return MENDELIAN_VIOLATION_QUAL_THRESHOLD; }

    public static String getAllSampleName() { return ALL_SAMPLE_NAME; }
    public static String getAllFamilyName() { return ALL_FAMILY_NAME; }

    public List<RodBinding<VariantContext>> getKnowns() { return knowns; }

    public List<RodBinding<VariantContext>> getEvals() { return evals; }

    public boolean isSubsettingToSpecificSamples() { return isSubsettingSamples; }
    public Set<String> getSampleNamesForEvaluation() { return sampleNamesForEvaluation; }

    public Set<String> getFamilyNamesForEvaluation() { return familyNamesForEvaluation; }

    public int getNumberOfSamplesForEvaluation() {
        if (sampleNamesForEvaluation!= null &&  !sampleNamesForEvaluation.isEmpty())
            return sampleNamesForEvaluation.size();
        else {
            return numSamplesFromArgument;
        }

    }
    public Set<String> getSampleNamesForStratification() { return sampleNamesForStratification; }

    public Set<String> getFamilyNamesForStratification() { return familyNamesForStratification; }

    public List<RodBinding<VariantContext>> getComps() { return comps; }

    public Set<SortableJexlVCMatchExp> getJexlExpressions() { return jexlExpressions; }

    public long getnProcessedLoci() {
        return nProcessedLoci;
    }

    public Set<String> getContigNames() {
        final TreeSet<String> contigs = new TreeSet<String>();
        for( final SAMSequenceRecord r :  getToolkit().getReferenceDataSource().getReference().getSequenceDictionary().getSequences()) {
            contigs.add(r.getSequenceName());
        }
        return contigs;
    }

    /**
     * getToolkit is protected, so we have to pseudo-overload it here so eval / strats can get the toolkit
     * @return
     */
    public GenomeAnalysisEngine getToolkit() {
        return super.getToolkit();
    }

    public boolean ignoreAC0Sites() {
        return ! keepSitesWithAC0;
    }
}