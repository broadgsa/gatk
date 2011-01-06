package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.evaluators.StandardEval;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications.RequiredStratification;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications.StandardStratification;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.tags.Analysis;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.tags.DataPoint;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.util.AnalysisModuleScanner;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.util.NewEvaluationContext;
import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.util.StateKey;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.*;

/**
 * General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ts/Tv ratios, and a lot more)
 */
@Reference(window=@Window(start=-50, stop=50))
public class NewVariantEvalWalker extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    // Output arguments
    @Output
    protected PrintStream out;

    // Help arguments
    @Argument(fullName="list", shortName="ls", doc="List the available eval modules and exit")
    protected Boolean LIST = false;

    // Partitioning the data arguments
    @Argument(shortName="select", doc="One or more stratifications to use when evaluating the data", required=false)
    protected ArrayList<String> SELECT_EXPS = new ArrayList<String>();

    @Argument(shortName="selectName", doc="Names to use for the list of stratifications (must be a 1-to-1 mapping)", required=false)
    protected ArrayList<String> SELECT_NAMES = new ArrayList<String>();

    @Argument(shortName="sample", doc="Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context", required=false)
    protected Set<String> SAMPLE_EXPRESSIONS;

    @Argument(shortName="knownName", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected String[] KNOWN_NAMES = {DbSNPHelper.STANDARD_DBSNP_TRACK_NAME};

    // Stratification arguments
    @Argument(fullName="stratificationModule", shortName="ST", doc="One or more specific stratification modules to apply to the eval track(s) (in addition to the standard stratifications, unless -noS is specified)", required=false)
    protected String[] STRATIFICATIONS_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardStratifications", shortName="noST", doc="Do not use the standard stratification modules by default (instead, only those that are specified with the -S option)")
    protected Boolean NO_STANDARD_STRATIFICATIONS = false;

    // Evaluator arguments
    @Argument(fullName="evalModule", shortName="EV", doc="One or more specific eval modules to apply to the eval track(s) (in addition to the standard modules, unless -noE is specified)", required=false)
    protected String[] MODULES_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardModules", shortName="noEV", doc="Do not use the standard modules by default (instead, only those that are specified with the -E option)")
    protected Boolean NO_STANDARD_MODULES = false;

    // Variables
    private Set<VariantContextUtils.JexlVCMatchExp> jexlExpressions = new TreeSet<VariantContextUtils.JexlVCMatchExp>();
    private Set<String> compNames = new TreeSet<String>();
    private Set<String> knownNames = new TreeSet<String>();
    private Set<String> evalNames = new TreeSet<String>();
    private Set<String> sampleNames = new TreeSet<String>();

    // The list of stratifiers and evaluators to use
    private TreeSet<VariantStratifier> stratificationObjects = null;
    private Set<Class<? extends VariantEvaluator>> evaluationObjects = null;

    // The set of all possible evaluation contexts
    private HashMap<StateKey, NewEvaluationContext> evaluationContexts = null;

    // Output report
    private GATKReport report = null;

    // Public constants
    public static String ALL_SAMPLE_NAME = "all";

    /**
     * List all of the available evaluation modules, then exit successfully
     */
    private void listModulesAndExit() {
        List<Class<? extends VariantStratifier>> vsClasses = new PluginManager<VariantStratifier>( VariantStratifier.class ).getPlugins();
        List<Class<? extends VariantEvaluator>> veClasses = new PluginManager<VariantEvaluator>( VariantEvaluator.class ).getPlugins();

        logger.info("Available stratifcation modules:");
        logger.info("(Standard modules are starred)");
        for (Class<? extends VariantStratifier> vsClass : vsClasses) {
            logger.info("\t" + vsClass.getSimpleName() + (RequiredStratification.class.isAssignableFrom(vsClass) || StandardStratification.class.isAssignableFrom(vsClass) ? "*" : ""));
        }
        logger.info("");

        logger.info("Available eval modules:");
        logger.info("(Standard modules are starred)");
        for (Class<? extends VariantEvaluator> veClass : veClasses) {
            logger.info("\t" + veClass.getSimpleName() + (StandardEval.class.isAssignableFrom(veClass) ? "*" : ""));
        }
        logger.info("");

        System.exit(0);
    }

    private TreeSet<VariantStratifier> initializeStratificationObjects(boolean noStandardStrats, String[] modulesToUse) {
        TreeSet<VariantStratifier> strats = new TreeSet<VariantStratifier>();
        Set<String> stratsToUse = new HashSet<String>();

        // Create a map for all stratification modules for easy lookup.
        HashMap<String, Class<? extends VariantStratifier>> classMap = new HashMap<String, Class<? extends VariantStratifier>>();
        for ( Class<? extends VariantStratifier> c : new PluginManager<VariantStratifier>( VariantStratifier.class ).getPlugins() ) {
            classMap.put(c.getSimpleName(), c);
        }

        // We must use all required stratification modules.
        for ( Class<? extends RequiredStratification> reqClass : new PluginManager<RequiredStratification>( RequiredStratification.class ).getPlugins() ) {
            if ( classMap.containsKey(reqClass.getSimpleName()) ) {
                stratsToUse.add(reqClass.getSimpleName());
            }
        }

        // By default, use standard stratification modules.
        if ( !noStandardStrats ) {
            for ( Class<? extends StandardStratification> stdClass : new PluginManager<StandardStratification>( StandardStratification.class ).getPlugins() ) {
                if ( classMap.containsKey(stdClass.getSimpleName()) ) {
                    stratsToUse.add(stdClass.getSimpleName());
                }
            }
        }

        // Now add the user-selected modules
        stratsToUse.addAll(Arrays.asList(modulesToUse));

        // Instantiate the stratifications
        for ( String module : stratsToUse ) {
            if ( !classMap.containsKey(module) ) {
                throw new UserException.CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if ( classMap.containsKey(module) ) {
                Class<? extends VariantStratifier> c = classMap.get(module);

                try {
                    VariantStratifier vs = c.newInstance();
                    vs.initialize(jexlExpressions, compNames, knownNames, evalNames, sampleNames);

                    strats.add(vs);
                } catch (InstantiationException e) {
                    throw new StingException("Unable to instantiate stratification module '" + c.getSimpleName() + "'");
                } catch (IllegalAccessException e) {
                    throw new StingException("Illegal access error when trying to instantiate stratification module '" + c.getSimpleName() + "'");
                }
            }
        }

        return strats;
    }

    private Set<Class<? extends VariantEvaluator>> initializeEvaluationObjects(boolean noStandardEvals, String[] modulesToUse) {
        Set<Class<? extends VariantEvaluator>> evals = new HashSet<Class<? extends VariantEvaluator>>();

        // Create a map for all eval modules for easy lookup.
        HashMap<String, Class<? extends VariantEvaluator>> classMap = new HashMap<String, Class<? extends VariantEvaluator>>();
        for ( Class<? extends VariantEvaluator> c : new PluginManager<VariantEvaluator>( VariantEvaluator.class ).getPlugins() ) {
            classMap.put(c.getSimpleName(), c);
        }

        // By default, use standard eval modules.
        if ( !noStandardEvals ) {
            for ( Class<? extends StandardEval> stdClass : new PluginManager<StandardEval>( StandardEval.class ).getPlugins() ) {
                if ( classMap.containsKey(stdClass.getSimpleName()) ) {
                    evals.add(classMap.get(stdClass.getSimpleName()));
                }
            }
        }

        // Get the specific classes provided.
        for ( String module : modulesToUse ) {
            if ( !classMap.containsKey(module) ) {
                throw new UserException.CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if ( classMap.containsKey(module) ) {
                evals.add(classMap.get(module));
            }
        }

        return evals;
    }

    private HashMap<StateKey, NewEvaluationContext> initializeEvaluationContexts(Set<VariantStratifier> stratificationObjects,  Set<Class<? extends VariantEvaluator>> evaluationObjects, Stack<VariantStratifier> stratStack, NewEvaluationContext ec) {
        HashMap<StateKey, NewEvaluationContext> ecs = new HashMap<StateKey, NewEvaluationContext>();

        if (stratStack == null) {
            stratStack = new Stack<VariantStratifier>();
            stratStack.addAll(stratificationObjects);
        }

        if (!stratStack.isEmpty()) {
            Stack<VariantStratifier> newStratStack = new Stack<VariantStratifier>();
            newStratStack.addAll(stratStack);

            NewEvaluationContext nec = new NewEvaluationContext();
            if (ec != null) {
                nec.putAll(ec);
            }

            VariantStratifier vs = newStratStack.pop();

            for ( String state : vs.getAllStates() ) {
                nec.put(vs, state);

                ecs.putAll(initializeEvaluationContexts(stratificationObjects, evaluationObjects, newStratStack, nec));
            }
        } else {
            HashMap<StateKey, NewEvaluationContext> necs = new HashMap<StateKey, NewEvaluationContext>();

            StateKey statekey = new StateKey();
            for ( VariantStratifier vs : ec.keySet() ) {
                String state = ec.get(vs);

                statekey.put(vs.getClass().getSimpleName(), state);
            }

            ec.addEvaluationClassList(evaluationObjects);

            necs.put(statekey, ec);

            return necs;
        }

        return ecs;
    }

    private GATKReport initializeGATKReport(Set<VariantStratifier> stratificationObjects, Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        GATKReport report = new GATKReport();

        for ( Class<? extends VariantEvaluator> ve : evaluationObjects ) {
            String tableName = ve.getSimpleName();
            String tableDesc = ve.getAnnotation(Analysis.class).description();

            report.addTable(tableName, tableDesc);

            GATKReportTable table = report.getTable(tableName);
            table.addPrimaryKey("entry", false);

            for ( VariantStratifier vs : stratificationObjects ) {
                String columnName = vs.getClass().getSimpleName();

                columnName = columnName.replace("Stratifier", "");
                columnName = columnName.replace("Status", "");

                table.addColumn(columnName, "unknown");
            }

            AnalysisModuleScanner scanner = new AnalysisModuleScanner(ve);
            Map<Field, DataPoint> datamap = scanner.getData();

            for (Field field : datamap.keySet()) {
                table.addColumn(field.getName(), 0.0);
            }
        }

        return report;
    }

    public void initialize() {
        // Just list the modules, and exit quickly.
        if (LIST) { listModulesAndExit(); }

        // Categorize each rod as an eval or a comp rod.
        for ( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if ( d.getName().startsWith("eval") ) {
                evalNames.add(d.getName());
            } else if ( d.getName().startsWith("comp") || d.getName().startsWith(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                compNames.add(d.getName());
            } else {
                logger.info(String.format("Not evaluating ROD binding '%s' because the name did not start with %s, comp, or eval", d.getName(), Utils.join(", ", KNOWN_NAMES)));
            }
        }

        // Set up set of known names
        knownNames.addAll(Arrays.asList(KNOWN_NAMES));

        // Now that we have all the rods categorized, determine the sample list from the eval rods.
        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), evalNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        sampleNames.addAll(SampleUtils.getSamplesFromCommandLineInput(vcfSamples, SAMPLE_EXPRESSIONS));

        // Initialize select expressions
        jexlExpressions.addAll(VariantContextUtils.initializeMatchExps(SELECT_NAMES, SELECT_EXPS));

        // Initialize the set of stratifications and evaluations to use
        stratificationObjects = initializeStratificationObjects(NO_STANDARD_STRATIFICATIONS, STRATIFICATIONS_TO_USE);
        evaluationObjects = initializeEvaluationObjects(NO_STANDARD_MODULES, MODULES_TO_USE);

        // Initialize the evaluation contexts
        evaluationContexts = initializeEvaluationContexts(stratificationObjects, evaluationObjects, null, null);

        // Initialize report table
        report = initializeGATKReport(stratificationObjects, evaluationObjects);
    }

    private EnumSet<VariantContext.Type> getAllowableVariationTypes(RefMetaDataTracker tracker, ReferenceContext ref, Set<String> evalNames) {
        EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.NO_VARIATION);

        Collection<VariantContext> vcs = tracker.getVariantContexts(ref, evalNames, null, ref.getLocus(), true, false);

        for ( VariantContext vc : vcs ) {
            allowableTypes.add(vc.getType());
        }

        return allowableTypes;
    }

    /**
     * For a list of track names, bind the variant contexts to a trackName->sampleName->VariantContext mapping.
     * Additional variant contexts per sample are automatically generated and added to the map unless the
     * sample name matches the ALL_SAMPLE_NAME constant.
     *
     * @param tracker         the metadata tracker
     * @param ref             the reference context
     * @param trackNames      the list of track names to process
     * @param sampleNames     the list of samples to include
     * @param allowableTypes  a set of allowable variation types
     * @return  a mapping of track names to a list of VariantContext objects
     */
    private HashMap<String, HashMap<String, VariantContext>> bindVariantContexts(RefMetaDataTracker tracker, ReferenceContext ref, Set<String> trackNames, Set<String> sampleNames, EnumSet<VariantContext.Type> allowableTypes) {
        HashMap<String, HashMap<String, VariantContext>> bindings = new HashMap<String, HashMap<String, VariantContext>>();

        for ( String trackName : trackNames ) {
            Collection<VariantContext> contexts = tracker.getVariantContexts(ref, trackName, allowableTypes, ref.getLocus(), true, true);

            VariantContext vc = contexts.size() == 1 ? contexts.iterator().next() : null;

            HashMap<String, VariantContext> vcs = new HashMap<String, VariantContext>();

            if ( vc != null ) {
                for ( String sample : sampleNames ) {
                    VariantContext vcsub = vc;

                    if (!sample.equals(ALL_SAMPLE_NAME)) {
                        vcsub = vcsub.subContextFromGenotypes(vcsub.getGenotype(sample));

                        HashMap<String,Object> newAts = new HashMap<String,Object>(vcsub.getAttributes());
                        VariantContextUtils.calculateChromosomeCounts(vcsub,newAts,true);
                        vcsub = VariantContext.modifyAttributes(vcsub,newAts);

                        logger.debug(String.format("VC %s subset to %s AC%n", vc.getSource(), vc.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY)));
                    }

                    vcs.put(sample, vcsub);
                }

                bindings.put(trackName, vcs);
            }
        }

        return bindings;
    }

    /**
     * Maps track names to sample name to VariantContext objects.  For eval tracks, VariantContexts per specified
     * sample are also included.
     *
     * @param tracker      the metadata tracker
     * @param ref          the reference context
     * @param compNames    the list of comp names to process
     * @param evalNames    the list of eval names to process
     * @param sampleNames  the list of samples to include
     * @return  a mapping of track names to a list of VariantContext objects
     */
    private HashMap<String, HashMap<String, VariantContext>> getVariantContexts(RefMetaDataTracker tracker, ReferenceContext ref, Set<String> compNames, Set<String> evalNames, Set<String> sampleNames) {
        HashMap<String, HashMap<String, VariantContext>> vcs = new HashMap<String, HashMap<String, VariantContext>>();

        Set<String> allSamplesList = new HashSet<String>();
        allSamplesList.add(ALL_SAMPLE_NAME);

        EnumSet<VariantContext.Type> allowableTypes = getAllowableVariationTypes(tracker, ref, evalNames);

        HashMap<String, HashMap<String, VariantContext>> compBindings = bindVariantContexts(tracker, ref, compNames, allSamplesList, allowableTypes);
        HashMap<String, HashMap<String, VariantContext>> evalBindings = bindVariantContexts(tracker, ref, evalNames, sampleNames, allowableTypes);

        //HashMap<String, HashMap<String, VariantContext>> evalBindings;
        //if (stratificationObjects.contains(SampleStratifier) {
        //    evalBindings = bindVariantContexts(tracker, ref, evalNames, sampleNames, allowableTypes);
        //} else {
        //    evalBindings = bindVariantContexts(tracker, ref, evalNames, allSamplesList, allowableTypes);
        //}

        vcs.putAll(compBindings);
        vcs.putAll(evalBindings);

        return vcs;
    }

    private ArrayList<StateKey> initializeStateKeys(HashMap<VariantStratifier, ArrayList<String>> stateMap, Stack<HashMap<VariantStratifier, ArrayList<String>>> stateStack, StateKey stateKey) {
        ArrayList<StateKey> stateKeys = new ArrayList<StateKey>();

        if (stateStack == null) {
            stateStack = new Stack<HashMap<VariantStratifier, ArrayList<String>>>();

            for ( VariantStratifier vs : stateMap.keySet() ) {
                HashMap<VariantStratifier, ArrayList<String>> oneSetOfStates = new HashMap<VariantStratifier, ArrayList<String>>();
                oneSetOfStates.put(vs, stateMap.get(vs));

                stateStack.add(oneSetOfStates);
            }
        }

        if (!stateStack.isEmpty()) {
            Stack<HashMap<VariantStratifier, ArrayList<String>>> newStateStack = new Stack<HashMap<VariantStratifier, ArrayList<String>>>();
            newStateStack.addAll(stateStack);

            StateKey newStateKey = new StateKey();
            if (stateKey != null) {
                newStateKey.putAll(stateKey);
            }

            HashMap<VariantStratifier, ArrayList<String>> oneSetOfStates = newStateStack.pop();
            VariantStratifier vs = oneSetOfStates.keySet().iterator().next();

            for ( String state : oneSetOfStates.get(vs)) {
                newStateKey.put(vs.getClass().getSimpleName(), state);

                stateKeys.addAll(initializeStateKeys(stateMap, newStateStack, newStateKey));
            }
        } else {
            ArrayList<StateKey> newStateKeys = new ArrayList<StateKey>();

            newStateKeys.add(stateKey);

            return newStateKeys;
        }

        return stateKeys;
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            //      track           sample  vc
            HashMap<String, HashMap<String, VariantContext>> vcs = getVariantContexts(tracker, ref, compNames, evalNames, sampleNames);

            for ( String compName : compNames ) {
                VariantContext comp = vcs.containsKey(compName) && vcs.get(compName) != null && vcs.get(compName).containsKey(ALL_SAMPLE_NAME) ? vcs.get(compName).get(ALL_SAMPLE_NAME) : null;

                for ( String evalName : evalNames ) {
                    for ( String sampleName : sampleNames ) {
                        VariantContext eval = vcs.get(evalName) == null ? null : vcs.get(evalName).get(sampleName);

                        HashMap<VariantStratifier, ArrayList<String>> stateMap = new HashMap<VariantStratifier, ArrayList<String>>();
                        for ( VariantStratifier vs : stratificationObjects ) {
                            ArrayList<String> states = vs.getRelevantStates(ref, comp, eval, sampleName);
                            stateMap.put(vs, states);
                        }

                        ArrayList<StateKey> stateKeys = initializeStateKeys(stateMap, null, null);

                        for ( StateKey stateKey : stateKeys ) {
                            NewEvaluationContext nec = evaluationContexts.get(stateKey);

                            nec.apply(tracker, ref, context, comp, eval);
                        }

                        //logger.info(ref.getLocus());
                        //logger.info("\tcomp: " + comp);
                        //logger.info("\teval: " + eval);
                    }
                }
            }
        }

        return null;
    }

    /**
     * A composite, 'reduce of reduces' function.
     *
     * @param lhs 'left-most' portion of data in the composite reduce.
     * @param rhs 'right-most' portion of data in the composite reduce.
     * @return The composite reduce type.
     */
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer result) {
        for ( StateKey stateKey : evaluationContexts.keySet() ) {
            NewEvaluationContext nec = evaluationContexts.get(stateKey);

            for ( VariantEvaluator ve : nec.getEvaluationClassList().values() ) {
                ve.finalizeEvaluation();

                AnalysisModuleScanner scanner = new AnalysisModuleScanner(ve);
                Map<Field, DataPoint> datamap = scanner.getData();

                GATKReportTable table = report.getTable(ve.getClass().getSimpleName());

                for ( VariantStratifier vs : stratificationObjects ) {
                    String columnName = vs.getClass().getSimpleName();

                    columnName = columnName.replace("Stratifier", "");
                    columnName = columnName.replace("Status", "");

                    table.set(stateKey.toString(), columnName, stateKey.get(vs.getClass().getSimpleName()));
                }

                for (Field field : datamap.keySet()) {
                    try {
                        field.setAccessible(true);
                        table.set(stateKey.toString(), field.getName(), field.get(ve));
                    } catch (IllegalAccessException e) {
                        throw new ReviewedStingException("Unable to access a data field in a VariantEval analysis module: " + e);
                    }
                }
            }
        }

        report.print(out);
    }
}