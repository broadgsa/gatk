package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.StandardEval;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.RequiredStratification;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.StandardStratification;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.lang.reflect.Field;
import java.util.*;

public class VariantEvalUtils {
    private final VariantEvalWalker variantEvalWalker;
    Logger logger;

    public VariantEvalUtils(VariantEvalWalker variantEvalWalker) {
        this.variantEvalWalker = variantEvalWalker;
        this.logger = variantEvalWalker.getLogger();
    }

    /**
     * List all of the available evaluation modules, then exit successfully
     */
    public void listModulesAndExit() {
        List<Class<? extends VariantStratifier>> vsClasses = new PluginManager<VariantStratifier>( VariantStratifier.class ).getPlugins();
        List<Class<? extends VariantEvaluator>> veClasses = new PluginManager<VariantEvaluator>( VariantEvaluator.class ).getPlugins();

        logger.info("Available stratification modules:");
        logger.info("(Standard modules are starred)");
        for (Class<? extends VariantStratifier> vsClass : vsClasses) {
            logger.info("\t" + vsClass.getSimpleName() + (RequiredStratification.class.isAssignableFrom(vsClass) || StandardStratification.class.isAssignableFrom(vsClass) ? "*" : ""));
        }
        logger.info("");

        logger.info("Available evaluation modules:");
        logger.info("(Standard modules are starred)");
        for (Class<? extends VariantEvaluator> veClass : veClasses) {
            logger.info("\t" + veClass.getSimpleName() + (StandardEval.class.isAssignableFrom(veClass) ? "*" : ""));
        }
        logger.info("");

        System.exit(0);
    }

    /**
     * Initialize required, standard and user-specified stratification objects
     *
     * @param noStandardStrats don't use the standard stratifications
     * @param modulesToUse     the list of stratification modules to use
     * @return set of stratifications to use
     */
    public TreeSet<VariantStratifier> initializeStratificationObjects(VariantEvalWalker variantEvalWalker, boolean noStandardStrats, String[] modulesToUse) {
        TreeSet<VariantStratifier> strats = new TreeSet<VariantStratifier>();
        Set<String> stratsToUse = new HashSet<String>();

        // Create a map for all stratification modules for easy lookup.
        HashMap<String, Class<? extends VariantStratifier>> classMap = new HashMap<String, Class<? extends VariantStratifier>>();
        for (Class<? extends VariantStratifier> c : new PluginManager<VariantStratifier>(VariantStratifier.class).getPlugins()) {
            classMap.put(c.getSimpleName(), c);
        }

        // We must use all required stratification modules.
        for (Class<? extends RequiredStratification> reqClass : new PluginManager<RequiredStratification>(RequiredStratification.class).getPlugins()) {
            if (classMap.containsKey(reqClass.getSimpleName())) {
                stratsToUse.add(reqClass.getSimpleName());
            }
        }

        // By default, use standard stratification modules.
        if (!noStandardStrats) {
            for (Class<? extends StandardStratification> stdClass : new PluginManager<StandardStratification>(StandardStratification.class).getPlugins()) {
                if (classMap.containsKey(stdClass.getSimpleName())) {
                    stratsToUse.add(stdClass.getSimpleName());
                }
            }
        }

        // Now add the user-selected modules
        stratsToUse.addAll(Arrays.asList(modulesToUse));

        // Instantiate the stratifications
        for (String module : stratsToUse) {
            if (!classMap.containsKey(module)) {
                throw new UserException.CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if (classMap.containsKey(module)) {
                Class<? extends VariantStratifier> c = classMap.get(module);

                try {
                    VariantStratifier vs = c.newInstance();
                    vs.setVariantEvalWalker(variantEvalWalker);
                    vs.initialize(variantEvalWalker.getJexlExpressions(), variantEvalWalker.getCompNames(), variantEvalWalker.getKnownNames(), variantEvalWalker.getEvalNames(), variantEvalWalker.getSampleNamesForStratification(), variantEvalWalker.getContigNames());

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

    /**
     * Initialize required, standard and user-specified evaluation objects
     *
     * @param noStandardEvals don't use the standard evaluations
     * @param modulesToUse    the list of evaluation modules to use
     * @return set of evaluations to use
     */
    public Set<Class<? extends VariantEvaluator>> initializeEvaluationObjects(boolean noStandardEvals, String[] modulesToUse) {
        Set<Class<? extends VariantEvaluator>> evals = new HashSet<Class<? extends VariantEvaluator>>();

        // Create a map for all eval modules for easy lookup.
        HashMap<String, Class<? extends VariantEvaluator>> classMap = new HashMap<String, Class<? extends VariantEvaluator>>();
        for (Class<? extends VariantEvaluator> c : new PluginManager<VariantEvaluator>(VariantEvaluator.class).getPlugins()) {
            classMap.put(c.getSimpleName(), c);
        }

        // By default, use standard eval modules.
        if (!noStandardEvals) {
            for (Class<? extends StandardEval> stdClass : new PluginManager<StandardEval>(StandardEval.class).getPlugins()) {
                if (classMap.containsKey(stdClass.getSimpleName())) {
                    evals.add(classMap.get(stdClass.getSimpleName()));
                }
            }
        }

        // Get the specific classes provided.
        for (String module : modulesToUse) {
            if (!classMap.containsKey(module)) {
                throw new UserException.CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if (classMap.containsKey(module)) {
                evals.add(classMap.get(module));
            }
        }

        return evals;
    }

    /**
     * Recursively initialize the evaluation contexts
     *
     * @param stratificationObjects the stratifications to use
     * @param evaluationObjects     the evaluations to use
     * @param stratStack            a stack of stratifications to apply
     * @param ec                    evaluation context
     * @return a map of all the evaluation contexts
     */
    public HashMap<StateKey, NewEvaluationContext> initializeEvaluationContexts(Set<VariantStratifier> stratificationObjects, Set<Class<? extends VariantEvaluator>> evaluationObjects, Stack<VariantStratifier> stratStack, NewEvaluationContext ec) {
        HashMap<StateKey, NewEvaluationContext> ecs = new HashMap<StateKey, NewEvaluationContext>();

        if (stratStack == null) {
            stratStack = new Stack<VariantStratifier>();
            stratStack.addAll(stratificationObjects);
        }

        if (!stratStack.isEmpty()) {
            Stack<VariantStratifier> newStratStack = new Stack<VariantStratifier>();
            newStratStack.addAll(stratStack);

            VariantStratifier vs = newStratStack.pop();

            for (String state : vs.getAllStates()) {
                NewEvaluationContext nec = new NewEvaluationContext();
                if (ec != null) {
                    nec.putAll(ec);
                }
                nec.put(vs, state);

                ecs.putAll(initializeEvaluationContexts(stratificationObjects, evaluationObjects, newStratStack, nec));
            }
        } else {
            HashMap<StateKey, NewEvaluationContext> necs = new HashMap<StateKey, NewEvaluationContext>();

            StateKey stateKey = new StateKey();
            for (VariantStratifier vs : ec.keySet()) {
                String state = ec.get(vs);

                stateKey.put(vs.getClass().getSimpleName(), state);
            }

            ec.addEvaluationClassList(variantEvalWalker, stateKey, evaluationObjects);

            necs.put(stateKey, ec);

            return necs;
        }

        return ecs;
    }

    /**
     * Initialize the output report
     *
     * @param stratificationObjects the stratifications to use
     * @param evaluationObjects     the evaluations to use
     * @return an initialized report object
     */
    public GATKReport initializeGATKReport(Set<VariantStratifier> stratificationObjects, Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        GATKReport report = new GATKReport();

        for (Class<? extends VariantEvaluator> ve : evaluationObjects) {
            String tableName = ve.getSimpleName();
            String tableDesc = ve.getAnnotation(Analysis.class).description();

            report.addTable(tableName, tableDesc);

            GATKReportTable table = report.getTable(tableName);
            table.addPrimaryKey("entry", false);
            table.addColumn(tableName, tableName);

            for (VariantStratifier vs : stratificationObjects) {
                String columnName = vs.getClass().getSimpleName();

                table.addColumn(columnName, "unknown");
            }

            try {
                VariantEvaluator vei = ve.newInstance();
                vei.initialize(variantEvalWalker);

                AnalysisModuleScanner scanner = new AnalysisModuleScanner(vei);
                Map<Field, DataPoint> datamap = scanner.getData();

                for (Field field : datamap.keySet()) {
                    field.setAccessible(true);

                    if (!(field.get(vei) instanceof TableType)) {
                        table.addColumn(field.getName(), 0.0);
                    }
                }
            } catch (InstantiationException e) {
                throw new StingException("InstantiationException: " + e);
            } catch (IllegalAccessException e) {
                throw new StingException("IllegalAccessException: " + e);
            }
        }

        return report;
    }

    /**
     * Figure out what the allowable variation types are based on the eval context
     *
     * @param tracker   the reference metadata tracker
     * @param ref       the reference context
     * @param compNames the comp track names
     * @param evalNames the evaluation track names
     * @return the set of allowable variation types
     */
    public EnumSet<VariantContext.Type> getAllowableVariationTypes(RefMetaDataTracker tracker,
                                                                   ReferenceContext ref,
                                                                   Set<String> compNames,
                                                                   Set<String> evalNames,
                                                                   boolean dynamicSelectTypes ) {
        if ( dynamicSelectTypes ) { // todo -- this code is really conceptually broken
            EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.NO_VARIATION);

            if (tracker != null) {
                Collection<VariantContext> evalvcs = tracker.getVariantContexts(ref, evalNames, null, ref.getLocus(), true, false);

                for (VariantContext vc : evalvcs) {
                    allowableTypes.add(vc.getType());
                }

                if (allowableTypes.size() == 1) {
                    // We didn't find any variation in the eval track, so now let's look at the comp track for allowable types
                    Collection<VariantContext> compvcs = tracker.getVariantContexts(ref, compNames, null, ref.getLocus(), true, false);

                    for (VariantContext vc : compvcs) {
                        allowableTypes.add(vc.getType());
                    }
                }
            }

            return allowableTypes;
        } else {
            return EnumSet.allOf(VariantContext.Type.class);
        }
    }

    /**
     * Subset a VariantContext to a single sample
     *
     * @param vc         the VariantContext object containing multiple samples
     * @param sampleName the sample to pull out of the VariantContext
     * @return a new VariantContext with just the requested sample
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, String sampleName) {
        ArrayList<String> sampleNames = new ArrayList<String>();
        sampleNames.add(sampleName);

        return getSubsetOfVariantContext(vc, sampleNames);
    }

    /**
     * Subset a VariantContext to a set of samples
     *
     * @param vc          the VariantContext object containing multiple samples
     * @param sampleNames the samples to pull out of the VariantContext
     * @return a new VariantContext with just the requested samples
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, Collection<String> sampleNames) {
        VariantContext vcsub = vc.subContextFromGenotypes(vc.getGenotypes(sampleNames).values());

        HashMap<String, Object> newAts = new HashMap<String, Object>(vcsub.getAttributes());

        int originalAlleleCount = vc.getHetCount() + 2 * vc.getHomVarCount();
        int newAlleleCount = vcsub.getHetCount() + 2 * vcsub.getHomVarCount();

        if (originalAlleleCount == newAlleleCount && newAlleleCount == 1) {
            newAts.put("ISSINGLETON", true);
        }

        VariantContextUtils.calculateChromosomeCounts(vcsub, newAts, true);
        vcsub = VariantContext.modifyAttributes(vcsub, newAts);

        //VariantEvalWalker.logger.debug(String.format("VC %s subset to %s AC%n", vc.getSource(), vc.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY)));

        return vcsub;
    }

    /**
     * For a list of track names, bind the variant contexts to a trackName->sampleName->VariantContext mapping.
     * Additional variant contexts per sample are automatically generated and added to the map unless the sample name
     * matches the ALL_SAMPLE_NAME constant.
     *
     * @param tracker        the metadata tracker
     * @param ref            the reference context
     * @param trackNames     the list of track names to process
     * @param allowableTypes a set of allowable variation types
     * @param byFilter       if false, only accept PASSing VariantContexts.  Otherwise, accept both PASSing and filtered
     *                       sites
     * @param subsetBySample if false, do not separate the track into per-sample VCs
     * @param trackPerSample if false, don't stratify per sample (and don't cut up the VariantContext like we would need
     *                       to do this)
     * @return a mapping of track names to a list of VariantContext objects
     */
    public HashMap<String, HashMap<String, VariantContext>> bindVariantContexts(RefMetaDataTracker tracker, ReferenceContext ref, Set<String> trackNames, EnumSet<VariantContext.Type> allowableTypes, boolean byFilter, boolean subsetBySample, boolean trackPerSample) {
        HashMap<String, HashMap<String, VariantContext>> bindings = new HashMap<String, HashMap<String, VariantContext>>();

        for (String trackName : trackNames) {
            HashMap<String, VariantContext> vcs = new HashMap<String, VariantContext>();

            Collection<VariantContext> contexts = tracker == null ? null : tracker.getVariantContexts(ref, trackName, allowableTypes, ref.getLocus(), true, true);
            VariantContext vc = contexts != null && contexts.size() == 1 ? contexts.iterator().next() : null;

            // First, filter the VariantContext to represent only the samples for evaluation
            if (vc != null) {
                VariantContext vcsub = vc;

                if (subsetBySample && vc.hasGenotypes() && vc.hasGenotypes(variantEvalWalker.getSampleNamesForEvaluation())) {
                    vcsub = getSubsetOfVariantContext(vc, variantEvalWalker.getSampleNamesForEvaluation());
                }

                if ((byFilter || !vcsub.isFiltered())) {
                    vcs.put(VariantEvalWalker.getAllSampleName(), vcsub);
                }

                // Now, if stratifying, split the subsetted vc per sample and add each as a new context
                if (vc.hasGenotypes() && trackPerSample) {
                    for (String sampleName : variantEvalWalker.getSampleNamesForEvaluation()) {
                        VariantContext samplevc = getSubsetOfVariantContext(vc, sampleName);

                        if ((byFilter || !samplevc.isFiltered())) {
                            vcs.put(sampleName, samplevc);
                        }
                    }
                }

                bindings.put(trackName, vcs);
            }
        }

        return bindings;
    }

    /**
     * Maps track names to sample name to VariantContext objects.  For eval tracks, VariantContexts per specified sample
     * are also included.
     *
     * @param tracker   the metadata tracker
     * @param ref       the reference context
     * @param compNames the list of comp names to process
     * @param evalNames the list of eval names to process
     * @return a mapping of track names to a list of VariantContext objects
     */
    public HashMap<String, HashMap<String, VariantContext>> getVariantContexts(RefMetaDataTracker tracker, ReferenceContext ref, Set<String> compNames, Set<String> evalNames, boolean dynamicSelectTypes) {
        HashMap<String, HashMap<String, VariantContext>> vcs = new HashMap<String, HashMap<String, VariantContext>>();

        EnumSet<VariantContext.Type> allowableTypes = getAllowableVariationTypes(tracker, ref, compNames, evalNames, dynamicSelectTypes);

        boolean byFilter = false;
        boolean perSampleIsEnabled = false;
        for (VariantStratifier vs : variantEvalWalker.getStratificationObjects()) {
            if (vs.getClass().getSimpleName().equals("Filter")) {
                byFilter = true;
            } else if (vs.getClass().getSimpleName().equals("Sample")) {
                perSampleIsEnabled = true;
            }
        }

        HashMap<String, HashMap<String, VariantContext>> evalBindings = bindVariantContexts(tracker, ref, evalNames, allowableTypes, byFilter, true, perSampleIsEnabled);
        HashMap<String, HashMap<String, VariantContext>> compBindings = bindVariantContexts(tracker, ref, compNames, allowableTypes, byFilter, false, false);

        vcs.putAll(compBindings);
        vcs.putAll(evalBindings);

        return vcs;
    }

    /**
     * Recursively initialize the state keys used to look up the right evaluation context based on the state of the
     * variant context
     *
     * @param stateMap   the map of allowable states
     * @param stateStack a stack of the states
     * @param stateKey   a state key object
     * @param stateKeys  all the state keys
     * @return a list of state keys
     */
    public ArrayList<StateKey> initializeStateKeys(HashMap<VariantStratifier, ArrayList<String>> stateMap, Stack<HashMap<VariantStratifier, ArrayList<String>>> stateStack, StateKey stateKey, ArrayList<StateKey> stateKeys) {
        if (stateStack == null) {
            stateStack = new Stack<HashMap<VariantStratifier, ArrayList<String>>>();

            for (VariantStratifier vs : stateMap.keySet()) {
                HashMap<VariantStratifier, ArrayList<String>> oneSetOfStates = new HashMap<VariantStratifier, ArrayList<String>>();
                oneSetOfStates.put(vs, stateMap.get(vs));

                stateStack.add(oneSetOfStates);
            }
        }

        if (!stateStack.isEmpty()) {
            Stack<HashMap<VariantStratifier, ArrayList<String>>> newStateStack = new Stack<HashMap<VariantStratifier, ArrayList<String>>>();
            newStateStack.addAll(stateStack);

            HashMap<VariantStratifier, ArrayList<String>> oneSetOfStates = newStateStack.pop();
            VariantStratifier vs = oneSetOfStates.keySet().iterator().next();

            for (String state : oneSetOfStates.get(vs)) {
                StateKey newStateKey = new StateKey();
                if (stateKey != null) {
                    newStateKey.putAll(stateKey);
                }

                newStateKey.put(vs.getClass().getSimpleName(), state);

                initializeStateKeys(stateMap, newStateStack, newStateKey, stateKeys);
            }
        } else {
            stateKeys.add(stateKey);

            return stateKeys;
        }

        return stateKeys;
    }
}