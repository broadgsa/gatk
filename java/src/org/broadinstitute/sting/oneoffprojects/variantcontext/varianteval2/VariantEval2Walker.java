package org.broadinstitute.sting.oneoffprojects.variantcontext.varianteval2;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.oneoffprojects.variantcontext.VariantContext;
import org.broadinstitute.sting.oneoffprojects.variantcontext.VariantContextAdaptors;
import org.broadinstitute.sting.oneoffprojects.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Test routine for new VariantContext object
 */
public class VariantEval2Walker extends RodWalker<Integer, Integer> {
    // todo -- add doc string
    @Argument(shortName="select", doc="", required=false)
    protected String[] SELECT_EXPS = {"QUAL > 500.0", "HARD_TO_VALIDATE==1", "GATK_STANDARD==1"};

    // todo -- add doc string
    @Argument(shortName="selectName", doc="", required=false)
    protected String[] SELECT_NAMES = {"q500plus", "low_mapq", "gatk_std_filters"};

    @Argument(shortName="known", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected String[] KNOWN_NAMES = {"dbsnp"};

    /** private class holding all of the information about a single evaluation group (e.g., for eval ROD) */
    private class EvaluationContext extends HashMap<String, Set<VariantEvaluator>> {
        // useful for typing
        public String trackName, contextName;
        VariantContextUtils.MatchExp selectExp;

        public EvaluationContext(String trackName, String contextName, VariantContextUtils.MatchExp selectExp) {
            this.trackName = trackName;
            this.contextName = contextName;
            this.selectExp = selectExp;
        }
    }

    private HashMap<String, EvaluationContext> contexts = new HashMap<String, EvaluationContext>();
    private Set<String> compNames = new HashSet<String>();

    private static String RAW_SET_NAME      = "raw";
    private static String RETAINED_SET_NAME = "called";
    private static String FILTERED_SET_NAME = "filtered";
    private static String ALL_SET_NAME      = "all";
    private static String KNOWN_SET_NAME    = "known";
    private static String NOVEL_SET_NAME    = "novel";

    // Dynamically determined variantEvaluation classes
    private List<String> variantEvaluationNames = new ArrayList<String>();
    private List<Class<? extends VariantEvaluator>> evaluationClasses = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    // --------------------------------------------------------------------------------------------------------------
    public void initialize() {
        determineAllEvalations();
        List<VariantContextUtils.MatchExp> selectExps = VariantContextUtils.initializeMatchExps(SELECT_NAMES, SELECT_EXPS);

        for ( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if ( d.getName().startsWith("eval") ) {
                for ( VariantContextUtils.MatchExp e : selectExps ) {
                    addNewContext(d.getName(), d.getName() + "." + e.name, e);
                }
                addNewContext(d.getName(), d.getName(), null);
            } else if ( d.getName().startsWith("dbsnp") || d.getName().startsWith("hapmap") || d.getName().startsWith("comp") ) {
                compNames.add(d.getName());
            } else {
                logger.info("Not evaluating ROD binding " + d.getName());
            }
        }
    }

    private void determineAllEvalations() {
        evaluationClasses = PackageUtils.getClassesImplementingInterface(VariantEvaluator.class);
        for ( VariantEvaluator e : instantiateEvalationsSet() ) {
            // for collecting purposes
            variantEvaluationNames.add(e.getName());
            logger.debug("Including VariantEvaluator " + e.getName() + " of class " + e.getClass());
        }

        Collections.sort(variantEvaluationNames);
    }

    private Set<VariantEvaluator> instantiateEvalationsSet() {
        Set<VariantEvaluator> evals = new HashSet<VariantEvaluator>();
        for ( Class c : evaluationClasses ) {
            try {
                evals.add((VariantEvaluator) c.newInstance());
            } catch (InstantiationException e) {
                throw new StingException(String.format("Cannot instantiate annotation class '%s': must be concrete class", c.getSimpleName()));
            } catch (IllegalAccessException e) {
                throw new StingException(String.format("Cannot instantiate annotation class '%s': must have no-arg constructor", c.getSimpleName()));
            }
        }

        return evals;
    }

    private void addNewContext(String trackName, String contextName, VariantContextUtils.MatchExp selectExp) {
        EvaluationContext group = new EvaluationContext(trackName, contextName, selectExp);

        for ( String filteredName : Arrays.asList(RAW_SET_NAME, RETAINED_SET_NAME, FILTERED_SET_NAME) ) {
            for ( String subname : Arrays.asList(ALL_SET_NAME, KNOWN_SET_NAME, NOVEL_SET_NAME) ) {
                group.put(subname + "." + filteredName, instantiateEvalationsSet());
            }
        }

        contexts.put(contextName, group);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    // --------------------------------------------------------------------------------------------------------------

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //System.out.printf("map at %s with %d skipped%n", context.getLocation(), context.getSkippedBases());

        Map<String, VariantContext> comps = getCompVariantContexts(tracker, context);

        // to enable walking over pairs where eval or comps have no elements
        for ( EvaluationContext group : contexts.values() ) {
            VariantContext vc = getVariantContext(group.trackName, tracker, context);

            //logger.debug(String.format("Updating %s of %s with variant", group.name, vc));
            for ( Map.Entry<String, Set<VariantEvaluator>> namedEvaluations : group.entrySet() ) {
                String evaluationName = namedEvaluations.getKey();
                Set<VariantEvaluator> evaluations = namedEvaluations.getValue();
                boolean evalWantsVC = applyVCtoEvaluation(evaluationName, vc, comps, group);

                for ( VariantEvaluator evaluation : evaluations ) {
                    // we always call update0 in case the evaluation tracks things like number of bases covered
                    evaluation.update0(tracker, ref, context);

                    // now call the single or paired update function
                    switch ( evaluation.getComparisonOrder() ) {
                        case 1:
                            if ( vc != null ) evaluation.update1(vc, tracker, ref, context);
                            break;
                        case 2:
                            for ( VariantContext comp : comps.values() ) {
                                evaluation.update2( evalWantsVC ? vc : null, comp, tracker, ref, context);
                            }
                            break;
                        default:
                            throw new StingException("BUG: Unexpected evaluation order " + evaluation);
                    }
                }
            }
        }

        return 0;
    }

    private boolean applyVCtoEvaluation(String evaluationName, VariantContext vc, Map<String, VariantContext> comps, EvaluationContext group) {
        if ( vc == null )
            return true;
        
        if ( evaluationName.contains(FILTERED_SET_NAME) && vc.isNotFiltered() )
            return false;

        if ( evaluationName.contains(RETAINED_SET_NAME) && vc.isFiltered() )
            return false;

        boolean vcKnown = vcIsKnown(vc, comps, KNOWN_NAMES);
        if ( evaluationName.contains(KNOWN_SET_NAME) && ! vcKnown )
            return false;
        else if ( evaluationName.contains(NOVEL_SET_NAME) && vcKnown )
            return false;

        if ( group.selectExp != null && ! VariantContextUtils.match(vc, group.selectExp) )
            return false;

        // nothing invalidated our membership in this set
        return true;
    }

    private Map<String, VariantContext> getCompVariantContexts(RefMetaDataTracker tracker, AlignmentContext context) {
        Map<String, VariantContext> comps = new HashMap<String, VariantContext>();

        for ( String compName : compNames ) {
            comps.put(compName, getVariantContext(compName, tracker, context));
        }

        return comps;
    }

    private boolean vcIsKnown(VariantContext vc, Map<String, VariantContext> comps, String[] knownNames ) {
        for ( String knownName : knownNames ) {
            VariantContext known = comps.get(knownName);
            if ( known != null && known.isNotFiltered() && known.getType() == vc.getType() )
                return true;
        }

        return false;
    }

// can't handle this situation
// todo -- warning, this leads to some missing SNPs at complex loci, such as:
// todo -- 591     1       841619  841620  rs4970464       0       -       A       A       -/C/T   genomic mixed   unknown 0       0       near-gene-3     exact   1
// todo -- 591     1       841619  841620  rs62677860      0       +       A       A       C/T     genomic single  unknown 0       0       near-gene-3     exact   1
//
//logger.info(String.format("Ignore second+ events at locus %s in rod %s => rec is %s", context.getLocation(), rodList.getName(), rec));

    private VariantContext getVariantContext(String name, RefMetaDataTracker tracker, AlignmentContext context) {
        if ( tracker != null ) {
            RODRecordList<ReferenceOrderedDatum> rodList = tracker.getTrackData(name, null);
            if ( rodList != null ) {
                for ( ReferenceOrderedDatum rec : rodList.getRecords() ) {
                    if ( rec.getLocation().getStart() == context.getLocation().getStart() ) {
                        VariantContext vc = VariantContextAdaptors.convertToVariantContext(rec);
                        if ( vc != null ) {
                            return vc;
                        }
                    }
                }
            }
        }

        return null;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    // --------------------------------------------------------------------------------------------------------------
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer point, Integer sum) {
        return point + sum;
    }

    public VariantEvaluator getEvalByName(String name, Set<VariantEvaluator> s) {
        for ( VariantEvaluator e : s )
            if ( e.getName().equals(name) )
                return e;
        return null;
    }

    public <T extends Comparable<T>> List<T> sorted(Collection<T> c ) {
        List<T> l = new ArrayList<T>(c);
        Collections.sort(l);
        return l;
    }

    public void onTraversalDone(Integer result) {
        for ( String evalName : variantEvaluationNames ) {
            boolean first = true;
            out.printf("%n%n");
            for ( String contextName : sorted(contexts.keySet()) ) {
                EvaluationContext group = contexts.get(contextName);

                out.printf("%s%n", Utils.dupString('-', 80));
                for ( String evalSubgroupName : sorted(group.keySet()) ) {
                    Set<VariantEvaluator> evalSet = group.get(evalSubgroupName);
                    VariantEvaluator eval = getEvalByName(evalName, evalSet);
                    String keyWord = contextName + "." + evalSubgroupName;
                    if ( first ) {
                        out.printf("%20s\t%40s\t%s%n", evalName, "context", Utils.join("\t\t", eval.getTableHeader()));
                        first = false;
                    }

                    for ( List<String> row : eval.getTableRows() )
                        out.printf("%20s\t%40s\t%s%n", evalName, keyWord, Utils.join("\t\t", row));
                }
            }
        }
    }
}
