package org.broadinstitute.sting.oneoffprojects.variantcontext.varianteval2;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
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
    protected String[] SELECT_STRINGS = {};
    // todo -- add doc string
    @Argument(shortName="selectName", doc="", required=false)
    protected String[] SELECT_NAMES = {};

    @Argument(shortName="known", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected String[] KNOWN_NAMES = {"dbsnp"};


    private class EvaluationGroup extends HashMap<String, Set<VariantEvaluator>> {
        // useful for typing
    }

    // todo -- generalize to multiple contexts, one for each eval
    private HashMap<String, EvaluationGroup> contexts = new HashMap<String, EvaluationGroup>();
    private List<VariantContextUtils.MatchExp> selectExps = null;

    public void initialize() {
        if ( SELECT_NAMES.length != SELECT_STRINGS.length )
            throw new StingException("Inconsistent number of provided filter names and expressions.");
        Map<String, String> map = new HashMap<String, String>();
        for ( int i = 0; i < SELECT_NAMES.length; i++ ) { map.put(SELECT_NAMES[i], SELECT_STRINGS[i]); }

        selectExps = VariantContextUtils.initializeMatchExps(map);

        // setup contexts
        // todo -- add selects
        contexts.put("eval", createEvaluationGroup());
        contexts.put("eval.filtered", createEvaluationGroup());
    }

    private EvaluationGroup createEvaluationGroup() {
        EvaluationGroup group = new EvaluationGroup();

        for ( String name : Arrays.asList("all", "known", "novel") ) {
            group.put(name, new HashSet<VariantEvaluator>(Arrays.asList(new CountVariants())));
        }

        return group;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //System.out.printf("map at %s with %d skipped%n", context.getLocation(), context.getSkippedBases());

        Map<String, VariantContext> allNamedVCs = getVariantContexts(context, tracker);
        Map<String, VariantContext> evalNamedVCs = getEvalVCs(allNamedVCs);
        Map<String, VariantContext> comps = getComparisonVCs(allNamedVCs);

        for ( VariantEvaluator eval : getAllEvaluations() ) {
            eval.update0(tracker, ref, context);
        }

        if ( evalNamedVCs.size() > 1 ) throw new StingException("VariantEval doesn't yet support for multiple independent eval tracks");
        for ( Map.Entry<String, VariantContext> evalNameVC: evalNamedVCs.entrySet() ) {
            String name = evalNameVC.getKey();
            VariantContext vc = evalNameVC.getValue();
            boolean isKnown = vcIsKnown(vc, allNamedVCs);

            if ( vc.isFiltered() ) name = name + ".filtered";
            EvaluationGroup group = contexts.get(name);

            for ( Set<VariantEvaluator> evaluations : Arrays.asList(group.get("all"), group.get(isKnown ? "known" : "novel")) ) {
                for ( VariantEvaluator evaluation : evaluations ) {
                    switch ( evaluation.getComparisonOrder() ) {
                        case 1:
                            evaluation.update1(vc, tracker, ref, context);
                            break;
                        case 2:
                            for ( VariantContext comp : comps.values() ) {
                                evaluation.update2(vc, comp, tracker, ref, context);
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

    private boolean vcIsKnown(VariantContext vc, Map<String, VariantContext> vcs) {
        for ( VariantContext known : getKnownVCs(vcs).values() ) {
            if ( known.isNotFiltered() && known.getType() == vc.getType() )
                return true;
        }

        return false;
    }


    private Map<String, VariantContext> getEvalVCs(Map<String, VariantContext> vcs) {
        return getVCsStartingWith(vcs, false);
    }

    private Map<String, VariantContext> getComparisonVCs(Map<String, VariantContext> vcs) {
        return getVCsStartingWith(vcs, true);
    }

    private Map<String, VariantContext> getVCsStartingWith(Map<String, VariantContext> vcs, boolean notStartsWith) {
        Map<String, VariantContext> map = new HashMap<String, VariantContext>();

        for ( Map.Entry<String, VariantContext> elt : vcs.entrySet() ) {
            boolean startP = elt.getKey().startsWith("eval");

            if ( (startP && ! notStartsWith) || (!startP && notStartsWith) ) {
                map.put(elt.getKey(), elt.getValue());
            }
        }

        return map;
    }

    private Map<String, VariantContext> getKnownVCs(Map<String, VariantContext> vcs) {
        Map<String, VariantContext> map = new HashMap<String, VariantContext>();

        for ( Map.Entry<String, VariantContext> elt : vcs.entrySet() ) {
            for ( String known1 : KNOWN_NAMES )  {
                if ( elt.getKey().equals(known1) ) {
                    map.put(elt.getKey(), elt.getValue());
                }
            }
        }

        return map;
    }

    private List<String> getAllEvaluationNames() {
        List<String> names = new ArrayList<String>();
        if ( contexts.size() == 0 ) throw new IllegalStateException("Contexts shouldn't be sized 0 when calling getAllEvaluationNames()");

        for ( VariantEvaluator eval : contexts.values().iterator().next().get("all") ) {
            names.add(eval.getName());
        }

        return names;
    }


    private Map<String, VariantContext> getVariantContexts(AlignmentContext context, RefMetaDataTracker tracker) {
        Map<String, VariantContext> map = new HashMap<String, VariantContext>();

        if ( tracker != null ) {
            for ( RODRecordList<ReferenceOrderedDatum> rodList : tracker.getBoundRodTracks() ) {
                boolean alreadyGrabbedOne = false;

                for ( ReferenceOrderedDatum rec : rodList.getRecords() ) {
                    if ( rec.getLocation().getStart() == context.getLocation().getStart() ) {
                        // ignore things that span this location but started earlier
                        if ( alreadyGrabbedOne ) {
                            // can't handle this situation
                            ;
                            //logger.info(String.format("Ignore second+ events at locus %s in rod %s => rec is %s", context.getLocation(), rodList.getName(), rec));
                        } else {
                            VariantContext vc = VariantContextAdaptors.convertToVariantContext(rec);
                            if ( vc != null ) {
                                alreadyGrabbedOne = true;
                                map.put(rec.getName(), vc);
                            }
                        }
                    }
                }
            }
        }

        return map;
    }

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


    public List<VariantEvaluator> getAllEvaluations() {
        List<VariantEvaluator> l = new ArrayList<VariantEvaluator>();

        for ( EvaluationGroup group : contexts.values() ) {
            for ( Set<VariantEvaluator> evals : group.values() ) {
                l.addAll(evals);
            }
        }

        return l;
    }


    public void onTraversalDone(Integer result) {
        for ( String evalName : getAllEvaluationNames() ) {
            boolean first = true;
            for ( Map.Entry<String, EvaluationGroup> elt : contexts.entrySet() ) {
                String contextName = elt.getKey();
                EvaluationGroup group = elt.getValue();

                for ( Map.Entry<String, Set<VariantEvaluator>> namedEvalGroup : group.entrySet() ) {
                    String evalSubgroupName = namedEvalGroup.getKey();
                    VariantEvaluator eval = getEvalByName(evalName, namedEvalGroup.getValue());
                    String keyWord = contextName + "." + evalSubgroupName;
                    if ( first ) {
                        out.printf("%s\t%s\t%s%n", evalName, "context", Utils.join("\t", eval.getTableHeader()));
                        first = false;
                    }

                    for ( List<String> row : eval.getTableRows() )
                        out.printf("%s\t%s\t%s%n", evalName, keyWord, Utils.join("\t", row));
                }
            }
        }
    }
}
