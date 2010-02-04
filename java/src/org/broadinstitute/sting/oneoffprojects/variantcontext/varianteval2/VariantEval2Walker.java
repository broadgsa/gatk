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
import org.apache.log4j.Logger;

import java.util.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;


//
// todo -- write a simple column table system and have the evaluators return this instead of the list<list<string>> objects
// todo -- site frequency spectrum eval (freq. of variants in eval as a function of their AC and AN numbers)
// todo -- multiple sample concordance tool (genotypes in eval vs. genotypes in truth)
// todo -- allele freqeuncy discovery tool (FREQ in true vs. discovery counts in eval).  Needs to process subset of samples in true (pools)
// todo -- clustered SNP counter
// todo -- HWEs
// todo -- Validation data analysis from VE1?  What is it and should we transition it over?
// todo -- indel metrics [count of sizes in/del should be in CountVariants]
// todo -- create JEXL context implementing object that simply looks up values for JEXL evaluations.  Throws error for unknown fields
//

//
// Todo -- should really include argument parsing @annotations from subclass in this walker.  Very
// todo -- useful general capability.  Right now you need to add arguments to VariantEval2 to handle new
// todo -- evaluation arguments (which is better than passing a string!)
//

//
// todo -- the whole organization only supports a single eval x comp evaluation.  We need to instantiate
// todo -- new contexts for each comparison object too!
//

/**
 * Test routine for new VariantContext object
 */
public class VariantEval2Walker extends RodWalker<Integer, Integer> {
    // --------------------------------------------------------------------------------------------------------------
    //
    // walker arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    // todo -- add doc string
    @Argument(shortName="select", doc="", required=false)
    protected String[] SELECT_EXPS = {"QUAL > 500.0", "HARD_TO_VALIDATE==1", "GATK_STANDARD==1"};

    // todo -- add doc string
    @Argument(shortName="selectName", doc="", required=false)
    protected String[] SELECT_NAMES = {"q500plus", "low_mapq", "gatk_std_filters"};

    @Argument(shortName="known", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected String[] KNOWN_NAMES = {"dbsnp"};

    //
    // Arguments for Mendelian Violation calculations
    //
    @Argument(shortName="family", doc="If provided, genotypes in will be examined for mendelian violations: this argument is a string formatted as dad+mom=child where these parameters determine which sample names are examined", required=false)
    protected String FAMILY_STRUCTURE;

    @Argument(shortName="MVQ", fullName="MendelianViolationQualThreshold", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation", required=false)
    protected double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 50;

    @Argument(shortName="PI", fullName="PrintInterestingSites", doc="If provided, interesting sites in the unselected, called set will be printed", required=false)
    protected boolean PRINT_INTERESTING_SITES = false;

//    @Argument(shortName="PIA", fullName="PrintAllInterestingSites", doc="If provided, interesting sites will be printed for all sets.  Verbose", required=false)
//    protected boolean PRINT_ALL_INTERESTING_SITES = false;
    
    // --------------------------------------------------------------------------------------------------------------
    //
    // private walker data
    //
    // --------------------------------------------------------------------------------------------------------------

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
                addNewContext(d.getName(), d.getName() + ".all", null);
            } else if ( d.getName().startsWith("dbsnp") || d.getName().startsWith("hapmap") || d.getName().startsWith("comp") ) {
                compNames.add(d.getName());
            } else {
                logger.info("Not evaluating ROD binding " + d.getName());
            }
        }

        determineContextNamePartSizes();
    }

    private void determineAllEvalations() {
        evaluationClasses = PackageUtils.getClassesImplementingInterface(VariantEvaluator.class);
        for ( VariantEvaluator e : instantiateEvalationsSet(false, null) ) {
            // for collecting purposes
            variantEvaluationNames.add(e.getName());
            logger.debug("Including VariantEvaluator " + e.getName() + " of class " + e.getClass());
        }

        Collections.sort(variantEvaluationNames);
    }

    private Set<VariantEvaluator> instantiateEvalationsSet(boolean baseline, String name) {
        Set<VariantEvaluator> evals = new HashSet<VariantEvaluator>();
        Object[] args = new Object[]{this};
        Class[] argTypes = new Class[]{this.getClass()};

        for ( Class c : evaluationClasses ) {
            try {
                Constructor constructor = c.getConstructor(argTypes);
                VariantEvaluator eval = (VariantEvaluator)constructor.newInstance(args);
                if ( baseline ) eval.printInterestingSites(name);
                evals.add(eval);
            } catch (InstantiationException e) {
                throw new StingException(String.format("Cannot instantiate class '%s': must be concrete class", c.getSimpleName()));
            } catch (NoSuchMethodException e) {
                throw new StingException(String.format("Cannot find expected constructor for class '%s': must have constructor accepting a single VariantEval2Walker object", c.getSimpleName()));
            } catch (IllegalAccessException e) {
                throw new StingException(String.format("Cannot instantiate class '%s':", c.getSimpleName()));
            } catch (InvocationTargetException e) {
                throw new StingException(String.format("Cannot instantiate class '%s':", c.getSimpleName()));
            }
        }

        return evals;
    }

    private void addNewContext(String trackName, String contextName, VariantContextUtils.MatchExp selectExp) {
        EvaluationContext group = new EvaluationContext(trackName, contextName, selectExp);

        for ( String filteredName : Arrays.asList(RAW_SET_NAME, RETAINED_SET_NAME, FILTERED_SET_NAME) ) {
            for ( String subname : Arrays.asList(ALL_SET_NAME, KNOWN_SET_NAME, NOVEL_SET_NAME) ) {
                String name = subname + "." + filteredName;
                group.put(name, instantiateEvalationsSet(subname == ALL_SET_NAME && filteredName == RETAINED_SET_NAME, trackName + "." + (selectExp == null ? "all" : selectExp.name) + "." + name));
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
                            if ( evalWantsVC && vc != null )
                                evaluation.update1(vc, tracker, ref, context);
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

        // todo -- remove me when the loop works correctly for comparisons of eval x comp for each comp
        if ( comps.size() > 1 ) throw new StingException("VariantEval2 currently only supports comparisons of N eval tracks vs. a single comparison track.  Yes, I know...");
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
                        VariantContext vc = VariantContextAdaptors.convertToVariantContext(name, rec);
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

    private <T extends Comparable<T>> List<T> sorted(Collection<T> c ) {
        List<T> l = new ArrayList<T>(c);
        Collections.sort(l);
        return l;
    }

    private final static String CONTEXT_HEADER = "track.subset.novelty.filter";
    private final static int N_CONTEXT_NAME_PARTS = CONTEXT_HEADER.split("\\.").length;
    private static int[] nameSizes = new int[N_CONTEXT_NAME_PARTS];
    static {
        int i = 0;
        for ( String elt : CONTEXT_HEADER.split("\\.") )
            nameSizes[i++] = elt.length();
    }

    private void determineContextNamePartSizes() {
        for ( String contextName : sorted(contexts.keySet()) ) {
            EvaluationContext group = contexts.get(contextName);
            for ( String evalSubgroupName : sorted(group.keySet()) ) {
                String keyWord = contextName + "." + evalSubgroupName;
                String[] parts = keyWord.split("\\.");
                if ( parts.length != N_CONTEXT_NAME_PARTS ) {
                    throw new StingException("Unexpected number of eval name parts " + keyWord + " length = " + parts.length + ", expected " + N_CONTEXT_NAME_PARTS);
                } else {
                    for ( int i = 0; i < parts.length; i++ )
                        nameSizes[i] = Math.max(nameSizes[i], parts[i].length());
                }
            }
        }
    }
    
    private String formatKeyword(String keyWord) {
        //System.out.printf("keyword %s%n", keyWord);

        StringBuilder s = new StringBuilder();
        int i = 0;
        for ( String part : keyWord.split("\\.") ) {
            //System.out.printf("part %s %d%n", part, nameSizes[i]);
            s.append(String.format("%"+nameSizes[i]+"s ", part));
            i++;
        }

        return s.toString();
    }

    public void onTraversalDone(Integer result) {
        // todo -- this really needs to be pretty printed; use some kind of table organization
        // todo -- so that we can load up all of the data in one place, analyze the widths of the columns
        // todo -- and pretty print it
        for ( String evalName : variantEvaluationNames ) {
            boolean first = true;
            out.printf("%n%n");
            // todo -- show that comp is dbsnp, etc. is columns
            for ( String contextName : sorted(contexts.keySet()) ) {
                EvaluationContext group = contexts.get(contextName);

                out.printf("%s%n", Utils.dupString('-', 80));
                for ( String evalSubgroupName : sorted(group.keySet()) ) {
                    Set<VariantEvaluator> evalSet = group.get(evalSubgroupName);
                    VariantEvaluator eval = getEvalByName(evalName, evalSet);
                    String keyWord = contextName + "." + evalSubgroupName;

                    if ( first ) {
                        out.printf("%20s %s %s%n", evalName, formatKeyword(CONTEXT_HEADER), Utils.join("\t", eval.getTableHeader()));
                        first = false;
                    }

                    for ( List<String> row : eval.getTableRows() )
                        out.printf("%20s %s %s%n", evalName, formatKeyword(keyWord), Utils.join("\t", row));
                }
            }
        }
    }

    protected Logger getLogger() { return logger; }
}
