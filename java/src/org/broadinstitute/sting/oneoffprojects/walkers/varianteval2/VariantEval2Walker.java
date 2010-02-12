package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.MutableVariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.apache.log4j.Logger;

import java.util.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.io.File;

// todo -- evalations should support comment lines
// todo -- add Mendelian variable explanations (nDeNovo and nMissingTransmissions)

//
// todo -- write a simple column table system and have the evaluators return this instead of the list<list<string>> objects
//

// todo -- site frequency spectrum eval (freq. of variants in eval as a function of their AC and AN numbers)
// todo -- multiple sample concordance tool (genotypes in eval vs. genotypes in truth)
// todo -- allele freqeuncy discovery tool (FREQ in true vs. discovery counts in eval).  Needs to process subset of samples in true (pools)
// todo -- clustered SNP counter
// todo -- HWEs
// todo -- Validation data analysis from VE1?  What is it and should we transition it over?
// todo -- indel metrics [count of sizes in/del should be in CountVariants]

//
// todo -- Performance
// todo -- create JEXL context implementing object that simply looks up values for JEXL evaluations.  Throws error for unknown fields
//

//
// todo -- port over SNP density evaluator.
// todo -- make it work with intervals correctly
//

// todo -- counts of snps per target [target name, gene, etc]

// todo -- add subgroup of known variants as to those at hapmap sites [it's in the dbSNP record]

// todo -- deal with performance issues with variant contexts

//
// Todo -- should really include argument parsing @annotations from subclass in this walker.  Very
// todo -- useful general capability.  Right now you need to add arguments to VariantEval2 to handle new
// todo -- evaluation arguments (which is better than passing a string!)
//

//
// todo -- the whole organization only supports a single eval x comp evaluation.  We need to instantiate
// todo -- new contexts for each comparison object too!  The output table should be clear as to what the "comp"
// todo -- variable is in the analysis
//

//
// todo -- write or find a simple way to organize the table like output of variant eval 2.  A generic table of strings?
//

// todo Extend VariantEval, our general-purpose tool for SNP evaluation, to differentiate Ti/Tv at CpG islands and also
// todo classify (and count) variants into coding, non-coding, synonomous/non-symonomous, 2/4 fold degenerate sites, etc.
// todo Assume that the incoming VCF has the annotations (you don't need to do this) but VE2 should split up results by
// todo these catogies automatically (using the default selects)
//
// todo -- We agreed to report two standard values for variant evaluation from Êhere out. ÊOne, we will continue to report
// todo -- the dbSNP 129 rate. Additionally, we will start to report the % of variants found that have already been seen in
// todo -- 1000 Genomes. ÊThis should be implemented as another standard comp_1kg binding, pointing to only variants
// todo -- discovered and released by 1KG.  Might need to make this data set ourselves and keep it in GATK/data like
// todo -- dbsnp rod
//
// todo -- aux. plotting routines for VE2
//
// todo -- Provide separate dbsnp rates for het only calls and any call where there is at least one hom-var genotype,
// todo -- since hets are much more likely to be errors

// todo -- Add Heng's hom run metrics


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

    @Argument(shortName="outputVCF", fullName="InterestingSitesVCF", doc="If provided, interesting sites emitted to this vcf and the INFO field annotated as to why they are interesting", required=false)
    protected String outputVCF = null;

    /** Right now we will only be looking at SNPS */
    EnumSet<VariantContext.Type> ALLOW_VARIANT_CONTEXT_TYPES = EnumSet.of(VariantContext.Type.SNP);

    // --------------------------------------------------------------------------------------------------------------
    //
    // private walker data
    //
    // --------------------------------------------------------------------------------------------------------------

    /** private class holding all of the information about a single evaluation group (e.g., for eval ROD) */
    private class EvaluationContext implements Comparable<EvaluationContext> {
        // useful for typing
        public String evalTrackName, compTrackName, novelty, filtered;
        public boolean enableInterestingSiteCaptures = false;
        VariantContextUtils.JexlVCMatchExp selectExp;
        Set<VariantEvaluator> evaluations;

        public boolean isIgnoringFilters()      { return filtered.equals(RAW_SET_NAME); }
        public boolean requiresFiltered()       { return filtered.equals(FILTERED_SET_NAME); }
        public boolean requiresNotFiltered()    { return filtered.equals(RETAINED_SET_NAME); }
        public boolean isIgnoringNovelty()      { return novelty.equals(ALL_SET_NAME); }
        public boolean requiresNovel()          { return novelty.equals(NOVEL_SET_NAME); }
        public boolean requiresKnown()          { return novelty.equals(KNOWN_SET_NAME); }

        public boolean isSelected() { return selectExp == null; }

        public String getDisplayName() {
            return Utils.join(".", Arrays.asList(evalTrackName, compTrackName, selectExp == null ? "all" : selectExp.name, filtered, novelty));
        }

        public int compareTo(EvaluationContext other) {
            return this.getDisplayName().compareTo(other.getDisplayName());
        }

        public EvaluationContext( String evalName, String compName, String novelty, String filtered, VariantContextUtils.JexlVCMatchExp selectExp ) {
            this.evalTrackName = evalName;
            this.compTrackName = compName;
            this.novelty = novelty;
            this.filtered = filtered;
            this.selectExp = selectExp;
            this.enableInterestingSiteCaptures = selectExp == null;
            this.evaluations = instantiateEvalationsSet();
        }
    }

    private List<EvaluationContext> contexts = null;

    // lists of all comp and eval ROD track names
    private Set<String> compNames = new HashSet<String>();
    private Set<String> evalNames = new HashSet<String>();

    private List<String> variantEvaluationNames = new ArrayList<String>();

    private static String RAW_SET_NAME      = "raw";
    private static String RETAINED_SET_NAME = "called";
    private static String FILTERED_SET_NAME = "filtered";
    private static String ALL_SET_NAME      = "all";
    private static String KNOWN_SET_NAME    = "known";
    private static String NOVEL_SET_NAME    = "novel";

    private final static String CONTEXT_HEADER = "eval.comp.select.filter.novelty";
    private final static int N_CONTEXT_NAME_PARTS = CONTEXT_HEADER.split("\\.").length;
    private static int[] nameSizes = new int[N_CONTEXT_NAME_PARTS];
    static {
        int i = 0;
        for ( String elt : CONTEXT_HEADER.split("\\.") )
            nameSizes[i++] = elt.length();
    }

    // Dynamically determined variantEvaluation classes
    private List<Class<? extends VariantEvaluator>> evaluationClasses = null;

    /** output writer for interesting sites */
    private VCFWriter writer = null;
    private boolean wroteHeader = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    // --------------------------------------------------------------------------------------------------------------

    public void initialize() {
        determineAllEvalations();
        List<VariantContextUtils.JexlVCMatchExp> selectExps = VariantContextUtils.initializeMatchExps(SELECT_NAMES, SELECT_EXPS);

        for ( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if ( d.getName().startsWith("eval") ) {
                evalNames.add(d.getName());
            } else if ( d.getName().startsWith("dbsnp") || d.getName().startsWith("hapmap") || d.getName().startsWith("comp") ) {
                compNames.add(d.getName());
            } else {
                logger.info("Not evaluating ROD binding " + d.getName());
            }
        }

        contexts = initializeEvaluationContexts(evalNames, compNames, selectExps);
        determineContextNamePartSizes();

        if ( outputVCF != null )
            writer = new VCFWriter(new File(outputVCF));
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

    private <T> List<T> append(List<T> selectExps, T elt) {
        List<T> l = new ArrayList<T>(selectExps);
        l.add(elt);
        return l;
    }

    private List<EvaluationContext> initializeEvaluationContexts(Set<String> evalNames, Set<String> compNames, List<VariantContextUtils.JexlVCMatchExp> selectExps) {
        List<EvaluationContext> contexts = new ArrayList<EvaluationContext>();

        selectExps = append(selectExps, null);
        for ( String evalName : evalNames ) {
            for ( String compName : compNames ) {
                for ( VariantContextUtils.JexlVCMatchExp e : selectExps ) {
                    for ( String filteredName : Arrays.asList(RAW_SET_NAME, RETAINED_SET_NAME, FILTERED_SET_NAME) ) {
                        for ( String novelty : Arrays.asList(ALL_SET_NAME, KNOWN_SET_NAME, NOVEL_SET_NAME) ) {
                            EvaluationContext context = new EvaluationContext(evalName, compName, novelty, filteredName, e);
                            contexts.add(context);
                        }
                    }
                }
            }
        }

        Collections.sort(contexts);
        return contexts;
    }

    private Set<VariantEvaluator> instantiateEvalationsSet() {
        Set<VariantEvaluator> evals = new HashSet<VariantEvaluator>();
        Object[] args = new Object[]{this};
        Class[] argTypes = new Class[]{this.getClass()};

        for ( Class c : evaluationClasses ) {
            try {
                Constructor constructor = c.getConstructor(argTypes);
                VariantEvaluator eval = (VariantEvaluator)constructor.newInstance(args);
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



    private boolean captureInterestingSitesOfEvalSet(EvaluationContext group) {
        //System.out.printf("checking %s%n", name);
        return group.requiresNotFiltered() && group.isIgnoringNovelty();
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    // --------------------------------------------------------------------------------------------------------------

    // todo -- call a single function to build a map from track name -> variant context / null for all
    //      -- eval + comp names.  Use this data structure to get data throughout rest of the loops here
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //System.out.printf("map at %s with %d skipped%n", context.getLocation(), context.getSkippedBases());

        if ( ref == null )
            return 0;

        Map<String, VariantContext> vcs = getVariantContexts(tracker, context);
        //Collection<VariantContext> comps = getCompVariantContexts(tracker, context);

        // to enable walking over pairs where eval or comps have no elements
        for ( EvaluationContext group : contexts ) {
            VariantContext vc = vcs.get(group.evalTrackName);

            //logger.debug(String.format("Updating %s of %s with variant", group.name, vc));
            Set<VariantEvaluator> evaluations = group.evaluations;
            boolean evalWantsVC = applyVCtoEvaluation(vc, vcs, group);
            List<String> interestingReasons = new ArrayList<String>();

            for ( VariantEvaluator evaluation : evaluations ) {
                if ( evaluation.enabled() ) {
                    // we always call update0 in case the evaluation tracks things like number of bases covered
                    evaluation.update0(tracker, ref, context);

                    // now call the single or paired update function
                    switch ( evaluation.getComparisonOrder() ) {
                        case 1:
                            if ( evalWantsVC && vc != null ) {
                                String interesting = evaluation.update1(vc, tracker, ref, context);
                                if ( interesting != null ) interestingReasons.add(interesting);
                            }
                            break;
                        case 2:
                            VariantContext comp = vcs.get(group.compTrackName);
                            String interesting = evaluation.update2( evalWantsVC ? vc : null, comp, tracker, ref, context );
                            if ( interesting != null ) interestingReasons.add(interesting);
                            break;
                        default:
                            throw new StingException("BUG: Unexpected evaluation order " + evaluation);
                    }
                }
            }

            if ( group.enableInterestingSiteCaptures && captureInterestingSitesOfEvalSet(group) )
                writeInterestingSite(interestingReasons, vc);
        }

        return 0;
    }

    private void writeInterestingSite(List<String> interestingReasons, VariantContext vc) {
        if ( writer != null && interestingReasons.size() > 0 ) {
            MutableVariantContext mvc = new MutableVariantContext(vc);

            for ( String why : interestingReasons ) {
                String key, value;
                String[] parts = why.split("=");

                switch ( parts.length ) {
                    case 1:
                        key = parts[0];
                        value = "1";
                        break;
                    case 2:
                        key = parts[0];
                        value = parts[1];
                        break;
                    default:
                        throw new IllegalStateException("BUG: saw a interesting site reason sting with multiple = signs " + why);
                }

                mvc.putAttribute(key, value);
            }

            if ( ! wroteHeader ) {
                writer.writeHeader(VariantContextAdaptors.createVCFHeader(null, vc));
                wroteHeader = true;
            }

            writer.addRecord(VariantContextAdaptors.toVCF(mvc));
            //interestingReasons.clear();
        }
    }

    private boolean applyVCtoEvaluation(VariantContext vc, Map<String, VariantContext> vcs, EvaluationContext group) {
        if ( vc == null )
            return true;

        if ( group.requiresFiltered() && vc.isNotFiltered() )
            return false;

        if ( group.requiresNotFiltered() && vc.isFiltered() )
            return false;

        boolean vcKnown = vcIsKnown(vc, vcs, KNOWN_NAMES);
        if ( group.requiresKnown() && ! vcKnown )
            return false;
        else if ( group.requiresNovel() && vcKnown )
            return false;

        if ( group.selectExp != null && ! VariantContextUtils.match(vc, group.selectExp) )
            return false;

        // nothing invalidated our membership in this set
        return true;
    }

    private boolean vcIsKnown(VariantContext vc, Map<String, VariantContext> vcs, String[] knownNames ) {
        for ( String knownName : knownNames ) {
            VariantContext known = vcs.get(knownName);
            if ( known != null && known.isNotFiltered() && known.getType() == vc.getType() ) {
                return true;
            }
        }

        return false;
    }

// can't handle this situation
// todo -- warning, this leads to some missing SNPs at complex loci, such as:
// todo -- 591     1       841619  841620  rs4970464       0       -       A       A       -/C/T   genomic mixed   unknown 0       0       near-gene-3     exact   1
// todo -- 591     1       841619  841620  rs62677860      0       +       A       A       C/T     genomic single  unknown 0       0       near-gene-3     exact   1
//
//logger.info(String.format("Ignore second+ events at locus %s in rod %s => rec is %s", context.getLocation(), rodList.getName(), rec));

    private Map<String, VariantContext> getVariantContexts(RefMetaDataTracker tracker, AlignmentContext context) {
        // todo -- we need to deal with dbSNP where there can be multiple records at the same start site.  A potential solution is to
        // todo -- allow the variant evaluation to specify the type of variants it wants to see and only take the first such record at a site
        Map<String, VariantContext> bindings = new HashMap<String, VariantContext>();
        bindVariantContexts(bindings, evalNames, tracker, context);
        bindVariantContexts(bindings, compNames, tracker, context);
        return bindings;
    }

    private void bindVariantContexts(Map<String, VariantContext> map, Collection<String> names,
                                     RefMetaDataTracker tracker, AlignmentContext context ) {
        for ( String name : names ) {
            Collection<VariantContext> contexts = tracker.getVariantContexts(name, ALLOW_VARIANT_CONTEXT_TYPES, context.getLocation(), true, true);
            if ( context.size() > 1 )
                throw new StingException("Found multiple variant contexts at " + context.getLocation());
            map.put(name, contexts.size() == 1 ? contexts.iterator().next() : null);
        }
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

    private void determineContextNamePartSizes() {
        for ( EvaluationContext group : contexts ) {
            String[] parts = group.getDisplayName().split("\\.");
            if ( parts.length != N_CONTEXT_NAME_PARTS ) {
                throw new StingException("Unexpected number of eval name parts " + group.getDisplayName() + " length = " + parts.length + ", expected " + N_CONTEXT_NAME_PARTS);
            } else {
                for ( int i = 0; i < parts.length; i++ )
                    nameSizes[i] = Math.max(nameSizes[i], parts[i].length());
            }
        }
    }

    private String formatKeyword(String keyWord) {
        //System.out.printf("keyword %s%n", keyWord);

        StringBuilder s = new StringBuilder();
        int i = 0;
        for ( String part : keyWord.split("\\.") ) {
            //System.out.printf("part %s %d%n", part, nameSizes[i]);
            s.append(String.format("%" + nameSizes[i] + "s ", part));
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
            String lastEvalTrack = null;
            for ( EvaluationContext group : contexts ) {
                if ( lastEvalTrack == null || ! lastEvalTrack.equals(group.evalTrackName) ) {
                    out.printf("%s%n", Utils.dupString('-', 80));
                    lastEvalTrack = group.evalTrackName;
                }

                VariantEvaluator eval = getEvalByName(evalName, group.evaluations);
                String keyWord = group.getDisplayName();

                if ( eval.enabled() ) {
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
