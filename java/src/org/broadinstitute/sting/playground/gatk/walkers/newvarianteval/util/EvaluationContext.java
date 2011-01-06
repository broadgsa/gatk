//package org.broadinstitute.sting.playground.gatk.walkers.varianteval.util;
//
//import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
//import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.VariantEvalWalker;
//import org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.evaluators.VariantEvaluator;
//import org.broadinstitute.sting.utils.Utils;
//import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
//
//import java.lang.reflect.Constructor;
//import java.util.Arrays;
//import java.util.HashSet;
//import java.util.Set;
//
///**
// * Created by IntelliJ IDEA.
// * User: kiran
// * Date: Dec 14, 2010
// * Time: 11:15:32 PM
// * To change this template use File | Settings | File Templates.
// */
//public class EvaluationContext implements Comparable<EvaluationContext> {
//    // useful for typing
//    public String evalTrackName, compTrackName, novelty, filtered, cpgStatus, functionalClass, sample;
//    public boolean enableInterestingSiteCaptures = false;
//    public VariantContextUtils.JexlVCMatchExp selectExp;
//    public Set<VariantEvaluator> evaluations;
//
//    private Set<Class<? extends VariantEvaluator>> evaluationClasses;
//
//    private static String RAW_SET_NAME      = "raw";
//    private static String RETAINED_SET_NAME = "called";
//    private static String FILTERED_SET_NAME = "filtered";
//    private static String ALL_SET_NAME      = "all";
//    private static String KNOWN_SET_NAME    = "known";
//    private static String NOVEL_SET_NAME    = "novel";
//    private final static String CONTEXT_SEPARATOR = "XXX";
//
//    public boolean isIgnoringFilters()      { return filtered.equals(RAW_SET_NAME); }
//    public boolean requiresFiltered()       { return filtered.equals(FILTERED_SET_NAME); }
//    public boolean requiresNotFiltered()    { return filtered.equals(RETAINED_SET_NAME); }
//    public boolean isIgnoringNovelty()      { return novelty.equals(ALL_SET_NAME); }
//    public boolean requiresNovel()          { return novelty.equals(NOVEL_SET_NAME); }
//    public boolean requiresKnown()          { return novelty.equals(KNOWN_SET_NAME); }
//
//    public boolean isSelected() { return selectExp == null; }
//
//    public String getDisplayName() {
//        return getName(CONTEXT_SEPARATOR);
//    }
//
//    public String getJexlName() {
//        return getName(".");
//    }
//
//    private String getName(String separator) {
//        return Utils.join(separator, Arrays.asList(evalTrackName, compTrackName, selectExp == null ? "all" : selectExp.name, filtered, novelty, cpgStatus, functionalClass, sample));
//    }
//
//    public String toString() { return getDisplayName(); }
//
//    public int compareTo(EvaluationContext other) {
//        return this.getDisplayName().compareTo(other.getDisplayName());
//    }
//
//    public EvaluationContext( String evalName, String compName, String novelty, String filtered, String cpgStatus, String functionalClass, String sample, VariantContextUtils.JexlVCMatchExp selectExp, Set<Class<? extends VariantEvaluator>> evaluationClasses ) {
//        this.evalTrackName = evalName;
//        this.compTrackName = compName;
//        this.novelty = novelty;
//        this.filtered = filtered;
//        this.selectExp = selectExp;
//        this.cpgStatus = cpgStatus;
//        this.functionalClass = functionalClass;
//        this.sample = sample;
//        this.enableInterestingSiteCaptures = selectExp == null;
//        this.evaluationClasses = evaluationClasses;
//        this.evaluations = instantiateEvalationsSet();
//    }
//
//    private Set<VariantEvaluator> instantiateEvalationsSet() {
//        Set<VariantEvaluator> evals = new HashSet<VariantEvaluator>();
//        Object[] args = new Object[]{this};
//        Class<?>[] argTypes = new Class<?>[]{VariantEvalWalker.class};
//
//        for ( Class<? extends VariantEvaluator> c : evaluationClasses ) {
//            try {
//                Constructor<? extends VariantEvaluator> constructor = c.getConstructor(argTypes);
//                VariantEvaluator eval = constructor.newInstance(args);
//                evals.add(eval);
//            } catch (Exception e) {
//                throw new DynamicClassResolutionException(c, e);
//            }
//        }
//
//        return evals;
//    }
//}
