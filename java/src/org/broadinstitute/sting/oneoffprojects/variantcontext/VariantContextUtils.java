package org.broadinstitute.sting.oneoffprojects.variantcontext;

import java.util.*;
import org.apache.commons.jexl.*;


public class VariantContextUtils {
    private static final String UNIQUIFIED_SUFFIX = ".unique";

    /**
     * @param other  another variant context
     *
     * throws an exception if there is a collision such that the same sample exists in both contexts
     * @return a context representing the merge of this context and the other provided context
     */
//    public VariantContext merge(VariantContext left, VariantContext other) {
//        return merge(left, other, false);
//    }

    /**
     * @param other            another variant context
     * @param uniquifySamples  if true and there is a collision such that the same sample exists in both contexts,
     *                           the samples will be uniquified(based on their sources);
     *                           otherwise, an exception will be thrown
     *
     * @return a context representing the merge of this context and the other provided context
     */
//    public VariantContext merge(VariantContext left, VariantContext other, boolean uniquifySamples) {
//        // todo -- make functional
//
//        if ( !left.getLocation().equals(other.getLocation()) )
//            throw new IllegalArgumentException("The locations must be identical for two contexts to be merged");
//
//        Set<String> samples = left.getSampleNames();
//        Set<Genotype> Gs = new HashSet<Genotype>(left.getGenotypes().values());
//
//        for ( Genotype g : other.getGenotypes().values() ) {
//            if ( samples.contains(g.getSampleName()) ) {
//                if ( uniquifySamples )
//                    g.setSampleName(g.getSampleName() + UNIQUIFIED_SUFFIX);
//                else
//                    throw new IllegalStateException("The same sample name exists in both contexts when attempting to merge");
//            }
//            Gs.add(g);
//        }
//
//        HashMap<Object, Object> attrs = new HashMap<Object, Object>(left.getAttributes());
//        attrs.putAll(other.getAttributes());
//
//        return new VariantContext(left, Gs, attrs);
//    }


    /**
     * @param subclass  the name of a subclass of variants to select
     *
     * @return a subset of this context which selects based on the given subclass
     */
//    public VariantContextUtils select(String subclass) {
//        HashSet<Genotype> Gs = new HashSet<Genotype>();
//        for ( Genotype g : genotypes ) {
//            if ( g.getAttribute(subclass) != null )
//                Gs.add(g);
//        }
//        return createNewContext(Gs, attributes);
//    }

    /**
     * @param expr  a jexl expression describing how to filter this context
     *
     * @return a subset of this context which is filtered based on the given expression
     */
//    public VariantContextUtils filter(String expr) {
//        HashSet<Genotype> Gs = new HashSet<Genotype>();
//        try {
//            Expression filterExpression = ExpressionFactory.createExpression(expr);
//
//            for ( Genotype g : genotypes ) {
//                JexlContext jContext = JexlHelper.createContext();
//                jContext.setVars(g.getAttributes());
//                if ( (Boolean)filterExpression.evaluate(jContext) )
//                    Gs.add(g);
//            }
//
//        } catch (Exception e) {
//            throw new StingException("JEXL error in VariantContext: " + e.getMessage());
//        }
//
//        return createNewContext(Gs, attributes);
//    }

    /**
     * @return a set of new variant contexts, one for each sample from this context
     */
//    public Set<VariantContextUtils> splitBySample() {
//        Set<VariantContextUtils> contexts = new HashSet<VariantContextUtils>();
//        for ( Genotype g : genotypes ) {
//            HashSet<Genotype> gAsSet = new HashSet<Genotype>();
//            gAsSet.add(g);
//            contexts.add(createNewContext(gAsSet, attributes));
//        }
//        return contexts;
//    }
}