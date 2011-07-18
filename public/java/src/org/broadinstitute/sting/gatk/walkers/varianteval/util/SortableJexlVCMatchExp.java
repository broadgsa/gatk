package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import org.apache.commons.jexl2.Expression;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

public class SortableJexlVCMatchExp extends VariantContextUtils.JexlVCMatchExp implements Comparable<SortableJexlVCMatchExp> {
    /**
     * Create a new matcher expression with name and JEXL expression exp
     *
     * @param name name
     * @param exp  expression
     */
    public SortableJexlVCMatchExp(String name, Expression exp) {
        super(name, exp);
    }

    public int compareTo(SortableJexlVCMatchExp sortableJexlVCMatchExp) {
        return this.name.compareTo(sortableJexlVCMatchExp.name);
    }
}
