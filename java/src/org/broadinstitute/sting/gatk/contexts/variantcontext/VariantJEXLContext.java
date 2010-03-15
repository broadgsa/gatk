/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.apache.commons.jexl.JexlContext;
import org.apache.commons.jexl.JexlHelper;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.util.*;

/**
 *
 * @author aaron
 *
 * Class VariantJEXLContext
 *
 * implements the JEXML context for VariantContext; this saves us from
 * having to generate a JEXML context lookup map everytime we want to evaluate an expression.
 *
 * This is package protected, only classes in variantcontext should have access to it.
 */

class VariantJEXLContext implements JexlContext {
    // our stored variant context
    private final JEXLMap map;

    public VariantJEXLContext(Collection<VariantContextUtils.JexlVCMatchExp> jexl, VariantContext vc) {
        map = new JEXLMap(jexl, vc);
    }

    @Override
    public void setVars(Map map) {
        throw new UnsupportedOperationException("this operation is unsupported");
    }

    @Override
    public Map getVars() {
        return map;
    }
}


/**
 * this is an implementation of a Map of JexlVCMatchExp to true or false values.  It lazy initializes each value
 * as requested to save as much processing time as possible.
 *
 * Compatible with JEXL 1.1 (this code will be easier if we move to 2.0, all of the functionality can go into the
 * JexlContext's get()
 * 
 */

class JEXLMap implements Map<VariantContextUtils.JexlVCMatchExp, Boolean> {
    // our variant context
    private final VariantContext vc;

    // our context
    private JexlContext jContext = null;

    // our mapping from JEXLVCMatchExp to Booleans, which will be set to NULL for previously uncached JexlVCMatchExp
    private final Map<VariantContextUtils.JexlVCMatchExp,Boolean> jexl;


    public JEXLMap(Collection<VariantContextUtils.JexlVCMatchExp> jexlCollection, VariantContext vc) {
        this.vc = vc;
        jexl = new HashMap<VariantContextUtils.JexlVCMatchExp,Boolean>();
        for (VariantContextUtils.JexlVCMatchExp exp: jexlCollection) {
            jexl.put(exp, null);
        }
    }

    /**
     * create the internal JexlContext, only when required.  This code is where new JEXL context variables
     * should get added.
     *
     * @param vc the VariantContext
     *
     */
    private void createContext(VariantContext vc) {
        // create a mapping of what we know about the variant context, its Chromosome, positions, etc.
        Map<String, String> infoMap = new HashMap<String, String>();
        infoMap.put("CHROM", vc.getLocation().getContig());
        infoMap.put("POS", String.valueOf(vc.getLocation().getStart()));
        infoMap.put("TYPE", vc.getType().toString());
        infoMap.put("QUAL", String.valueOf(10 * vc.getNegLog10PError()));

        // add alleles
        infoMap.put("ALLELES", Utils.join(";", vc.getAlleles()));
        infoMap.put("N_ALLELES", String.valueOf(vc.getNAlleles()));

        // add attributes
        addAttributesToMap(infoMap, vc.getAttributes());

        // add filter fields
        infoMap.put("FILTER", String.valueOf(vc.isFiltered() ? "1" : "0"));
        for ( Object filterCode : vc.getFilters() ) {
            infoMap.put(String.valueOf(filterCode), "1");
        }

        // add genotypes
        // todo -- comment this back in when we figure out how to represent it nicely
//        for ( Genotype g : vc.getGenotypes().values() ) {
//            String prefix = g.getSampleName() + ".";
//            addAttributesToMap(infoMap, g.getAttributes(), prefix);
//            infoMap.put(prefix + "GT", g.getGenotypeString());
//        }


        // create the internal context that we can evaluate expressions against
        jContext = JexlHelper.createContext();
        jContext.setVars(infoMap);
    }

    /**
     * @return the size of the internal data structure
     */
    @Override
    public int size() {
        return jexl.size();
    }

    /**
     * @return true if we're empty
     */
    @Override
    public boolean isEmpty() { return this.jexl.isEmpty(); }

    /**
     * do we contain the specified key
     * @param o the key
     * @return true if we have a value for that key
     */
    @Override
    public boolean containsKey(Object o) { return jexl.containsKey(o); }

    @Override
    public Boolean get(Object o) {
        // if we've already determined the value, return it
        if (jexl.containsKey(o) && jexl.get(o) != null) return jexl.get(o);

        // try and cast the expression
        VariantContextUtils.JexlVCMatchExp e = (VariantContextUtils.JexlVCMatchExp) o;
        evaulateExpression(e);
        return jexl.get(e);
    }

    /**
     * get the keyset of map
     * @return a set of keys of type JexlVCMatchExp
     */
    @Override
    public Set<VariantContextUtils.JexlVCMatchExp> keySet() {
        return jexl.keySet();
    }

    /**
     * get all the values of the map.  This is an expensive call, since it evaluates all keys that haven't
     * been evaluated yet.  This is fine if you truely want all the keys, but if you only want a portion, or  know
     * the keys you want, you would be better off using get() to get them by name.
     * @return a collection of boolean values, representing the results of all the variants evaluated
     */
    @Override
    public Collection<Boolean> values() {
        // this is an expensive call
        for (VariantContextUtils.JexlVCMatchExp exp : jexl.keySet())
            if (jexl.get(exp) == null)
                evaulateExpression(exp);
        return jexl.values();
    }

    /**
     * evaulate a JexlVCMatchExp's expression, given the current context (and setup the context if it's null)
     * @param exp the JexlVCMatchExp to evaluate
     */
    private void evaulateExpression(VariantContextUtils.JexlVCMatchExp exp) {
        // if the context is null, we need to create it to evaluate the JEXL expression
        if (this.jContext == null) createContext(vc);
        try {
            jexl.put (exp, (Boolean) exp.exp.evaluate(jContext));
        } catch (Exception e) {
            throw new StingException(e.getMessage());
        }
    }

    /**
     * helper function: adds the list of attributes to the information map we're building
     * @param infoMap the map
     * @param attributes the attributes
     */
    private static void addAttributesToMap(Map<String, String> infoMap, Map<String, ?> attributes ) {
        for (Map.Entry<String, ?> e : attributes.entrySet()) {
            infoMap.put(String.valueOf(e.getKey()), String.valueOf(e.getValue()));
        }
    }

    @Override
    public Boolean put(VariantContextUtils.JexlVCMatchExp jexlVCMatchExp, Boolean aBoolean) {
        return jexl.put(jexlVCMatchExp,aBoolean);
    }

    @Override
    public void putAll(Map<? extends VariantContextUtils.JexlVCMatchExp, ? extends Boolean> map) {
        jexl.putAll(map);
    }

    // //////////////////////////////////////////////////////////////////////////////////////
    // The Following are unsupported at the moment
    // //////////////////////////////////////////////////////////////////////////////////////

    // this doesn't make much sense to implement, boolean doesn't offer too much variety to deal
    // with evaluating every key in the internal map.
    @Override
    public boolean containsValue(Object o) {
        throw new UnsupportedOperationException("containsValue() not supported on a JEXLMap");
    }

    // this doesn't make much sense
    @Override
    public Boolean remove(Object o) {
        throw new UnsupportedOperationException("remove() not supported on a JEXLMap");
    }


    @Override
    public Set<Entry<VariantContextUtils.JexlVCMatchExp, Boolean>> entrySet() {
        throw new UnsupportedOperationException("clear() not supported on a JEXLMap");
    }

    // nope
    @Override
    public void clear() {
        throw new UnsupportedOperationException("clear() not supported on a JEXLMap");
    }
}
