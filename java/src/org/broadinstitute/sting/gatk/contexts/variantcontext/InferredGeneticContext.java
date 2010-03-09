package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.broadinstitute.sting.utils.StingException;

import java.util.*;


/**
 * Common utility routines for VariantContext and Genotype
 *
 * @author depristo
 */
final class InferredGeneticContext {
    public static final double NO_NEG_LOG_10PERROR = -1.0;

    private double negLog10PError = NO_NEG_LOG_10PERROR;
    private String name = null;
    private Set<String> filters = new HashSet<String>();
    private Map<String, Object> attributes = new HashMap<String, Object>();

//    public InferredGeneticContext(String name) {
//        this.name = name;
//    }
//
//    public InferredGeneticContext(String name, double negLog10PError) {
//        this(name);
//        setNegLog10PError(negLog10PError);
//    }

    public InferredGeneticContext(String name, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        this.name = name;
        setNegLog10PError(negLog10PError);
        if ( filters != null )
            setFilters(filters);
        if ( attributes != null )
            setAttributes(attributes);
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * Sets the name
     *
     * @param name    the name associated with this information
     */
    public void setName(String name) {
        if ( name == null ) throw new IllegalArgumentException("Name cannot be null " + this);
        this.name = name;
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Filter
    //
    // ---------------------------------------------------------------------------------------------------------

    public Set<String> getFilters() {
        return Collections.unmodifiableSet(filters);
    }

    public boolean isFiltered() {
        return filters.size() > 0;
    }

    public boolean isNotFiltered() {
        return ! isFiltered();
    }

    public void addFilter(String filter) {
        if ( filter == null ) throw new IllegalArgumentException("BUG: Attempting to add null filter " + this);
        if ( getFilters().contains(filter) ) throw new IllegalArgumentException("BUG: Attempting to add duplicate filter " + filter + " at " + this);
        filters.add(filter);
    }

    public void addFilters(Collection<String> filters) {
        if ( filters == null ) throw new IllegalArgumentException("BUG: Attempting to add null filters at" + this);
        for ( String f : filters )
            addFilter(f);
    }

    public void clearFilters() {
        filters.clear();
    }

    public void setFilters(Collection<String> filters) {
        clearFilters();
        addFilters(filters);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with log error rates
    //
    // ---------------------------------------------------------------------------------------------------------

    public boolean hasNegLog10PError() {
        return getNegLog10PError() != NO_NEG_LOG_10PERROR;
    }

    /**
     * @return the -1 * log10-based error estimate
     */
    public double getNegLog10PError() { return negLog10PError; }
    public double getPhredScaledQual() { return getNegLog10PError() * 10; }

    public void setNegLog10PError(double negLog10PError) {
        if ( negLog10PError < 0 && negLog10PError != NO_NEG_LOG_10PERROR ) throw new IllegalArgumentException("BUG: negLog10PError cannot be < than 0 : " + negLog10PError);
        if ( Double.isInfinite(negLog10PError) ) throw new IllegalArgumentException("BUG: negLog10PError should not be Infinity");
        if ( Double.isNaN(negLog10PError) ) throw new IllegalArgumentException("BUG: negLog10PError should not be NaN");

        this.negLog10PError = negLog10PError;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with attributes
    //
    // ---------------------------------------------------------------------------------------------------------
    public void clearAttributes() {
        this.attributes.clear();
    }

    /**
     * @return the attribute map
     */
    public Map<String, Object> getAttributes() {
        return Collections.unmodifiableMap(attributes);
    }

    // todo -- define common attributes as enum

    public void setAttributes(Map<String, ?> map) {
        this.attributes.clear();
        putAttributes(map);
    }

    public void putAttribute(String key, Object value) {
        putAttribute(key, value, false);
    }

    public void putAttribute(String key, Object value, boolean allowOverwrites) {
        if ( hasAttribute(key) && ! allowOverwrites )
            throw new StingException("Attempting to overwrite key->value binding: key = " + key + " this = " + this);

        this.attributes.put(key, value);
    }

    public void removeAttribute(String key) {
        this.attributes.remove(key);
    }

    public void putAttributes(Map<String, ?> map) {
        if ( map != null ) {
            for ( Map.Entry<String, ?> elt : map.entrySet() ) {
                putAttribute(elt.getKey(), elt.getValue());
            }
        }
    }

    public boolean hasAttribute(String key) {
        return attributes.containsKey(key);
    }

    public int getNumAttributes() {
        return attributes.size();
    }

    /**
     * @param key    the attribute key
     *
     * @return the attribute value for the given key (or null if not set)
     */
    public Object getAttribute(String key) {
        return attributes.get(key);
    }

    public Object getAttribute(String key, Object defaultValue) {
        if ( hasAttribute(key) )
            return attributes.get(key);
        else
            return defaultValue;
    }

//    public AttributedObject getAttributes(Collection<Object> keys) {
//        AttributedObject selected = new AttributedObject();
//
//        for ( Object key : keys )
//            selected.putAttribute(key, this.getAttribute(key));
//
//        return selected;
//    }

    public String getAttributeAsString(String key)      { return (String)getAttribute(key); }
    public int getAttributeAsInt(String key)            { return (Integer)getAttribute(key); }
    public double getAttributeAsDouble(String key)      { return (Double)getAttribute(key); }

    public String getAttributeAsString(String key, String defaultValue)   { return (String)getAttribute(key, defaultValue); }
    public int getAttributeAsInt(String key, int defaultValue)            { return (Integer)getAttribute(key, defaultValue); }
    public double getAttributeAsDouble(String key, double defaultValue)   { return (Double)getAttribute(key, defaultValue); }
}