package org.broadinstitute.sting.utils.variantcontext;


import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;

import java.util.*;


/**
 * Common utility routines for VariantContext and Genotype
 *
 * @author depristo
 */
public final class InferredGeneticContext {
    public static final double NO_NEG_LOG_10PERROR = -1.0;

    private static Set<String> NO_FILTERS = Collections.unmodifiableSet(new HashSet<String>());
    private static Map<String, Object> NO_ATTRIBUTES = Collections.unmodifiableMap(new HashMap<String, Object>());

    private double negLog10PError = NO_NEG_LOG_10PERROR;
    private String name = null;
    private Set<String> filters = NO_FILTERS;
    private Map<String, Object> attributes = NO_ATTRIBUTES;

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
        if ( filters == NO_FILTERS ) // immutable -> mutable
            filters = new HashSet<String>(filters);

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
        filters = new HashSet<String>();
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
        attributes = new HashMap<String, Object>();
    }

    /**
     * @return the attribute map
     */
    public Map<String, Object> getAttributes() {
        return Collections.unmodifiableMap(attributes);
    }

    // todo -- define common attributes as enum

    public void setAttributes(Map<String, ?> map) {
        clearAttributes();
        putAttributes(map);
    }

    public void putAttribute(String key, Object value) {
        putAttribute(key, value, false);
    }

    public void putAttribute(String key, Object value, boolean allowOverwrites) {
        if ( ! allowOverwrites && hasAttribute(key) )
            throw new IllegalStateException("Attempting to overwrite key->value binding: key = " + key + " this = " + this);

        if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
            attributes = new HashMap<String, Object>();
        
        attributes.put(key, value);
    }

    public void removeAttribute(String key) {
        if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
            attributes = new HashMap<String, Object>();
        attributes.remove(key);
    }

    public void putAttributes(Map<String, ?> map) {
        if ( map != null ) {
            // for efficiency, we can skip the validation if the map is empty
            if ( attributes.size() == 0 ) {
                if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
                    attributes = new HashMap<String, Object>();
                attributes.putAll(map);
            } else {
                for ( Map.Entry<String, ?> elt : map.entrySet() ) {
                    putAttribute(elt.getKey(), elt.getValue(), false);
                }
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

    public String getAttributeAsString(String key, String defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof String ) return (String)x;
        return String.valueOf(x); // throws an exception if this isn't a string
    }

    public int getAttributeAsInt(String key, int defaultValue) {
        Object x = getAttribute(key);
        if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
        if ( x instanceof Integer ) return (Integer)x;
        return Integer.valueOf((String)x); // throws an exception if this isn't a string
    }

    public double getAttributeAsDouble(String key, double defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Double ) return (Double)x;
        return Double.valueOf((String)x); // throws an exception if this isn't a string
    }

    public boolean getAttributeAsBoolean(String key, boolean defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Boolean ) return (Boolean)x;
        return Boolean.valueOf((String)x); // throws an exception if this isn't a string
    }

//    public String getAttributeAsString(String key)      { return (String.valueOf(getAttribute(key))); } // **NOTE**: will turn a null Object into the String "null"
//    public int getAttributeAsInt(String key)            { Object x = getAttribute(key); return x instanceof Integer ? (Integer)x : Integer.valueOf((String)x); }
//    public double getAttributeAsDouble(String key)      { Object x = getAttribute(key); return x instanceof Double ? (Double)x : Double.valueOf((String)x); }
//    public boolean getAttributeAsBoolean(String key)      { Object x = getAttribute(key); return x instanceof Boolean ? (Boolean)x : Boolean.valueOf((String)x); }
//    public Integer getAttributeAsIntegerNoException(String key)  { try {return getAttributeAsInt(key);} catch (Exception e) {return null;} }
//    public Double getAttributeAsDoubleNoException(String key)    { try {return getAttributeAsDouble(key);} catch (Exception e) {return null;} }
//    public String getAttributeAsStringNoException(String key)    { if (getAttribute(key) == null) return null; return getAttributeAsString(key); }
//    public Boolean getAttributeAsBooleanNoException(String key)  { try {return getAttributeAsBoolean(key);} catch (Exception e) {return null;} }
}