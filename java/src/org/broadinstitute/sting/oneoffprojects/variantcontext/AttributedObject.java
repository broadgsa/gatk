package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;


/**
 * @author depristo
 *         <p/>
 *         Class AttributedObject
 *         <p/>
 *         Common functions in VariantContext
 */
public class AttributedObject {
    public static final double NO_NEG_LOG_10PERROR = 0.0;
    private double negLog10PError = NO_NEG_LOG_10PERROR;

    private Map<Object, Object> attributes = new HashMap<Object, Object>();

    public AttributedObject() { }

    public AttributedObject(double negLog10PError) {
        setNegLog10PError(negLog10PError);
    }

    public AttributedObject(Map<? extends Object, ? extends Object> attributes) {
        setAttributes(attributes);
    }

    public AttributedObject(Map<? extends Object, ? extends Object> attributes, double negLog10PError) {
        this(attributes);

        setNegLog10PError(negLog10PError);
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with log error rates
    //
    // ---------------------------------------------------------------------------------------------------------

    public boolean hasNegLog10PError() {
        return getNegLog10PError() == NO_NEG_LOG_10PERROR;
    }

    /**
     * @return the -1 * log10-based error estimate
     */
    public double getNegLog10PError() { return negLog10PError; }

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
    public Map<Object, Object> getAttributes() {
        return attributes;
    }

    // todo -- define common attributes as enum

    public void setAttributes(Map<? extends Object, ? extends Object> map) {
        this.attributes.clear();
        putAttributes(map);
    }

    public void putAttribute(Object key, Object value) {
        putAttribute(key, value, false);
    }

    public void putAttribute(Object key, Object value, boolean allowOverwrites) {
        if ( hasAttribute(key) && ! allowOverwrites )
            throw new StingException("Attempting to overwrite key->value binding: key = " + key + " this = " + this);

        this.attributes.put(key, value);
    }

    public void removeAttribute(Object key) {
        this.attributes.remove(key);
    }

    public void putAttributes(Map<? extends Object, ? extends Object> map) {
        if ( map != null ) {
            for ( Map.Entry<? extends Object, ? extends Object> elt : map.entrySet() ) {
                putAttribute(elt.getKey(), elt.getValue());
            }
        }
    }

    public boolean hasAttribute(Object key) {
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
    public Object getAttribute(Object key) {
        return attributes.get(key);
    }

    public Object getAttribute(Object key, Object defaultValue) {
        if ( hasAttribute(key) )
            return attributes.get(key);
        else
            return defaultValue;
    }

    public AttributedObject getAttributes(Collection<Object> keys) {
        AttributedObject selected = new AttributedObject();

        for ( Object key : keys )
            selected.putAttribute(key, this.getAttribute(key));

        return selected;
    }


    public String getAttributeAsString(Object key)      { return (String)getAttribute(key); }
    public int getAttributeAsInt(Object key)            { return (Integer)getAttribute(key); }
    public double getAttributeAsDouble(Object key)      { return (Double)getAttribute(key); }

    public String getAttributeAsString(Object key, String defaultValue)   { return (String)getAttribute(key, defaultValue); }
    public int getAttributeAsInt(Object key, int defaultValue)            { return (Integer)getAttribute(key, defaultValue); }
    public double getAttributeAsDouble(Object key, double defaultValue)   { return (Double)getAttribute(key, defaultValue); }
}