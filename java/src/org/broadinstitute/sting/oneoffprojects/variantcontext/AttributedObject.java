package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;


/**
 * @author ebanks
 *         <p/>
 *         Class VariantContext
 *         <p/>
 *         This class represents a context that unifies one or more variants
 */
public class AttributedObject {
    private Map<Object, Object> attributes = new HashMap<Object, Object>();

    public AttributedObject() {
        ;
    }

    public AttributedObject(Map<Object, Object> attributes) {
        setAttributes(attributes);
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

    public void setAttributes(Map<? extends Object, Object> map) {
        this.attributes.clear();
        putAttributes(attributes);
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

    public void putAttributes(Map<? extends Object, Object> map) {
        for ( Map.Entry<Object, Object> elt : attributes.entrySet() ) {
            putAttribute(elt.getKey(), elt.getValue());
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