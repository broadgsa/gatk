/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;

import java.util.*;

/**
 * A class to encapsulate the raw data for classes compatible with the ReducibleAnnotation interface
 */
public class ReducibleAnnotationData<T> {
    protected String rawData;
    protected Map<Allele, T> attributeMap;

    /**
     * Create a new ReducibleAnnotationData using the raw data string from a VCF
     * @param inputData the raw data as read in from a VCF
     */
    public ReducibleAnnotationData(final String inputData) {
        rawData = inputData; attributeMap = new HashMap<>();
        attributeMap.put(Allele.NO_CALL, null);
    }

    /**
     *
     * @return the string of raw data as represented in the VCF
     */
    public String getRawData() {return rawData;}

    /**
     * Note: parent class ReducibleAnnotationData is non-allele specific and stores all values with the no-call allele
     * @return the list of alleles for which we have raw annotation data
     */
    public List<Allele> getAlleles() {
        List ret = new ArrayList<Allele>();
        ret.addAll(attributeMap.keySet());
        return ret;
    }

    /**
     *
     * @param key   the allele of interest
     * @return  do we have data for the allele of interest?
     */
    public boolean hasAttribute(Allele key) {
        return attributeMap.containsKey(key);
    }

    /**
     *
     * @param key the allele of interest
     * @return  data for the allele of interest
     */
    public T getAttribute(Allele key) {
        return attributeMap.get(key);
    }

    /**
     *
     * @param key   the allele of interest
     * @param value raw data corresponding to the allele of interest
     */
    public void putAttribute(Allele key, T value) {
        attributeMap.put(key, value);
    }

    /**
     * Assign all of the per-allele raw data at once
     * @param inputMap  the pre-calculated per-allele data
     */
    public void setAttributeMap(Map<Allele, T> inputMap) {
        attributeMap = inputMap;
    }

    /**
     * Get the stored raw per-allele data
     * @return
     */
    public Map<Allele, T> getAttributeMap() {return Collections.unmodifiableMap(attributeMap);}

}
