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
import org.broadinstitute.gatk.utils.exceptions.GATKException;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class to encapsulate the raw data for allele-specific classes compatible with the ReducibleAnnotation interface
 * @param <T> the type of raw data to be stored for later annotation calculation
 */
public class AlleleSpecificAnnotationData<T> extends ReducibleAnnotationData<T>{
    final private List<Allele> alleleList;
    private Allele refAllele;

    public AlleleSpecificAnnotationData(final List<Allele> inputAlleles, final String inputData) {
        super(inputData);
        attributeMap = new HashMap<>();
        for(final Allele a : inputAlleles) {
            attributeMap.put(a, null);
        }
        alleleList = inputAlleles;
        for(Allele a : alleleList) {
            if(a.isReference()) {
                refAllele = a;
            }
        }
    }

    @Override
    public List<Allele> getAlleles() {return Collections.unmodifiableList(alleleList);}

    /**
     * Get the reference allele for this allele-specific data.
     * (Used in cases where annotations compare some attribute of the alt alleles to that of the reference.)
     * @return  the reference allele for this data
     */
    public Allele getRefAllele() {return refAllele;}

    public void setAttributeMap(Map<Allele, T> inputMap) {
        super.setAttributeMap(inputMap);
        checkRefAlleles();
    }

    private void checkRefAlleles() {
        boolean foundRef = false;
        for (Allele a : alleleList) {
            if (a.isReference()) {
                if (foundRef)
                    throw new GATKException("ERROR: multiple reference alleles found in annotation data\n");
                foundRef = true;
            }
        }
        if (!foundRef)
            throw new GATKException("ERROR: no reference alleles found in annotation data\n");
    }

    public String makeRawAnnotationString(String printDelim) {
        String annotationString = "";
        for (final Allele current : alleleList) {
            if (!annotationString.isEmpty())
                annotationString += printDelim;
            if(attributeMap.get(current) != null)
                annotationString += attributeMap.get(current).toString();
        }
        return annotationString.replaceAll("[\\[\\]\\s]", "");
    }
}
