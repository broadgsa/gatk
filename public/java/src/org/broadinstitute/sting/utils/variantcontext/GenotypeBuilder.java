/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext;


import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.*;

/**
 * A builder class for genotypes
 *
 * @author Mark DePristo
 */
public final class GenotypeBuilder {
    private String sampleName = null;
    private List<Allele> alleles = null;

    private boolean isPhased = false;
    private int GQ = -1;
    private int DP = -1;
    private int[] AD = null;
    private int[] PL = null;
    private Map<String, Object> extendedAttributes = null;

    private final static Map<String, Object> NO_ATTRIBUTES =
            new HashMap<String, Object>(0);

    /**
     * Create a empty builder.  Both a sampleName and alleles must be provided
     * before trying to make a Genotype from this builder.
     */
    public GenotypeBuilder() {}

    /**
     * Create a builder using sampleName.  Alleles must be provided
     * before trying to make a Genotype from this builder.
     * @param sampleName
     */
    public GenotypeBuilder(final String sampleName) {
        name(sampleName);
    }

    /**
     * Make a builder using sampleName and alleles for starting values
     * @param sampleName
     * @param alleles
     */
    public GenotypeBuilder(final String sampleName, final List<Allele> alleles) {
        name(sampleName);
        alleles(alleles);
    }

    /**
     * Create a new builder starting with the values in Genotype g
     * @param g
     */
    public GenotypeBuilder(final FastGenotype g) {
        copy(g);
    }

    /**
     * Copy all of the values for this builder from Genotype g
     * @param g
     * @return
     */
    public GenotypeBuilder copy(final FastGenotype g) {
        name(g.getSampleName());
        alleles(g.getAlleles());
        phased(g.isPhased());
        GQ(g.getGQ());
        DP(g.getDP());
        AD(g.getAD());
        PL(g.getPL());
        attributes(g.getExtendedAttributes());
        return this;
    }

    /**
     * Reset all of the builder attributes to their defaults.  After this
     * function you must provide sampleName and alleles before trying to
     * make more Genotypes.
     */
    public final void reset() {
        sampleName = null;
        alleles = null;
        isPhased = false;
        GQ = -1;
        DP = -1;
        AD = null;
        PL = null;
        extendedAttributes = null;
    }

    /**
     * Create a new Genotype object using the values set in this builder.
     *
     * After creation the values in this builder can be modified and more Genotypes
     * created, althrough the contents of array values like PL should never be modified
     * inline as they are not copied for efficiency reasons.
     *
     * @return a newly minted Genotype object with values provided from this builder
     */
    @Ensures({"result != null"})
    public FastGenotype make() {
        final Map<String, Object> ea = extendedAttributes == null ? NO_ATTRIBUTES : extendedAttributes;
        return new FastGenotype(sampleName, alleles, isPhased, GQ, DP, AD, PL, ea);
    }

    @Requires({"sampleName != null"})
    @Ensures({"this.sampleName != null"})
    public GenotypeBuilder name(final String sampleName) {
        this.sampleName = sampleName;
        return this;
    }

    @Ensures({"this.alleles != null"})
    public GenotypeBuilder alleles(final List<Allele> alleles) {
        if ( alleles == null )
            this.alleles = Collections.emptyList();
        else
            this.alleles = alleles;
        return this;
    }

    public GenotypeBuilder phased(final boolean phased) {
        isPhased = phased;
        return this;
    }

    @Requires({"GQ >= -1"})
    @Ensures({"this.GQ == GQ"})
    public GenotypeBuilder GQ(final int GQ) {
        this.GQ = GQ;
        return this;
    }

    @Requires({"DP >= -1"})
    @Ensures({"this.DP == DP"})
    public GenotypeBuilder DP(final int DP) {
        this.DP = DP;
        return this;
    }

    @Requires({"AD == null || AD.length > 0"})
    @Ensures({"this.AD == AD"})
    public GenotypeBuilder AD(final int[] AD) {
        this.AD = AD;
        return this;
    }

    @Requires("PL == null || PL.length > 0")
    @Ensures({"this.PL == PL"})
    public GenotypeBuilder PL(final int[] PL) {
        this.PL = PL;
        return this;
    }

    @Requires("attributes != null")
    @Ensures("extendedAttributes != null")
    public GenotypeBuilder attributes(final Map<String, Object> attributes) {
        for ( Map.Entry<String, Object> pair : attributes.entrySet() )
            attribute(pair.getKey(), pair.getValue());
        return this;
    }

    @Requires({"key != null"})
    @Ensures({"extendedAttributes != null", "extendedAttributes.containsKey(key)"})
    public GenotypeBuilder attribute(final String key, final Object value) {
        if ( extendedAttributes == null )
            extendedAttributes = new HashMap<String, Object>(5);
        extendedAttributes.put(key, value);
        return this;
    }
}