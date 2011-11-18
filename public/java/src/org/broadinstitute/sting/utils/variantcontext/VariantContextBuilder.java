/*
 * Copyright (c) 2011, The Broad Institute
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

import com.google.java.contract.Requires;
import org.broad.tribble.Feature;
import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Builder class for VariantContext
 *
 * @author depristo
 */
public class VariantContextBuilder {
    // required fields
    private String source = null;
    private String contig = null;
    private long start = -1;
    private long stop = -1;
    private Collection<Allele> alleles = null;

    // optional -> these are set to the appropriate default value
    private String ID = VCFConstants.EMPTY_ID_FIELD;
    private GenotypesContext genotypes = GenotypesContext.NO_GENOTYPES;
    private double negLog10PError = VariantContext.NO_NEG_LOG_10PERROR;
    private Set<String> filters = null;
    private Map<String, Object> attributes = null;
    private boolean attributesCanBeModified = false;
    private Byte referenceBaseForIndel = null;
    private boolean genotypesAreUnparsed = false;

    /** enum of what must be validated */
    final private EnumSet<VariantContext.Validation> toValidate = EnumSet.noneOf(VariantContext.Validation.class);

    public VariantContextBuilder() {

    }

    public VariantContextBuilder(String source, String contig, long start, long stop, Collection<Allele> alleles) {
        this.source = source;
        this.contig = contig;
        this.start = start;
        this.stop = stop;
        this.alleles = alleles;
        toValidate.add(VariantContext.Validation.ALLELES);
    }

    /**
     * Returns a new builder based on parent -- the new VC will have all fields initialized
     * to their corresponding values in parent.  This is the best way to create a derived VariantContext
     *
     * @param parent
     */
    public VariantContextBuilder(VariantContext parent) {
        this.alleles = parent.alleles;
        this.attributes = parent.getAttributes();
        this.attributesCanBeModified = false;
        this.contig = parent.contig;
        this.filters = parent.getFiltersMaybeNull();
        this.genotypes = parent.genotypes;
        this.genotypesAreUnparsed = parent.hasAttribute(VariantContext.UNPARSED_GENOTYPE_MAP_KEY);
        this.ID = parent.getID();
        this.negLog10PError = parent.getNegLog10PError();
        this.referenceBaseForIndel = parent.getReferenceBaseForIndel();
        this.source = parent.getSource();
        this.start = parent.getStart();
        this.stop = parent.getEnd();
    }

    @Requires({"alleles != null", "!alleles.isEmpty()"})
    public VariantContextBuilder alleles(final Collection<Allele> alleles) {
        this.alleles = alleles;
        toValidate.add(VariantContext.Validation.ALLELES);
        return this;
    }

    /**
     * Attributes can be null -> meaning there are no attributes.  After
     * calling this routine the builder assumes it can modify the attributes
     * object here, if subsequent calls are made to set attribute values
     * @param attributes
     */
    public VariantContextBuilder attributes(final Map<String, Object> attributes) {
        this.attributes = attributes;
        this.attributesCanBeModified = true;
        return this;
    }

    public VariantContextBuilder attribute(final String key, final Object value) {
        if ( ! attributesCanBeModified ) {
            this.attributesCanBeModified = true;
            this.attributes = new HashMap<String, Object>();
        }
        attributes.put(key, value);
        return this;
    }

    /**
     * filters can be null -> meaning there are no filters
     * @param filters
     */
    public VariantContextBuilder filters(final Set<String> filters) {
        this.filters = filters;
        return this;
    }

    public VariantContextBuilder filters(final String ... filters) {
        filters(new HashSet<String>(Arrays.asList(filters)));
        return this;
    }

    public VariantContextBuilder passFilters() {
        return filters(VariantContext.PASSES_FILTERS);
    }

    public VariantContextBuilder unfiltered() {
        this.filters = null;
        return this;
    }

    /**
     * genotypes can be null -> meaning there are no genotypes
     * @param genotypes
     */
    public VariantContextBuilder genotypes(final GenotypesContext genotypes) {
        this.genotypes = genotypes;
        if ( genotypes != null )
            toValidate.add(VariantContext.Validation.GENOTYPES);
        return this;
    }

    public VariantContextBuilder genotypes(final Collection<Genotype> genotypes) {
        return genotypes(GenotypesContext.copy(genotypes));
    }

    public VariantContextBuilder genotypes(final Genotype ... genotypes) {
        return genotypes(GenotypesContext.copy(Arrays.asList(genotypes)));
    }

    public VariantContextBuilder noGenotypes() {
        this.genotypes = null;
        return this;
    }

    public VariantContextBuilder genotypesAreUnparsed(final boolean genotypesAreUnparsed) {
        this.genotypesAreUnparsed = genotypesAreUnparsed;
        return this;
    }

    @Requires("ID != null")
    public VariantContextBuilder id(final String ID) {
        this.ID = ID;
        return this;
    }

    public VariantContextBuilder noID() {
        return id(VCFConstants.EMPTY_ID_FIELD);
    }

    @Requires("negLog10PError <= 0")
    public VariantContextBuilder negLog10PError(final double negLog10PError) {
        this.negLog10PError = negLog10PError;
        return this;
    }

    /**
     * Null means no refBase is available
     * @param referenceBaseForIndel
     */
    public VariantContextBuilder referenceBaseForIndel(final Byte referenceBaseForIndel) {
        this.referenceBaseForIndel = referenceBaseForIndel;
        toValidate.add(VariantContext.Validation.REF_PADDING);
        return this;
    }

    @Requires("source != null")
    public VariantContextBuilder source(final String source) {
        this.source = source;
        return this;
    }

    @Requires({"contig != null", "start >= 0", "stop >= 0"})
    public VariantContextBuilder loc(final String contig, final long start, final long stop) {
        this.contig = contig;
        this.start = start;
        this.stop = stop;
        toValidate.add(VariantContext.Validation.ALLELES);
        toValidate.add(VariantContext.Validation.REF_PADDING);
        return this;
    }

    @Requires({"contig != null", "start >= 0", "stop >= 0"})
    public VariantContextBuilder chr(final String contig) {
        this.contig = contig;
        return this;
    }

    @Requires({"start >= 0"})
    public VariantContextBuilder start(final long start) {
        this.start = start;
        toValidate.add(VariantContext.Validation.ALLELES);
        toValidate.add(VariantContext.Validation.REF_PADDING);
        return this;
    }

    @Requires({"stop >= 0"})
    public VariantContextBuilder stop(final long stop) {
        this.stop = stop;
        return this;
    }

    public VariantContext make() {
        return new VariantContext(source, ID, contig, start, stop, alleles,
                genotypes, negLog10PError, filters, attributes,
                referenceBaseForIndel, genotypesAreUnparsed, toValidate);
    }
}
