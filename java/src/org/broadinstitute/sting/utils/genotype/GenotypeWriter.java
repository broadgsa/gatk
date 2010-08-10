package org.broadinstitute.sting.utils.genotype;

import org.broad.tribble.util.variantcontext.VariantContext;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron, ebanks
 *         <p/>
 *         Class GenotypeWriter
 *         <p/>
 *         The interface for writing genotype calls.
 */
public interface GenotypeWriter {
    /**
     * Add a record, given a variant context, with the genotype fields restricted to what is defined in the header
     * @param vc  the variant context representing the call to add
     * @param refBase This is required for VCF writers, as the VCF format explicitly requires (previous) ref base for an indel. 
     */
    public void add(VariantContext vc, byte refBase);

    /** finish writing, closing any open files. */
    public void close();

}
