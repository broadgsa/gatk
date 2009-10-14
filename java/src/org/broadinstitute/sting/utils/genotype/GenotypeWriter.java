package org.broadinstitute.sting.utils.genotype;

import java.util.List;

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
 * @author aaron
 *         <p/>
 *         Class Genotype
 *         <p/>
 *         The interface for storing genotype calls.
 */
public interface GenotypeWriter {

    /**
     * Add a genotype, given a genotype locus
     * @param call the locus to add
     */
    public void addGenotypeCall(Genotype call);

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position the position to add the no call at
     */
    public void addNoCall(int position);

    /** finish writing, closing any open files. */
    public void close();

    /**
     * add a multi-sample call if we support it
     * @param genotypes the list of genotypes, that are backed by sample information
     */
    public void addMultiSampleCall(List<Genotype> genotypes, GenotypeMetaData metadata);

    /**
     * @return true if we support multisample, false otherwise
     */
    public boolean supportsMultiSample();

}
