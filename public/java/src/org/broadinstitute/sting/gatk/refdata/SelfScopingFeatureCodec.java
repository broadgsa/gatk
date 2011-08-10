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

package org.broadinstitute.sting.gatk.refdata;

import java.io.File;

/**
 * An interface marking that a given Tribble codec can look at the file and determine whether the
 * codec specifically parsing the contents of the file.
 */
public interface SelfScopingFeatureCodec {
    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     *
     * The GATK assumes that there's never a situation where two SelfScopingFeaetureCodecs
     * return true for the same file.  If this occurs the GATK splits out an error.
     *
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param potentialInput the file to test for parsiability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(final File potentialInput);
}
