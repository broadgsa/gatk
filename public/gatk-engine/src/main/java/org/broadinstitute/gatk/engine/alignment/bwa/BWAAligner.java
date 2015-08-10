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

package org.broadinstitute.gatk.engine.alignment.bwa;

import org.broadinstitute.gatk.engine.alignment.Aligner;

/**
 * Align reads using BWA.
 *
 * @author mhanna
 * @version 0.1
 */
public abstract class BWAAligner implements Aligner {
    /**
     * The supporting files used by BWA.
     */
    protected BWTFiles bwtFiles;

    /**
     * The current configuration for the BWA aligner.
     */
    protected BWAConfiguration configuration;

    /**
     * Create a new BWAAligner.  Purpose of this call is to ensure that all BWA constructors accept the correct
     * parameters.
     * @param bwtFiles The many files representing BWTs persisted to disk.
     * @param configuration Configuration parameters for the alignment.
     */
    public BWAAligner(BWTFiles bwtFiles, BWAConfiguration configuration) {
        this.bwtFiles = bwtFiles;
        this.configuration = configuration;
    }

    /**
     * Update the configuration passed to the BWA aligner.
     * @param configuration New configuration to set.
     */    
    public abstract void updateConfiguration(BWAConfiguration configuration);
}
