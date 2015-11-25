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

/**
 * Configuration for the BWA/C aligner.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAConfiguration {
    /**
     * The maximum edit distance used by BWA.
     */
    public Float maximumEditDistance = null;

    /**
     * How many gap opens are acceptable within this alignment?
     */
    public Integer maximumGapOpens = null;

    /**
     * How many gap extensions are acceptable within this alignment?
     */
    public Integer maximumGapExtensions = null;

    /**
     * Do we disallow indels within a certain range from the start / end?
     */
    public Integer disallowIndelWithinRange = null;

    /**
     * What is the scoring penalty for a mismatch?
     */
    public Integer mismatchPenalty = null;

    /**
     * What is the scoring penalty for a gap open?
     */
    public Integer gapOpenPenalty = null;

    /**
     * What is the scoring penalty for a gap extension?
     */
    public Integer gapExtensionPenalty = null;

    /**
     * Enter bwa's 'non-stop' mode (equivalent to bwa aln -N parameter).
     */
    public Boolean nonStopMode = false;

    /**
     * Set the max queue size that bwa will use when searching for matches (equivalent to bwa aln -m parameter).
     */
    public Integer maxEntriesInQueue = null;
}
