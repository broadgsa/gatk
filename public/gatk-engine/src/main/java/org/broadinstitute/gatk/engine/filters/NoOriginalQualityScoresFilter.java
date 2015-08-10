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

package org.broadinstitute.gatk.engine.filters;

import htsjdk.samtools.SAMRecord;


/**
 * Filter out reads that do not have an original quality quality score (OQ) tag
 *
 * <p>The OQ tag can be added during the base recalibration process in order to preserve original information.</p>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         -rf NoOriginalQualityScores
 * </pre>
 *
 * @author rpoplin
 * @since Nov 19, 2009
 */
public class NoOriginalQualityScoresFilter extends ReadFilter {
    public boolean filterOut( final SAMRecord read ) {
        return (read.getAttribute("OQ") == null);
    }
}
