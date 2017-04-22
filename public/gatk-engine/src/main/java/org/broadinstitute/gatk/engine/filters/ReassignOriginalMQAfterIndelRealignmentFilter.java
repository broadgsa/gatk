/*
* Copyright 2012-2016 Broad Institute, Inc.
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
import org.broadinstitute.gatk.utils.commandline.Argument;

/**
 * Revert the MQ of reads that were modified by IndelRealigner
 *
 * <p>IndelRealigner systematically adds +10 to the MQ of the reads it realigns. In some cases, that brings the resulting MQ to a value higher than MQ 60, which is the normal cap for MQ values. Since many downstream tools assume that MQ is <= 60, this may potentially cause problems.</p>
 *
 * <p>This read filter makes it possible to revert the MQ values of all the reads touched by IndelRealigner. It works by subtracting 10 from the MQ of all reads that have an "OC" tag, which stands for "original CIGAR" and is added by IndelRealigner to any read that it realigns.</p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Enable the filter</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         -rf ReassignOriginalMQAfterIndelRealignmentFilter
 * </pre>
 *
 *
 *
 * <h3>Caveat</h3>
 *
 * <p>There is currently no way to tell programmatically that a file has already been processed with this filter, so you should check the header manually before running this tool. Running it multiple times on the same BAM file would levy an unjustified penalty on realigned reads.</p>
 *
 * @author ami
 * @since 1/30/15.
 */
public class ReassignOriginalMQAfterIndelRealignmentFilter extends ReadFilter {

    public boolean filterOut(SAMRecord rec) {
        final String ORIGINAL_CIGAR_TAG = "OC";

        if (rec.getAttribute(ORIGINAL_CIGAR_TAG) != null)
            rec.setMappingQuality(rec.getMappingQuality() - 10);
        return false;
    }
}
