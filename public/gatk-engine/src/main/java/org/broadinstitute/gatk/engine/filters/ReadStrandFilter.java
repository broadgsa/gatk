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
import org.broadinstitute.gatk.utils.commandline.Argument;

/**
 * Filter out reads based on strand orientation
 *
 * <p>This filter is useful for isolating reads from only forward or reverse strands. By default, it filters out reads
 * from the negative (reverse) strand. This logic can be reversed by using the -filterPositive flag.</p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Set the read strand filter to filter out positive (forward) strand reads</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         -rf ReadStrand \
 *         -filterPositive
 * </pre>
 *
 * @author chartl
 * @version 0.1
 */
public class ReadStrandFilter extends ReadFilter {
    @Argument(fullName = "filterPositive", shortName = "fp", doc="Discard reads on the forward strand",required=false)
	boolean filterForward = false;
    
    public boolean filterOut(SAMRecord read) {
        // check the length
        return read.getReadNegativeStrandFlag() != filterForward;
    }

}
