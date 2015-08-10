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
 * Filter out reads based on length
 *
 * <p>This filter is useful for running on only reads that are longer (or shorter) than the given threshold sizes. </p>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         -rf ReadLength \
 *         -minRead 50 \
 *         -maxRead 101
 * </pre>
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadLengthFilter extends ReadFilter {
    @Argument(fullName = "maxReadLength", shortName = "maxRead", doc="Discard reads with length greater than the specified value", required=true)
    private int maxReadLength;

    @Argument(fullName = "minReadLength", shortName = "minRead", doc="Discard reads with length shorter than the specified value", required=true)
    private int minReadLength = 1;
    public boolean filterOut(SAMRecord read) {
        // check the length
        return read.getReadLength() > maxReadLength || read.getReadLength() < minReadLength;
    }

}
