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

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.commandline.Argument;

import java.util.Set;

/**
 * Only use reads belonging to a specific sample
 *
 * <p>This filter is useful for isolating data from one particular sample in a multisample file.</p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Use only reads from the sample named NA12878</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         -rf Sample \
 *         -goodSM NA12878
 * </pre>
 *
 */
public class SampleFilter extends ReadFilter {
    @Argument(fullName = "sample_to_keep", shortName = "goodSM", doc="The name of the sample(s) to keep, filtering out all others", required=true)
    private Set SAMPLES_TO_KEEP = null;

    public boolean filterOut( final SAMRecord read ) {
        final SAMReadGroupRecord readGroup = read.getReadGroup();
        return !( readGroup != null && SAMPLES_TO_KEEP.contains(readGroup.getSample()) );
    }
}
