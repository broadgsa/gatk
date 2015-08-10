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
 * Set the mapping quality of all reads to a given value.
 *
 *  <p>
 *     If a BAM file contains erroneous or missing mapping qualities (MAPQ), this read transformer will set all your
 *     mapping qualities to a given value (see arguments list for default value).
 *  </p>
 *
 * <h3>See also</h3>
 *
 * <p>ReassignOneMappingQualityFilter: reassigns a single MAPQ value, as opposed to all those found in the BAM file.</p>
 *
 * <h3>Caveats</h3>
 *
 * <p>Note that due to the order of operations involved in applying filters, it is possible that other read filters
 * (determined either at command-line or internally by the tool you are using) will be applied to your data before
 * this read transformation can be applied. If one of those other filters acts on the read mapping quality (MAPQ),
 * then you may not obtain the expected results. Unfortunately it is currently not possible to change the order of
 * operations from command line. To avoid the problem, we recommend applying this filter separately from any other
 * analysis, using PrintReads.</p>
 *
 *
 * <h3>Input</h3>
 *  <p>
 *	    BAM file(s)
 *  </p>
 *
 * <h3>Output</h3>
 *  <p>
 *      BAM file(s) with the mapping qualities of all reads reassigned to the specified value
 *  </p>
 *
 * <h3>Usage example</h3>
 *  <pre>
 *  java -jar GenomeAnalysisTK.jar \
 *      -T PrintReads \
 *      -R reference.fasta \
 *      -I input.bam \
 *      -o output.file \
 *      -rf ReassignMappingQuality \
 *      -DMQ 35
 *  </pre>
 *
 * @author carneiro
 * @since 8/8/11
 */

public class ReassignMappingQualityFilter extends ReadFilter {

    @Argument(fullName = "default_mapping_quality", shortName = "DMQ", doc = "Default read mapping quality to assign to all reads", required = false)
    public int defaultMappingQuality = 60;

    public boolean filterOut(SAMRecord rec) {
        rec.setMappingQuality(defaultMappingQuality);
        return false;
    }
}

