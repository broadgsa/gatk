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
 * Filter out reads with low mapping qualities
 *
 * <p>This filter is intended to ensure that only reads that are likely
 * to be mapped in the right place, and therefore to be informative, will be used in analysis.</p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Set the mapping quality filter to filter out reads that have MAPQ < 15</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T HaplotypeCaller \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.vcf \
 *         -rf MappingQuality \
 *         -mmq 15
 * </pre>
 *
 * @author ebanks
 * @version 0.1
 */

public class MappingQualityFilter extends ReadFilter {

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 10;

    public boolean filterOut(SAMRecord rec) {
        return (rec.getMappingQuality() < MIN_MAPPING_QUALTY_SCORE);
    }
}