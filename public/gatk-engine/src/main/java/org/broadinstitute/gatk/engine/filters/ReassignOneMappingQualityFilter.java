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
 * Set the mapping quality of reads with a given value to another given value.
 *
 *  <p>
 *     This read transformer will change a certain read mapping quality to a different value without affecting reads that
 *     have other mapping qualities. This is intended primarily for users of RNA-Seq data handling programs such
 *     as TopHat, which use MAPQ = 255 to designate uniquely aligned reads. According to convention, 255 normally
 *     designates "unknown" quality, and most GATK tools automatically ignore such reads. By reassigning a different
 *     mapping quality to those specific reads, users of TopHat and other tools can circumvent this problem without
 *     affecting the rest of their dataset.
 *  </p>
 *
 *  <p>
 *     This differs from the ReassignMappingQuality filter by its selectivity -- only one mapping quality is targeted.
 *     ReassignMappingQuality will change ALL mapping qualities to a single one, and is typically used for datasets
 *     that have no assigned mapping qualities.
 *  </p>
 *
 * <h3>Input</h3>
 *  <p>
 *	    BAM file(s)
 *  </p>
 *
 *
 * <h3>Output</h3>
 *  <p>
 *      BAM file(s) with one read mapping quality selectively reassigned as desired
 *  </p>
 *
 * <h3>Usage example</h3>
 *  <pre>
 *    java -jar GenomeAnalysisTK.jar \
 *      -T PrintReads \
 *      -R reference.fasta \
 *      -I input.bam \
 *      -o output.file \
 *      -rf ReassignOneMappingQuality \
 *      -RMQF 255 \
 *      -RMQT 60
 *  </pre>
 *
 * @author vdauwera
 * @since 2/19/13
 */

public class ReassignOneMappingQualityFilter extends ReadFilter {

    @Argument(fullName = "reassign_mapping_quality_from", shortName = "RMQF", doc = "Original mapping quality", required = false)
    public int reassignMappingQualityFrom = 255;

    @Argument(fullName = "reassign_mapping_quality_to", shortName = "RMQT", doc = "Desired mapping quality", required = false)
    public int reassignMappingQualityTo = 60;

    public boolean filterOut(SAMRecord rec) {
        if (rec.getMappingQuality() == reassignMappingQualityFrom)
        rec.setMappingQuality(reassignMappingQualityTo);
        return false;
    }
}

