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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.filters.ReadFilter;

/**
 * Filter out reads with low mapping qualities for HaplotypeCaller
 *
 * <p>This filter is applied by default for HaplotypeCaller and is designed to ensure that only reads that are likely
 * to be informative will be used in the reassembly process. It performs the same basic function as the regular
 * MappingQualityFilter, but it is used at specific points in the operation of HC where it is helpful
 * to be able to apply a different quality threshold from the general case.</p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Set the HC-specific mapping quality filter to filter out reads with MAPQ < 10</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T HaplotypeCaller \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.vcf \
 *         -mmq 10
 * </pre>
 *
 * <p>Note that the HCMappingQuality filter itself does not need to be specified in the command line because it is set
 * automatically for HaplotypeCaller.</p>
 *
 * @author mdepristo
 */
public class HCMappingQualityFilter extends ReadFilter {
    private final static Logger logger = Logger.getLogger(HCMappingQualityFilter.class);

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for analysis with the HaplotypeCaller", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 20;

    @Override
    public void initialize(GenomeAnalysisEngine engine) {
        if ( MIN_MAPPING_QUALTY_SCORE > 0 )
            logger.info("Filtering out reads with MAPQ < " + MIN_MAPPING_QUALTY_SCORE);
    }

    public boolean filterOut(SAMRecord rec) {
        return (rec.getMappingQuality() < MIN_MAPPING_QUALTY_SCORE);
    }
}