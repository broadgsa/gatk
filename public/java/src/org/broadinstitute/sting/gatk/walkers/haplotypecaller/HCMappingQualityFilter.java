/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.ReadFilter;

/**
 * Filter out reads with low mapping qualities.
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