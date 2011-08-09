/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;

/**
 * A read filter (transformer) that sets all reads mapping quality to a given value.
 *
 *  <p>
 *     If a BAM file contains erroneous or missing mapping qualities, this 'filter' will set
 *     all your mapping qualities to a given value. Default being 60.
 *  </p>
 *
 *
 * <h2>Input</h2>
 *  <p>
 *	    BAM file(s)
 *  </p>
 *
 *
 * <h2>Output</h2>
 *  <p>
 *      BAM file(s) with all reads mapping qualities reassigned
 *  </p>
 *
 * <h2>Examples</h2>
 *  <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -rf ReassignMappingQuality
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

