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

package org.broadinstitute.sting.commandline;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;

import java.util.List;

public class IntervalArgumentCollection {
    /**
     * Using this option one can instruct the GATK engine to traverse over only part of the genome.  This argument can be specified multiple times.
     * One may use samtools-style intervals either explicitly (e.g. -L chr1 or -L chr1:100-200) or listed in a file (e.g. -L myFile.intervals).
     * Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf).
     * To specify the completely unmapped reads in the BAM file (i.e. those without a reference contig) use -L unmapped.
     */
    @Input(fullName = "intervals", shortName = "L", doc = "One or more genomic intervals over which to operate. Can be explicitly specified on the command line or in a file (including a rod file)", required = false)
    public List<IntervalBinding<Feature>> intervals = null;

    /**
     * Using this option one can instruct the GATK engine NOT to traverse over certain parts of the genome.  This argument can be specified multiple times.
     * One may use samtools-style intervals either explicitly (e.g. -XL chr1 or -XL chr1:100-200) or listed in a file (e.g. -XL myFile.intervals).
     * Additionally, one may specify a rod file to skip over the positions for which there is a record in the file (e.g. -XL file.vcf).
     */
    @Input(fullName = "excludeIntervals", shortName = "XL", doc = "One or more genomic intervals to exclude from processing. Can be explicitly specified on the command line or in a file (including a rod file)", required = false)
    public List<IntervalBinding<Feature>> excludeIntervals = null;

    /**
     * How should the intervals specified by multiple -L or -XL arguments be combined?  Using this argument one can, for example, traverse over all of the positions
     * for which there is a record in a VCF but just in chromosome 20 (-L chr20 -L file.vcf -isr INTERSECTION).
     */
    @Argument(fullName = "interval_set_rule", shortName = "isr", doc = "Indicates the set merging approach the interval parser should use to combine the various -L or -XL inputs", required = false)
    public IntervalSetRule intervalSetRule = IntervalSetRule.UNION;

    /**
     * Should abutting (but not overlapping) intervals be treated as separate intervals?
     */
    @Argument(fullName = "interval_merging", shortName = "im", doc = "Indicates the interval merging rule we should use for abutting intervals", required = false)
    public IntervalMergingRule intervalMerging = IntervalMergingRule.ALL;

    /**
     * For example, '-L chr1:100' with a padding value of 20 would turn into '-L chr1:80-120'.
     */
    @Argument(fullName = "interval_padding", shortName = "ip", doc = "Indicates how many basepairs of padding to include around each of the intervals specified with the -L/--intervals argument", required = false)
    public int intervalPadding = 0;
}
