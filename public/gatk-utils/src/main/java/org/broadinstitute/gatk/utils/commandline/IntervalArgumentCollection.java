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

package org.broadinstitute.gatk.utils.commandline;

import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.interval.IntervalSetRule;

import java.util.List;

public class IntervalArgumentCollection {
    /**
     * Use this option to perform the analysis over only part of the genome. This argument can be specified multiple times.
     * You can use samtools-style intervals either explicitly on the command line (e.g. -L chr1 or -L chr1:100-200) or
     * by loading in a file containing a list of intervals (e.g. -L myFile.intervals).
     *
     * Additionally, you can also specify a ROD file (such as a VCF file) in order to perform the analysis at specific
     * positions based on the records present in the file (e.g. -L file.vcf).
     *
     * Finally, you can also use this to perform the analysis on the reads that are completely unmapped in the BAM file
     * (i.e. those without a reference contig) by specifying -L unmapped.
     */
    @Input(fullName = "intervals", shortName = "L", doc = "One or more genomic intervals over which to operate", required = false)
    public List<IntervalBinding<Feature>> intervals = null;

    /**
     * Use this option to exclude certain parts of the genome from the analysis (like -L, but the opposite).
     * This argument can be specified multiple times. You can use samtools-style intervals either explicitly on the
     * command line (e.g. -XL chr1 or -XL chr1:100-200) or by loading in a file containing a list of intervals
     * (e.g. -XL myFile.intervals).
     *
     * Additionally, you can also specify a ROD file (such as a VCF file) in order to exclude specific
     * positions from the analysis based on the records present in the file (e.g. -XL file.vcf).
     * */
    @Input(fullName = "excludeIntervals", shortName = "XL", doc = "One or more genomic intervals to exclude from processing", required = false)
    public List<IntervalBinding<Feature>> excludeIntervals = null;

    /**
     * By default, the program will take the UNION of all intervals specified using -L and/or -XL. However, you can
     * change this setting for -L, for example if you want to take the INTERSECTION of the sets instead. E.g. to perform the
     * analysis on positions for which there is a record in a VCF, but restrict this to just those on chromosome 20,
     * you would do -L chr20 -L file.vcf -isr INTERSECTION. However, it is not possible to modify the merging approach
     * for intervals passed using -XL (they will always be merged using UNION).
     *
     * Note that if you specify both -L and -XL, the -XL interval set will be subtracted from the -L interval set.
     */
    @Argument(fullName = "interval_set_rule", shortName = "isr", doc = "Set merging approach to use for combining interval inputs", required = false)
    public IntervalSetRule intervalSetRule = IntervalSetRule.UNION;

    /**
     * By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not
     * actually overlap) into a single continuous interval. However you can change this behavior if you want them to be
     * treated as separate intervals instead.
     */
    @Argument(fullName = "interval_merging", shortName = "im", doc = "Interval merging rule for abutting intervals", required = false)
    public IntervalMergingRule intervalMerging = IntervalMergingRule.ALL;

    /**
     * Use this to add padding to the intervals specified using -L and/or -XL. For example, '-L chr1:100' with a
     * padding value of 20 would turn into '-L chr1:80-120'. This is typically used to add padding around exons when
     * analyzing exomes. The general Broad exome calling pipeline uses 100 bp padding by default.
     */
    @Argument(fullName = "interval_padding", shortName = "ip", doc = "Amount of padding (in bp) to add to each interval", required = false, minValue = 0)
    public int intervalPadding = 0;
}
