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

/**
 * Filter out reads whose mate maps to a different contig
 *
 * <p>This filter is intended to ensure that only reads that are likely to be mapped in the right place, and therefore
 * to be informative, will be used in analysis. If mates in a pair are mapping to different contigs, it is likely that
 * at least one of them is in the wrong place. One exception is you are using a draft genome assembly in which the
 * chromosomes are fragmented into many contigs; then you may legitimately have reads that are correctly mapped but are
 * on different contigs than their mate. This read filter can be disabled from the command line using the -drf argument.
 * </p>
 *
 * <h4>Enable the bad mate filter</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         -rf BadMate
 * </pre>
 *
 * <h4>Disable the bad mate filter</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         <b>-drf</b> BadMate
 * </pre>
 *
 * @author ebanks
 * @version 0.1
 */

public class BadMateFilter extends DisableableReadFilter {

    public boolean filterOut(final SAMRecord rec) {
        return hasBadMate(rec);
    }

    public static boolean hasBadMate(final SAMRecord rec) {
        return (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag() && !rec.getReferenceIndex().equals(rec.getMateReferenceIndex()));
    }

}