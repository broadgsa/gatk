/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;


/**
 * Left-align indels within reads in a bam file
 *
 * <p>This tool left-aligns any indels within read cigars in order to standardize representation when there are multiple valid
 * representations possible (i.e. where the same indel can be placed at multiple positions and still represent the same haplotype).
 * The standard convention is to place an indel at the left-most position possible, but this is not always followed, so
 * this tool can be used to correct the representation of indels.</p>
 *
 * <h3>Note</h3>
 * <p>This is only really needed when calling variants with legacy locus-based tools such as UnifiedGenotyper. With more
 * sophisticated tools (like HaplotypeCaller) that involve reconstructing haplotypes (eg through haplotype assembly), the problem
 * of multiple valid representations is handled internally and does not need to be corrected explicitly.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A bam file with mapped reads.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A bam file in which indels have been left-aligned where appropriate.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R reference.fasta \
 *   -T LeftAlignIndels \
 *   -I reads.bam \
 *   -o output_with_leftaligned_indels.bam
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class} )
public class LeftAlignIndels extends ReadWalker<Integer, Integer> {

    @Output(required=false, doc="Output bam")
    protected GATKSAMFileWriter writer = null;

    public void initialize() {}

    private void emit(final SAMRecord read) {
        if ( writer != null )
            writer.addAlignment(read);
    }

    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        // we can not deal with screwy records
        if ( read.getReadUnmappedFlag() || read.getCigar().numCigarElements() == 0 ) {
            emit(read);
            return 0;
        }

        // move existing indels (for 1 indel reads only) to leftmost position within identical sequence
        int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
        if ( numBlocks == 2 ) {
            Cigar newCigar = AlignmentUtils.leftAlignIndel(IndelRealigner.unclipCigar(read.getCigar()), ref.getBases(), read.getReadBases(), 0, 0, true);
            newCigar = IndelRealigner.reclipCigar(newCigar, read);
            read.setCigar(newCigar);
        }

        emit(read);
        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {}
}