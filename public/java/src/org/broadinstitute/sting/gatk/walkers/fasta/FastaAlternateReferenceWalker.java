/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Collections;
import java.util.List;


/**
 * Generates an alternative reference sequence over the specified interval.
 *
 * <p>
 * Given variant tracks, it replaces the reference bases at variation sites with the bases supplied by the ROD(s).
 * Additionally, allows for one or more "snpmask" VCFs to set overlapping bases to 'N'.
 * Note that if there are multiple variants at a site, it takes the first one seen.
 * Reference bases for each interval will be output as a separate fasta sequence (named numerically in order).
 *
 * <h2>Input</h2>
 * <p>
 * The reference, requested intervals, and any number of variant rod files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A fasta file representing the requested intervals.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T FastaAlternateReferenceMaker \
 *   -o output.fasta \
 *   -L input.intervals \
 *   --variant input.vcf \
 *   [--snpmask mask.vcf]
 * </pre>
 *
 */
@WalkerName("FastaAlternateReferenceMaker")
@Reference(window=@Window(start=-1,stop=50))
@Requires(value={DataSource.REFERENCE})
public class FastaAlternateReferenceWalker extends FastaReferenceWalker {

    /**
     * Variants from these input files are used by this tool to construct an alternate reference.
     */
    @Input(fullName = "variant", shortName = "V", doc="variants to model", required=false)
    public List<RodBinding<VariantContext>> variants = Collections.emptyList();

    /**
     * Snps from this file are used as a mask when constructing the alternate reference.
     */
    @Input(fullName="snpmask", shortName = "snpmask", doc="SNP mask VCF file", required=false)
    public RodBinding<VariantContext> snpmask;

    private int deletionBasesRemaining = 0;

    public Pair<GenomeLoc, String> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (deletionBasesRemaining > 0) {
            deletionBasesRemaining--;
            return new Pair<GenomeLoc, String>(context.getLocation(), "");
        }

        String refBase = String.valueOf((char)ref.getBase());

        // Check to see if we have a called snp
        for ( VariantContext vc : tracker.getValues(variants) ) {
            if ( vc.isFiltered() )
                continue;

            if ( vc.isSimpleDeletion()) {
                deletionBasesRemaining = vc.getReference().length();
                // delete the next n bases, not this one
                return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
            } else if ( vc.isSimpleInsertion()) {
                return new Pair<GenomeLoc, String>(context.getLocation(), refBase.concat(vc.getAlternateAllele(0).toString()));
            } else if (vc.isSNP()) {
                return new Pair<GenomeLoc, String>(context.getLocation(), vc.getAlternateAllele(0).toString());
            }
        }

        // if we don't have a called site, and we have a mask at this site, mask it
        for ( VariantContext vc : tracker.getValues(snpmask) ) {
            if ( vc.isSNP()) {
                return new Pair<GenomeLoc, String>(context.getLocation(), "N");
            }
        }


        // if we got here then we're just ref
        return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
    }
}