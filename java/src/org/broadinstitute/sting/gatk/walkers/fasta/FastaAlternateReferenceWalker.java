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

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;


/**
 * Generates an alternative reference sequence over the specified interval.  Given variant ROD tracks,
 * it replaces the reference bases at variation sites with the bases supplied by the ROD(s).  Additionally,
 * allows for a "snpmask" ROD to set overlapping bases to 'N'.
 */
@WalkerName("FastaAlternateReferenceMaker")
@Reference(window=@Window(start=0,stop=50))
@Requires(value={DataSource.REFERENCE})
public class FastaAlternateReferenceWalker extends FastaReferenceWalker {

    private int deletionBasesRemaining = 0;

    public Pair<GenomeLoc, String> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (deletionBasesRemaining > 0) {
            deletionBasesRemaining--;
            return new Pair<GenomeLoc, String>(context.getLocation(), "");
        }

        String refBase = String.valueOf(ref.getBase());

        for ( VariantContext vc : tracker.getAllVariantContexts(ref) ) {
            // if we have multiple variants at a locus, just take the first one we see
            if (!vc.getName().startsWith("snpmask") && vc.isDeletion()) {
                deletionBasesRemaining = vc.getReference().length();
                // delete the next n bases, not this one
                return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
            } else if (!vc.getName().startsWith("snpmask") && vc.isInsertion()) {
                return new Pair<GenomeLoc, String>(context.getLocation(), refBase.concat(vc.getAlternateAllele(0).toString()));
            } else if (vc.isSNP()) {
                return new Pair<GenomeLoc, String>(context.getLocation(), (vc.getName().startsWith("snpmask") ? "N" : vc.getAlternateAllele(0).toString()));
            }
        }

        // if we got here then we're just ref
        return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
    }
}