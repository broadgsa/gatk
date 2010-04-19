package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;


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