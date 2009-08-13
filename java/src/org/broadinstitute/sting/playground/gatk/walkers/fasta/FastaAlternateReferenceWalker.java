package org.broadinstitute.sting.playground.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.Iterator;

// create a fasta sequence file from a reference, intervals, and rod(s) of variants
// if there are multiple variants at a site, we take the first one seen

@WalkerName("FastaAlternateReferenceMaker")
@Requires(value={DataSource.REFERENCE})
public class FastaAlternateReferenceWalker extends FastaReferenceWalker {

    @Argument(fullName="maskSNPs", shortName="mask", doc="print 'N' at SNP sites instead of the alternate allele", required=false)
    private Boolean MASK_SNPS = false;
    @Argument(fullName="outputSequenomFormat", shortName="sequenom", doc="output results in sequenom format (overrides 'maskSNPs' argument)", required=false)
    private Boolean SEQUENOM = false;

    int deletionBasesRemaining = 0;

    public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        String refBase = String.valueOf(ref.getBase());

        if ( deletionBasesRemaining > 0 ) {
            deletionBasesRemaining--;
            return new Pair<GenomeLoc, String>(context.getLocation(), (SEQUENOM ? (deletionBasesRemaining == 0 ? refBase.concat("]") : refBase) : ""));
        }

        Iterator<ReferenceOrderedDatum> rods = rodData.getAllRods().iterator();
        while ( rods.hasNext() ) {
            ReferenceOrderedDatum rod = rods.next();
            if ( !(rod instanceof AllelicVariant) )
                continue;

            // if we have multiple variants at a locus, just take the first damn one we see for now
            AllelicVariant variant = (AllelicVariant)rod;
            if ( variant.isDeletion() ) {
                deletionBasesRemaining = variant.length();
                // delete the next n bases, not this one
                return new Pair<GenomeLoc, String>(context.getLocation(), (SEQUENOM ? refBase.concat("[-/") : refBase));
            } else if ( variant.isInsertion() ) {
                return new Pair<GenomeLoc, String>(context.getLocation(), (SEQUENOM ? refBase.concat("[+/"+variant.getAltBasesFWD()+"]") : refBase.concat(variant.getAltBasesFWD())));
            } else if ( variant.isSNP() ) {
                return new Pair<GenomeLoc, String>(context.getLocation(), (SEQUENOM || MASK_SNPS ? "N" : variant.getAltBasesFWD()));
            }
        }

        // if we got here then we're just ref
        return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
	}
}