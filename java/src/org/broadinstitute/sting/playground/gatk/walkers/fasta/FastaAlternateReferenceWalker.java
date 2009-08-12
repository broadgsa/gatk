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
public class FastaAlternateReferenceWalker extends RefWalker<Pair<GenomeLoc, String>, Pair<GenomeLoc, String>> {

    @Argument(fullName="maskSNPs", shortName="mask", doc="print 'N' at SNP sites instead of the alternate allele", required=false)
    private Boolean MASK_SNPS = false;

	private StringBuffer sb = new StringBuffer();
    int deletionBasesRemaining = 0;

    public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        if ( deletionBasesRemaining > 0 ) {
            deletionBasesRemaining--;
            return new Pair<GenomeLoc, String>(context.getLocation(), "");
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
                return new Pair<GenomeLoc, String>(context.getLocation(), String.valueOf(ref.getBase()));
            } else if ( variant.isInsertion() ) {
                return new Pair<GenomeLoc, String>(context.getLocation(), String.valueOf(ref.getBase()).concat(variant.getAltBasesFWD()));
            } else if ( variant.isSNP() ) {
                if ( MASK_SNPS )
                    return new Pair<GenomeLoc, String>(context.getLocation(), "N");
                else
                    return new Pair<GenomeLoc, String>(context.getLocation(), variant.getAltBasesFWD());
            }
        }

        // if we got here then we're just ref
        return new Pair<GenomeLoc, String>(context.getLocation(), String.valueOf(ref.getBase()));
	}

    public Pair<GenomeLoc, String> reduceInit() {
        return new Pair<GenomeLoc, String>(null, "");
    }

	public Pair<GenomeLoc, String> reduce(Pair<GenomeLoc, String> value, Pair<GenomeLoc, String> sum) {
        // if there is no interval to the left, then this is the first one
        if ( sum.first == null ) {
            sum.first = value.first;
            sum.second = value.second;
        }
        // if the intervals don't overlap, print out the leftmost one and start a new one
        // (end of contig or new interval)
        else if ( value.first.getStart() != sum.first.getStop() + 1 ) {
            printFasta(sum.first, sum.second);
            sum.first = value.first;
            sum.second = value.second;
        }
        // otherwise, merge them
        else {
            sum.first = GenomeLocParser.setStop(sum.first,value.first.getStop());
            sum.second = sum.second.concat(value.second);
        }
		return sum;
	}

    public void onTraversalDone(Pair<GenomeLoc, String> sum) {
        if (sum.second != null)
            printFasta(sum.first, sum.second);
    }

    private void printFasta(GenomeLoc loc, String s) {
        out.println(">" + loc);
        int lines = s.length() / 60;
        int currentStart = 0;
        for (int i=0; i < lines; i++) {
            out.println(s.substring(currentStart, currentStart+60));
            currentStart += 60;
        }
        out.println(s.substring(currentStart));
    }
}