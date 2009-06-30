package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Pair;

// create a fasta sequence file from a reference and intervals

@WalkerName("FastaReferenceMaker")
public class FastaReferenceWalker extends RefWalker<Pair<GenomeLoc, Character>, Pair<GenomeLoc, String>> {

	public Pair<GenomeLoc, Character> map(RefMetaDataTracker rodData, char ref, LocusContext context) {
        return new Pair<GenomeLoc, Character>(context.getLocation(), ref);
	}

    public Pair<GenomeLoc, String> reduceInit() {
        return new Pair<GenomeLoc, String>(null, "");
    }

	public Pair<GenomeLoc, String> reduce(Pair<GenomeLoc, Character> value, Pair<GenomeLoc, String> sum) {
        // if there is no interval to the left, then this is the first one
        if ( sum.first == null ) {
            sum.first = value.first;
            sum.second = value.second.toString();
        }
        // if the intervals don't overlap, print out the leftmost one and start a new one
        else if ( value.first.getStart() != sum.first.getStop() + 1 ) {
            printFasta(sum.first, sum.second);
            sum.first = value.first;
            sum.second = value.second.toString();
        }
        // otherwise, merge them
        else {
            sum.first = GenomeLocParser.setStop(sum.first,value.first.getStop());
            sum.second = new String(sum.second + value.second);
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