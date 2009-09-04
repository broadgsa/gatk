package org.broadinstitute.sting.playground.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Pair;

// create a fasta sequence file from a reference and intervals

@WalkerName("FastaReferenceMaker")
public class FastaReferenceWalker extends RefWalker<Pair<GenomeLoc, String>, GenomeLoc> {

    protected FastaSequence fasta;

    public void initialize() {
        fasta = new FastaSequence(out);
    }

	public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        return new Pair<GenomeLoc, String>(context.getLocation(), String.valueOf(ref.getBase()));
	}

    public GenomeLoc reduceInit() {
        return null;
    }

	public GenomeLoc reduce(Pair<GenomeLoc, String> value, GenomeLoc sum) {
        // if there is no interval to the left, then this is the first one
        if ( sum == null ) {
            sum = value.first;
            fasta.append(value.second.toString());
        }
        // if the intervals don't overlap, print out the leftmost one and start a new one
        // (end of contig or new interval)
        else if ( value.first.getStart() != sum.getStop() + 1 ) {
            fasta.flush();
            sum = value.first;
            fasta.append(value.second.toString());
        }
        // otherwise, merge them
        else {
            sum = GenomeLocParser.setStop(sum, value.first.getStop());
            fasta.append(value.second.toString());
        }
		return sum;
	}

    public void onTraversalDone(GenomeLoc sum) {
        fasta.flush();
    }
}