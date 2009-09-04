package org.broadinstitute.sting.playground.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.*;
import java.util.Iterator;


@WalkerName("FastaAlternateReferenceMaker")
@Requires(value={DataSource.REFERENCE})
public class FastaAlternateReferenceWalker extends FastaReferenceWalker {

    @Argument(fullName="outputSequenomFormat", shortName="sequenom", doc="output results in sequenom format (overrides 'maskSNPs' argument)", required=false)
    private Boolean SEQUENOM = false;
    @Argument(fullName="outputIndelPositions", shortName="indelPositions", doc="output the positions of the indels in the new reference", required=false)
    String indelsFile = null;

    private int deletionBasesRemaining = 0;
    private long basesSeen = 0;
    private PrintWriter indelsWriter = null;

    public void initialize() {
        super.initialize();
        if ( indelsFile != null ) {
            try {
                indelsWriter = new PrintWriter(indelsFile);
            } catch (IOException e) {
                throw new RuntimeException("Unable to open indel positions output file: " + indelsFile);
            }
        }
    }

    public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        String refBase = String.valueOf(ref.getBase());

        if ( deletionBasesRemaining > 0 ) {
            deletionBasesRemaining--;
            return new Pair<GenomeLoc, String>(context.getLocation(), (SEQUENOM ? (deletionBasesRemaining == 0 ? refBase.concat("/-]") : refBase) : ""));
        }

        Iterator<ReferenceOrderedDatum> rods = rodData.getAllRods().iterator();
        while ( rods.hasNext() ) {
            ReferenceOrderedDatum rod = rods.next();
            if ( !(rod instanceof AllelicVariant) )
                continue;

            // if we have multiple variants at a locus, just take the first damn one we see for now
            AllelicVariant variant = (AllelicVariant)rod;
            if ( !rod.getName().startsWith("snpmask") && variant.isDeletion() ) {
                deletionBasesRemaining = variant.length();
                basesSeen++;
                if ( indelsWriter != null )
                    indelsWriter.println(fasta.getCurrentID() + ":" + basesSeen + "-" + (basesSeen + variant.length()));
                // delete the next n bases, not this one
                return new Pair<GenomeLoc, String>(context.getLocation(), (SEQUENOM ? refBase.concat("[") : refBase));
            } else if ( !rod.getName().startsWith("snpmask") && variant.isInsertion() ) {
                basesSeen++;
                if ( indelsWriter != null )
                    indelsWriter.println(fasta.getCurrentID() + ":" + basesSeen + "-" + (basesSeen + variant.length()));
                basesSeen += variant.length();
                return new Pair<GenomeLoc, String>(context.getLocation(), (SEQUENOM ? refBase.concat("[-/"+variant.getAltBasesFWD()+"]") : refBase.concat(variant.getAltBasesFWD())));
            } else if ( variant.isSNP() ) {
                basesSeen++;
                return new Pair<GenomeLoc, String>(context.getLocation(), (rod.getName().startsWith("snpmask") ? "N" : variant.getAltBasesFWD()));
            }
        }

        // if we got here then we're just ref
        basesSeen++;
        return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
	}

    public void onTraversalDone(GenomeLoc sum) {
        super.onTraversalDone(sum);
        if ( indelsWriter != null )
            indelsWriter.close();
    }

}