package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

/**
 * Generates an alternative reference sequence over the specified interval.  Given variant ROD tracks,
 * it replaces the reference bases at variation sites with the bases supplied by the ROD(s).  Additionally,
 * allows for a "snpmask" ROD to set overlapping bases to 'N'.
 */
@WalkerName("FastaAlternateReferenceMaker")
@Requires(value={DataSource.REFERENCE})
public class FastaAlternateReferenceWalker extends FastaReferenceWalker {

    @Argument(fullName="outputIndelPositions", shortName="indelPositions", doc="output the positions of the indels in the new reference", required=false)
    String indelsFile = null;

    private int deletionBasesRemaining = 0;
    private long basesSeen = 0;
    private PrintWriter indelsWriter = null;

    public void initialize() {
        super.initialize();
        if (indelsFile != null) {
            try {
                indelsWriter = new PrintWriter(indelsFile);
            } catch (IOException e) {
                throw new RuntimeException("Unable to open indel positions output file: " + indelsFile);
            }
        }
    }

    public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        String refBase = String.valueOf(ref.getBase());

        if (deletionBasesRemaining > 0) {
            deletionBasesRemaining--;
            return new Pair<GenomeLoc, String>(context.getLocation(), "");
        }

        Iterator<GATKFeature> rods = rodData.getAllRods().iterator();
        while (rods.hasNext()) {
            GATKFeature rod = rods.next();
            if (!(rod.getUnderlyingObject() instanceof Variation))
                continue;
            // if we have multiple variants at a locus, just take the first damn one we see for now
            Variation variant = (Variation) rod.getUnderlyingObject();
            if (!rod.getName().startsWith("snpmask") && variant.isDeletion()) {
                deletionBasesRemaining = variant.getAlleleList().get(0).length();
                basesSeen++;
                if (indelsWriter != null)
                    indelsWriter.println(fasta.getCurrentID() + ":" + basesSeen + "-" + (basesSeen + variant.getAlleleList().get(0).length()));
                // delete the next n bases, not this one
                return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
            } else if (!rod.getName().startsWith("snpmask") && variant.isInsertion()) {
                basesSeen++;
                if (indelsWriter != null)
                    indelsWriter.println(fasta.getCurrentID() + ":" + basesSeen + "-" + (basesSeen + variant.getAlleleList().get(0).length()));
                basesSeen += variant.getAlleleList().get(0).length();
                return new Pair<GenomeLoc, String>(context.getLocation(), refBase.concat(Utils.join("",variant.getAlleleList())));
            } else if (variant.isSNP()) {
                basesSeen++;
                return new Pair<GenomeLoc, String>(context.getLocation(), (rod.getName().startsWith("snpmask") ? "N" : String.valueOf(variant.getAlternativeBaseForSNP())));
            }
        }

        // if we got here then we're just ref
        basesSeen++;
        return new Pair<GenomeLoc, String>(context.getLocation(), refBase);
    }

    public void onTraversalDone(GenomeLoc sum) {
        super.onTraversalDone(sum);
        if (indelsWriter != null)
            indelsWriter.close();
    }

}