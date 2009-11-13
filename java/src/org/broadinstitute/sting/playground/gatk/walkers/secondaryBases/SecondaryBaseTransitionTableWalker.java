package org.broadinstitute.sting.playground.gatk.walkers.secondaryBases;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.IntervalRod;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodGenotypeChipAsGFF;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import javax.xml.transform.Result;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;


/**
 * Created by IntelliJ IDEA.
 * User: michael
 * Date: Nov 2, 2009
 * Time: 8:45:32 PM
 * To change this template use File | Settings | File Templates.
 */
@Reference(window=@Window(start=-1,stop=1))
public class SecondaryBaseTransitionTableWalker extends LocusWalker<Integer, Integer> {

    HashMap<String,Long> counts = new HashMap<String,Long>();
    public IndexedFastaSequenceFile refSeq;

    /*public void initialize() {
        File refFile = this.getToolkit().getArguments().referenceFile;

        try {
            refSeq = new IndexedFastaSequenceFile(refFile);
        } catch (IOException e) {
            refSeq = null;
        }
    }*/

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String refBase = Character.toString(ref.getBase());
        ReadBackedPileup pileup = new ReadBackedPileup(ref.getBase(),context);
        String primaryBases = pileup.getBasesWithStrand();
        String secondaryBases = pileup.getSecondaryBasePileup();
        String contextBases = new String(ref.getBases());

        /*if (refSeq != null) {
            long startPos = context.getPosition() - 1;
            long stopPos = context.getPosition() + 1;

            ReferenceSequence prevRefSequence = refSeq.getSubsequenceAt(context.getContig(), startPos, context.getPosition() - 1);
            ReferenceSequence nextRefSequence = refSeq.getSubsequenceAt(context.getContig(), context.getPosition() + 1, stopPos);

            String prev = new String(prevRefSequence.getBases());
            String next = new String(nextRefSequence.getBases());
            String total = new String(refSeq.getSubsequenceAt(context.getContig(),startPos,stopPos).getBases());
            out.println(total + "  " + prev + " " + refBase + " " + next);
        }*/

        String precedingBase = contextBases.substring(0,1);
        String nextBase = contextBases.substring(2);
        /*out.println(contextBases + "  " + precedingBase + " " + refBase + " " + nextBase);*/
        /*out.println(" ");*/

        boolean rods = false;
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) {
            if (datum != null && !(datum instanceof IntervalRod)) {
                rods = true;}
        }

        if (!rods && precedingBase != null && secondaryBases != null) {
            for (int i = 0; i < primaryBases.length(); i ++) {
                if (secondaryBases.charAt(i) != 'N' && secondaryBases.charAt(i) != '.' && Character.toUpperCase(primaryBases.charAt(i)) != 'N') {
                    String quenchingBase;
                    if (primaryBases.charAt(i) == Character.toUpperCase(primaryBases.charAt(i))) {
                        quenchingBase = precedingBase;
                    }
                    else {
                        quenchingBase = nextBase;
                    }
                    String reference = Character.toString(Character.toUpperCase(refBase.charAt(0)));
                    String primary = Character.toString(Character.toUpperCase(primaryBases.charAt(i)));
                    String quencher = Character.toString(Character.toUpperCase(quenchingBase.charAt(0)));
                    String secondary = Character.toString(Character.toUpperCase(secondaryBases.charAt(i)));
                    String key = reference + primary + quencher + secondary;
                    if (counts.containsKey(key)) {
                        counts.put(key, counts.get(key) + Long.valueOf(1));
                    }
                    else {
                        counts.put(key, Long.valueOf(1));
                    }
                }
            }
        }
        return 1;
    }

    public Integer reduceInit() {return 0;}

    public Integer reduce(Integer value, Integer sum) {return sum + value;}

    public void onTraversalDone(Integer result) {
        out.println("ReferenceBase \tPrimaryBase \tPreviousBase \tSecondaryBase \tCount");
        for (String key : counts.keySet()) {
            out.println(key.charAt(0)+"\t"+key.charAt(1)+"\t"+key.charAt(2)+"\t"+key.charAt(3)+"\t"+counts.get(key));
        }
    }
}
