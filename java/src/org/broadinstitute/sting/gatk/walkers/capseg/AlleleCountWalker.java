package org.broadinstitute.sting.gatk.walkers.capseg;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;
import java.util.Set;

/**
 * get the allele counts at loci that overlap both the bait list and DbSNP
 */
//@By(DataSource.REFERENCE_ORDERED_DATA)
public class AlleleCountWalker extends LocusWalker<Integer, Integer> {
    @Argument(doc = "our output file name", shortName = "fl")
    File output;

    String dbTrack = "DbSNP";
    String callTrack = "calls";
    Random generator = new Random();
    private PrintWriter out;

    public void initialize() {
        try {
            out = new PrintWriter(new FileWriter(output));
        } catch (IOException e) {
            throw new IllegalArgumentException("Uncaught exception!");
        }
        out.println("contig,position,A,B,aBase,bBase,dbSNPOrientation");
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        // check that the tracker and the context are not null, and we a variant, and some amount of bases
        if (tracker == null || context == null || tracker.getAllVariantContexts(ref).size() < 1 || context.getBasePileup().size() == 0)
            return 0;

        // get the het
        VariantContext het = tracker.getVariantContext(ref, callTrack, context.getLocation());
        VariantContext dbSNP = tracker.getVariantContext(ref, dbTrack, context.getLocation());
        if (het == null || het.getHetCount() == 0) return 0;


        byte[] bases = context.getBasePileup().getBases();
        int[] counts = new int[5];
        for (byte b : bases) {
            if (b == 'a' || b == 'A') counts[0]++;
            else if (b == 'c' || b == 'C') counts[1]++;
            else if (b == 'g' || b == 'G') counts[2]++;
            else if (b == 't' || b == 'T') counts[3]++;
            else counts[4]++;
        }
        int aIndex = 0;
        byte a = 'a';

        // find the maximumvalue
        for (int index = 0; index < counts.length; index++)
            if (counts[index] > counts[aIndex]) {
                aIndex = index;
                a = (byte) ((index < 1) ? 'a' : (index < 2) ? 'c' : (index < 3) ? 'g' : 't');
            }

        int bIndex = (aIndex != 0) ? 0 : 1;
        byte b = 'a';
        for (int index = 0; index < counts.length; index++)
            if (counts[index] > counts[bIndex] && index != aIndex) {
                bIndex = index;
                b = (byte) ((index < 1) ? 'a' : (index < 2) ? 'c' : (index < 3) ? 'g' : 't');
            }

        boolean usedDbSNP = false;
        // check if there is an ordering that we'd like to subscribe to
        if (dbSNP != null) {
            // get the alt from the DbSNP file
            Set<Allele> alts = dbSNP.getAlternateAlleles();

            // if true, swap the two
            if (alts.size() == 1 && alts.iterator().next().basesMatch(String.valueOf((char) a))) {
                byte temp = a;
                a = b;
                b = temp;

                int tmp = aIndex;
                aIndex = bIndex;
                bIndex = tmp;
                usedDbSNP = true;
            } else if (alts.size() == 1 && alts.iterator().next().basesMatch(String.valueOf((char) b))) {
                usedDbSNP = true;
            }
        }

      if (!usedDbSNP && generator.nextDouble() > 0.5) {

            byte temp = a;
            a = b;
            b = temp;

            int tmp = aIndex;
            aIndex = bIndex;
            bIndex = tmp;
        }
        if (counts[aIndex] == 0 && counts[bIndex] == 0)
            return 0;
        out.println(context.getLocation().getContig() + "," +
                context.getLocation().getStart() + "," + counts[aIndex] + "," + counts[bIndex] + "," + (char) a + "," + (char) b + "," + usedDbSNP);
        return 1;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer result) {
        out.close();
        logger.info("[REDUCE RESULT] Traversal result is: " + result);
    }
}
