package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;

/**
 * filter an indel callset based on given criteria
 */
public class IndelSubsets implements ConcordanceType {

    public final static int HOMOPOLYMER_MAX = 20;

    // Hint: to disable either of these features, set the value to zero (via command-line arguments);
    // then empty files will be output for the "below threshold" set, effectively disabling it
    private int homopolymerCutoff = 1;
    private int sizeCutoff = 2;

    private PrintWriter[][][][] writers = new PrintWriter[2][2][2][2];

    public IndelSubsets() {}

    public void initialize(String prefix, HashMap<String,String> args) {        
        if ( args.get("sizeCutoff") != null )
            sizeCutoff = Integer.valueOf(args.get("sizeCutoff"));
        if ( args.get("homopolymerCutoff") != null )
            homopolymerCutoff = Integer.valueOf(args.get("homopolymerCutoff"));

        try {
            for (int i = 0; i < 2; i++) {
                String name1 = i == 0 ? "set1" : "noSet1";
                for (int j = 0; j < 2; j++) {
                    String name2 = j == 0 ? "set2" : "noSet2";
                    for (int k = 0; k < 2; k++) {
                        String name3 = "size" + (k == 0 ? Integer.toString(sizeCutoff)+"AndUnder" : "GreaterThan"+Integer.toString(sizeCutoff));
                        for (int m = 0; m < 2; m++) {
                            String name4 = "homopolymer" + (m == 0 ? Integer.toString(homopolymerCutoff)+"AndUnder" : "GreaterThan"+Integer.toString(homopolymerCutoff));
                            writers[i][j][k][m] = new PrintWriter(prefix + "." + name1 + "." + name2 + "." + name3 + "." + name4 + ".calls");
                        }
                    }
                }
            }
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file(s) for writing"));
        }
    }

    public void computeConcordance(RefMetaDataTracker tracker, ReferenceContext ref) {
        Variation indel1 = (Variation)tracker.lookup("callset1", null);
        Variation indel2 = (Variation)tracker.lookup("callset2", null);

        int set1 = ( indel1 != null && indel1.isIndel() ? 0 : 1 );
        int set2 = ( indel2 != null && indel2.isIndel() ? 0 : 1 );

        // skip this locus if they're both not valid indels
        if ( set1 == 1 && set2 == 1 )
            return;

        // only deal with a valid indel
        Variation indel = ( indel1 != null ? indel1 : indel2 );

        // we only deal with the first allele
        int size = ( indel.getAlternateAlleleList().get(0).length() <= sizeCutoff ? 0 : 1 );
        int homopol = ( homopolymerRunSize(ref, indel) <= homopolymerCutoff ? 0 : 1 );

        writers[set1][set2][size][homopol].println(indel.toString());
    }

    private int homopolymerRunSize(ReferenceContext ref, Variation indel) {
        char[] bases = ref.getBases();
        GenomeLoc window = ref.getWindow();
        GenomeLoc locus = ref.getLocus();

        int refBasePos = (int)(locus.getStart() - window.getStart());
        char indelBase = indel.isDeletion() ? bases[refBasePos+1] : indel.getAlternateAlleleList().get(0).charAt(0);
        int leftRun = 0;
        for ( int i = refBasePos; i >= 0; i--) {
            if ( bases[i] != indelBase )
                break;
            leftRun++;
        }

        indelBase = indel.isDeletion() ? bases[Math.min(refBasePos+indel.getAlternateAlleleList().get(0).length(),bases.length-1)] : indel.getAlternateAlleleList().get(0).charAt(indel.getAlternateAlleleList().get(0).length()-1);
        int rightRun = 0;
        for ( int i = refBasePos + (indel.isDeletion() ? 1+indel.getAlternateAlleleList().get(0).length() : 1); i < bases.length; i++) {
            if ( bases[i] != indelBase )
                break;
            rightRun++;
        }

        //System.out.println(String.valueOf(bases) + ": " + leftRun + " / " + rightRun);
        return Math.max(leftRun, rightRun);
    }

    public void cleanup() {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    for (int m = 0; m < 2; m++)
                        writers[i][j][k][m].close();
    }
}