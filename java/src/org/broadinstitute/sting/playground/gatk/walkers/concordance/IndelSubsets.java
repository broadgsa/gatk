package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeCall;

import java.util.*;

/**
 * filter an indel callset based on given criteria
 */
public class IndelSubsets implements ConcordanceType {

    public final static int HOMOPOLYMER_MAX = 20;

    // Hint: to disable either of these features, set the value to zero (via command-line arguments);
    // then empty files will be output for the "below threshold" set, effectively disabling it
    private int homopolymerCutoff = 1;
    private int sizeCutoff = 2;

    private String[][][][] tags = new String[2][2][2][2];
    private String sample1, sample2;

    public IndelSubsets() {}

    public void initialize(Map<String, String> args, Set<String> samples) {
        if ( samples.size() != 2 )
            throw new StingException("IndelSubsets concordance test cannot handle anything other than 2 VCF records");

        if ( args.get("sizeCutoff") != null )
            sizeCutoff = Integer.valueOf(args.get("sizeCutoff"));
        if ( args.get("homopolymerCutoff") != null )
            homopolymerCutoff = Integer.valueOf(args.get("homopolymerCutoff"));

        Iterator<String> iter = samples.iterator();
        sample1 = iter.next();
        sample2 = iter.next();

        for (int i = 0; i < 2; i++) {
            String name1 = i == 0 ? sample1 : "no_" + sample1;
            for (int j = 0; j < 2; j++) {
                String name2 = j == 0 ? sample2 : "no_" + sample2;
                for (int k = 0; k < 2; k++) {
                    String name3 = "size" + (k == 0 ? Integer.toString(sizeCutoff)+"AndUnder" : "GreaterThan"+Integer.toString(sizeCutoff));
                    for (int m = 0; m < 2; m++) {
                        String name4 = "homopolymer" + (m == 0 ? Integer.toString(homopolymerCutoff)+"AndUnder" : "GreaterThan"+Integer.toString(homopolymerCutoff));
                        tags[i][j][k][m] = name1 + "." + name2 + "." + name3 + "." + name4;
                    }
                }
            }
        }
    }

    public String computeConcordance(Map<String, VCFGenotypeCall> samplesToRecords, ReferenceContext ref) {

        VCFGenotypeCall indel1 = samplesToRecords.get(sample1);
        VCFGenotypeCall indel2 = samplesToRecords.get(sample2);

        int set1 = ( indel1 != null && !indel1.isPointGenotype() ? 0 : 1 );
        int set2 = ( indel2 != null && !indel2.isPointGenotype() ? 0 : 1 );

        // skip this locus if they're both not valid indels
        if ( set1 == 1 && set2 == 1 )
            return null;

        // only deal with a valid indel
        Variation indel = ( indel1 != null ? indel1.toVariation() : indel2.toVariation() );

        // we only deal with the first allele
        int size = ( indel.getAlternateAlleleList().get(0).length() <= sizeCutoff ? 0 : 1 );
        int homopol = ( homopolymerRunSize(ref, indel) <= homopolymerCutoff ? 0 : 1 );

        return tags[set1][set2][size][homopol];
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

    public String getInfoName() { return "IndelSubsets"; }    
}