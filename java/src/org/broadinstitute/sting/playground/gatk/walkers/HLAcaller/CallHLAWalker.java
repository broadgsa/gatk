/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import java.io.*;
import java.util.*;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
/**
 *
 * @author shermanjia
 */
public class CallHLAWalker extends LocusWalker<Integer, Pair<Long, Long>>{
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false)
    public boolean suppressPrinting = false;

    //String HLAdatabaseFile = "/Users/shermanjia/Work/HLA.sam";
    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.sam";


    ArrayList<String> HLAreads = new ArrayList<String>();
    ArrayList<String> HLAcigars = new ArrayList<String>();
    ArrayList<String> HLAnames = new ArrayList<String>();
    ArrayList<String> HLApositions = new ArrayList<String>();
    int[] HLAstartpos;
    int[] HLAstoppos;
    int numHLAlleles = 0;
    double[][] Prob; String[][] Alleles;
    int j1, k1, j2, k2;

    Hashtable Scores = new Hashtable();

    boolean DatabaseLoaded = false;
    boolean PrintedOutput = false;

    public Pair<Long, Long> reduceInit() {
        //Load sequences corresponding to HLA alleles from sam file
        if (!DatabaseLoaded){
            try{
                out.printf("Reading HLA database ...");
                FileInputStream fstream = new FileInputStream(HLAdatabaseFile);
                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine; String [] s = null;
                //Read File Line By Line
                while ((strLine = br.readLine()) != null)   {
                    s = strLine.split("\\t");
                    if (s.length>=10){
                        //Parse the reads with cigar parser
                        HLAreads.add(CigarFormatted(s[5],s[9]));
                        HLAcigars.add(s[5]);
                        HLAnames.add(s[0]);
                        HLApositions.add(s[3]);
                    }
                }
                in.close();
                int n = HLApositions.size(); numHLAlleles = n;
                HLAstartpos = new int[n]; HLAstoppos = new int[n];
                Prob = new double[n][n];

                for (int i = 0; i < n; i++){
                    //Find start and stop positions for each allele
                    HLAstartpos[i]=Integer.parseInt(HLApositions.get(i));
                    HLAstoppos[i]=HLAstartpos[i]+HLAreads.get(i).length()-1;
                    //Initialize matrix of probabilities / likelihoods
                    for (int j = 0; j <n; j++){
                        Prob[i][j]=0;
                    }
                    //For debugging: get index for specific alleles
                    if (HLAnames.get(i).equals("HLA_A*310102"))
                        j1 = i;
                    if (HLAnames.get(i).equals("HLA_A*320101"))
                        k1 = i;
                    if (HLAnames.get(i).equals("HLA_A*330301"))
                        j2 = i;
                    if (HLAnames.get(i).equals("HLA_A*0265"))
                        k2 = i;
                }
                out.printf("DONE!\n");
            }catch (Exception e){//Catch exception if any
              System.err.println("Error: " + e.getMessage());
            }

            DatabaseLoaded = true;
            out.printf("Comparing reads to database ...\n");

            //For debugging: prints which HLA alleles were indexed before
            //out.printf("%s,%s\t%s,%s\n",HLAnames.get(j1),HLAnames.get(k1),HLAnames.get(j2),HLAnames.get(k2));
        }
        return new Pair<Long,Long>(0l,0l);
    }

    private String CigarFormatted(String cigar, String read){
        // returns a cigar-formatted sequence (removes insertions, inserts 'D' to where deletions occur
        String formattedRead = ""; char c; String count;
        int cigarPlaceholder = 0; int subcigarLength = 0;
        int readPlaceholder = 0; int subreadLength = 0;
        
        //reads cigar string
        for (int i = 0; i < cigar.length(); i++){
            c = cigar.charAt(i);
            if (c == 'M'){
                //If reach M for match/mismatch, get number immediately preceeding 'M' and tack on that many characters to sequence
                subcigarLength = i-cigarPlaceholder;
                count = cigar.substring(cigarPlaceholder, i);
                
                subreadLength = Integer.parseInt(count);
                formattedRead = formattedRead + read.substring(readPlaceholder, readPlaceholder+subreadLength);

                //increment placeholders
                cigarPlaceholder = i+1;
                readPlaceholder = readPlaceholder + subreadLength;                
            } else if (c == 'I'){
                //***NOTE: To be modified later if needed (insertions removed here)***

                //If reaches I for insertion, get number before 'I' and skip that many characters in sequence
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);

                //increment placeholders without adding inserted bases to sequence (effectively removes insertion).
                cigarPlaceholder = i+1;
                readPlaceholder = readPlaceholder + subreadLength;               
            } else if (c == 'D'){
                //If reaches D for deletion, insert 'D' into sequence as placeholder
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);

                //Add one 'D' for each deleted base
                String deletion = "";
                for (int j = 1; j <= subreadLength; j++){
                    deletion = deletion + "D";
                }

                //update placeholders
                formattedRead = formattedRead + deletion;
                cigarPlaceholder = i+1;
            }

        }
        return formattedRead;
    }

    private static int unsignedByteToInt(byte b) {
        //converts base quality from byte to int (not really needed)
        return (int) b & 0xFF;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        GenomeLoc Gloc = context.getLocation();

        //Create pileup of reads at this locus
        ReadBackedPileup pileup = new ReadBackedPileup(ref.getBase(), context);

        long loc = context.getPosition();

        if( context.getReads().size() > 0 ) {
            //out.printf("RG for first read: %s%n",context.getReads().get(0).getReadName());
            int numAs = 0;
            int numCs = 0;
            int numGs = 0;
            int numTs = 0;
            String c1 = ""; String c2 = "";
            long pos_k = 0; long pos_j = 0;
            
            //Debugging purposes: print location, reference base, pileup, and count (before quality filtering)
            //out.printf("%s\t", context.getLocation());
            //out.printf("ref=%s\t", ref.getBase());
            //out.printf("%s\t",pileup.getBases().toString());
            //out.printf("%s\t",pileup.getBasePileupAsCountsString());


            //Calculate posterior probabilities!
            NewHotnessGenotypeLikelihoods G = new NewHotnessGenotypeLikelihoods();

            //I've tried simply adding the entire pileup to G, but quality scores are not checked, and 'N' bases throw an error
            //(NewHotnessGenotypeLikelihoods potentially has bug on line 57)
            //G.add(pileup);

            //Check for bad bases and ensure mapping quality myself. This works.
            for (int i = 0; i < context.getReads().size(); i++) {
                SAMRecord read = context.getReads().get(i);
                int offset = context.getOffsets().get(i);
                char base = read.getReadString().charAt(offset);
                byte qual = read.getBaseQualities()[offset];
                int mapquality = read.getMappingQuality();
                if (mapquality >= 20 && BaseUtils.simpleBaseToBaseIndex(base) != -1) {
                    if (base == 'A'){numAs++;}
                    if (base == 'C'){numCs++;}
                    if (base == 'T'){numTs++;}
                    if (base == 'G'){numGs++;}
                    //consider base in likelihood calculations if it looks good and has high mapping score
                    G.add(base, qual);
                }
            }

            //Debugging purposes
            //out.printf("A[%s]\tC[%s]\tT[%s]\tG[%s]\t",numAs,numCs,numTs,numGs);


            //Store confidence scores - this is a local hash that we use to get likelihood given a particular genotype
            Scores = new Hashtable();
            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                Scores.put(g.toString(), G.getLikelihood(g));
            }

            //Get likelihood score for homozygous ref: used to normalize likelihoood scores at 0.
            String homref = String.valueOf(ref.getBase())+String.valueOf(ref.getBase());
            Double homreflikelihood = Double.parseDouble((String) Scores.get(homref).toString());


            //Update likelihood for each combinations of alleles
            Double likelihood = 0.0;
            String r1 = "", r2 = "", s1 = "", s2 = "";
            for (int j = 0; j < numHLAlleles; j++){

                //check if allele 1 overlaps current position
                if (loc >= HLAstartpos[j] && loc <= HLAstoppos[j]){
                    pos_j = loc - HLAstartpos[j];
                    c1 = Character.toString(Character.toUpperCase(HLAreads.get(j).charAt((int) pos_j)));

                    //Extract bases for HLA alleles indicated in reduceInit (for debugging)
                    if (j == j1) r1 = c1;
                    if (j == k1) r2 = c1;
                    if (j == j2) s1 = c1;
                    if (j == k2) s2 = c1;

                    //Updade likelihoods
                    for (int k = 0; k < numHLAlleles; k++){

                        //check if allele 2 overlaps current position
                        if (loc >= HLAstartpos[k] && loc <= HLAstoppos[k]){
                            pos_k = loc - HLAstartpos[k];
                            c2 = Character.toString(Character.toUpperCase(HLAreads.get(k).charAt((int) pos_k)));

                            //updates likelihoods for both permutations of the alleles, normalized to the likelihood for homozygous reference
                            
                            if (Scores.containsKey(c1 + c2)){

                                likelihood = Double.parseDouble((String) Scores.get(c1 + c2).toString());
                                Prob[j][k]= Prob[j][k]+ likelihood - homreflikelihood;

                            }else if(Scores.containsKey(c2 + c1)){

                                likelihood = Double.parseDouble((String) Scores.get(c2 + c1).toString());
                                Prob[j][k]= Prob[j][k]+ likelihood - homreflikelihood;

                            }
                        }
                    }
                }
            }

            if ( !suppressPrinting ){
                //Debugging: print updated likelihoods for 2 sets of HLA alleles, as well as normalized likelihoods for all 10 genotypes
                //out.printf("%s%s:%5.0f\t%s%s:%5.0f\t",r1,r2,Prob[j1][k1],s1,s2,Prob[j2][k2]);
                for ( DiploidGenotype g : DiploidGenotype.values() ) {
                    //out.printf("%s %5.0f\t",g.toString(),likelihood - homreflikelihood);
                }
                //out.printf("\n");
            }           
        }
        
        return context.getReads().size();
    }

    public boolean isReduceByInterval() {
        return true;
    }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = value.longValue() + sum.getFirst();
        long right = sum.getSecond() + 1l;
        return new Pair<Long,Long>(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {

        //Print HLA allele combinations with highest likelihood sums
        if (!PrintedOutput){
            double max = 0; int i_max =0; int j_max = 0;
            double m2 = 0; int i_max_2 =0; int j_max_2 = 0;
            double m3 = 0; int i_max_3 =0; int j_max_3 = 0;

            out.printf("\n");

            //Find max scores.
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = 0; j < numHLAlleles; j++){
                    if (i >= j){
                        //Print likelihoods for all alleles
                        out.printf("%s\t%s\t%5.0f",HLAnames.get(i),HLAnames.get(j),Prob[i][j]);
                        if (Prob[i][j] > max){
                            m3 = m2;  i_max_3 = i_max_2; j_max_3 = j_max_2;
                            m2 = max; i_max_2 = i_max;   j_max_2 = j_max;
                            max = Prob[i][j]; i_max = i; j_max = j;
                        }
                    }
                }
            }

            //Highest likelihoods
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = 0; j < numHLAlleles; j++){
                    if (i >= j){
                        if (Prob[i][j] == max){
                            out.printf("Highest likelihood: %5.0f in %s and %s; i=%s, j=%s\n",max,HLAnames.get(i),HLAnames.get(j),i,j);

                        }
                    }
                }
            }

            //2nd Highest likelihoods
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = 0; j < numHLAlleles; j++){
                    if (i >= j){
                        if (Prob[i][j] == m2){
                            out.printf("2nd Highest likelihood: %5.0f in %s and %s; i=%s, j=%s\n",m2,HLAnames.get(i),HLAnames.get(j),i,j);
                        }
                    }
                }
            }

            PrintedOutput = true;
            //out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n",((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
        }
    }
}
