/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;


import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
/**
 *
 * @author shermanjia
 */
public class CallHLAWalker extends LocusWalker<Integer, Pair<Long, Long>>{
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false)
    public boolean suppressPrinting = false;

    //String HLAdatabaseFile = "/Users/shermanjia/Work/HLA.sam";
    //String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.sam";
    //String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.sam";
    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.4digitUnique.sam";



    ArrayList<String> HLAreads = new ArrayList<String>();
    ArrayList<String> HLAcigars = new ArrayList<String>();
    ArrayList<String> HLAnames = new ArrayList<String>();
    ArrayList<String> HLApositions = new ArrayList<String>();
    
    ArrayList<SAMRecord> AllReads = new ArrayList<SAMRecord>();
    ArrayList<String> AllReadNames = new ArrayList<String>();

    int[] HLAstartpos;
    int[] HLAstoppos;
    int numHLAlleles = 0;
    int numInterval = 1;
    double[][] Prob; String[][] Alleles;
    int j1, k1, j2, k2;
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;


    Hashtable Scores = new Hashtable();

    boolean DatabaseLoaded = false;
    boolean PrintedOutput = false;
    boolean DEBUG = true;

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
                int i = 0;
                while ((strLine = br.readLine()) != null)   {
                    s = strLine.split("\\t");
                    if (s.length>=10){
                        //Parse the reads with cigar parser
                        HLAreads.add(CigarFormatted(s[5],s[9]));
                        HLAcigars.add(s[5]);
                        HLAnames.add(s[0]);
                        HLApositions.add(s[3]);
                        if (s[0].indexOf("HLA_A") > -1){
                            if (iAstart < 0){iAstart=i;}
                            iAstop = i; i++;
                        }else if (s[0].indexOf("HLA_B") > -1){
                            if (iBstart < 0){iBstart=i;}
                            iBstop = i; i++;
                        }else if (s[0].indexOf("HLA_C") > -1){
                            if (iCstart < 0){iCstart=i;}
                            iCstop = i; i++;
                        }
                    }
                }
                in.close();
                int n = HLApositions.size(); numHLAlleles = n;
                HLAstartpos = new int[n]; HLAstoppos = new int[n];
                Prob = new double[n][n];

                for (i = 0; i < n; i++){
                    //Find start and stop positions for each allele
                    HLAstartpos[i]=Integer.parseInt(HLApositions.get(i));
                    HLAstoppos[i]=HLAstartpos[i]+HLAreads.get(i).length()-1;
                    //Initialize matrix of probabilities / likelihoods
                    for (int j = 0; j <n; j++){
                        Prob[i][j]=0;
                    }
                    //For debugging: get index for specific alleles
                    if (HLAnames.get(i).equals("HLA_A*110101"))
                        j1 = i;
                    if (HLAnames.get(i).equals("HLA_A*01010101"))
                        k1 = i;
                    if (HLAnames.get(i).equals("HLA_A*110201"))
                        j2 = i;
                    if (HLAnames.get(i).equals("HLA_A*0109"))
                        k2 = i;
                }
                out.printf("DONE!\n");
            }catch (Exception e){//Catch exception if any
              System.err.println("Error: " + e.getMessage());
            }

            DatabaseLoaded = true;
            out.printf("Comparing reads to database ...\n");

            //For debugging: prints which HLA alleles were indexed before
            if (j1 > k1){int tmp = k1; k1 = j1; j1 = tmp;}
            if (j2 > k2){int tmp = k2; k2 = j2; j2 = tmp;}

            if (DEBUG){
                out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
                out.printf("%s,%s\t%s,%s\n",HLAnames.get(j1),HLAnames.get(k1),HLAnames.get(j2),HLAnames.get(k2));
            }
        }
        out.printf("Computing for interval %s...\n",numInterval);
        numInterval++;
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
            int depth = 0;
            String c1 = ""; String c2 = "";
            long pos_k = 0; long pos_j = 0;
            
            //Debugging purposes: print location, reference base, pileup, and count (before quality filtering)
            if (DEBUG){
                out.printf("%s\t", context.getLocation());
                out.printf("ref=%s\t", ref.getBase());
                //out.printf("%s\t",pileup.getBases().toString());
                //out.printf("%s\t",pileup.getBasePileupAsCountsString());
            }

            //Calculate posterior probabilities!


            GenotypeLikelihoods G = new ThreeStateErrorGenotypeLikelihoods();


            //I've tried simply adding the entire pileup to g, but quality scores are not checked, and 'N' bases throw an error
            //(GenotypeLikelihoods potentially has bug on line 57)
            //g.add(pileup);

            //Check for bad bases and ensure mapping quality myself. This works.
            for (int i = 0; i < context.getReads().size(); i++) {
                
                SAMRecord read = context.getReads().get(i);
                int offset = context.getOffsets().get(i);
                char base = read.getReadString().charAt(offset);
                byte qual = read.getBaseQualities()[offset];
                int mapquality = read.getMappingQuality();
                if (mapquality >= 5 && BaseUtils.simpleBaseToBaseIndex(base) != -1) {

                    //Adds read to list of reads for final phasing comparison
                    String name = read.getReadName();
                    if (!AllReadNames.contains(name)){
                        AllReadNames.add(name);
                        AllReads.add(read);
                    }

                    if (base == 'A'){numAs++; depth++;}
                    if (base == 'C'){numCs++; depth++;}
                    if (base == 'T'){numTs++; depth++;}
                    if (base == 'G'){numGs++; depth++;}
                    //consider base in likelihood calculations if it looks good and has high mapping score
                    G.add(base, qual, read, offset);
                }
            }

            //Debugging purposes
            if (DEBUG) {out.printf("A[%s]\tC[%s]\tT[%s]\tg[%s]\t",numAs,numCs,numTs,numGs);}

            if (depth > 0){
                //Store confidence scores - this is a local hash that we use to get likelihood given a particular genotype
                Scores = new Hashtable();
                for ( DiploidGenotype g : DiploidGenotype.values() ) {
                    Scores.put(g.toString(), G.getLikelihood(g));
                    //also hash other combination not stored by DiploidGenotype
                    if (g.toString().equals("AC")) Scores.put("CA", G.getLikelihood(g));
                    if (g.toString().equals("AG")) Scores.put("GA", G.getLikelihood(g));
                    if (g.toString().equals("AT")) Scores.put("TA", G.getLikelihood(g));
                    if (g.toString().equals("CG")) Scores.put("GC", G.getLikelihood(g));
                    if (g.toString().equals("CT")) Scores.put("TC", G.getLikelihood(g));
                    if (g.toString().equals("GT")) Scores.put("TG", G.getLikelihood(g));

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

                        //Only check A-A, B-C, C-C combinations
                        int kStart = 0, kStop = 0;
                        if (j >= iAstart && j <= iAstop){
                            kStart = iAstart; kStop = iAstop;
                        } else if (j >= iBstart && j <= iBstop){
                            kStart = iBstart; kStop = iBstop;
                        } else if (j >= iCstart && j <= iCstop){
                            kStart = iCstart; kStop = iCstop;
                        }

                        //Fill half-matrix only to speed up process
                        if (j > kStart){kStart = j;}

                        if (DEBUG){
                            //out.printf("j[%s],k[%s,%s]\t",j,kStart,kStop);
                        }

                        //Update likelihoods
                        for (int k = kStart; k <= kStop; k++){

                            //check if allele 2 overlaps current position
                            if (loc >= HLAstartpos[k] && loc <= HLAstoppos[k]){
                                pos_k = loc - HLAstartpos[k];
                                c2 = Character.toString(Character.toUpperCase(HLAreads.get(k).charAt((int) pos_k)));

                                //updates likelihoods for both permutations of the alleles, normalized to the likelihood for homozygous reference
                                if (!homref.equals(c1+c2)){
                                    if (Scores.containsKey(c1 + c2)){
                                        //out.printf("j[%s],k[%s],g[%s],s[%.0f]\t",j,k,c1+c2,Scores.get(c1 + c2));
                                        likelihood = Double.parseDouble((String) Scores.get(c1 + c2).toString());
                                        Prob[j][k]= Prob[j][k]+ likelihood - homreflikelihood;
                                    } else{
                                        if (DEBUG){
                                        //out.printf("\nCharacters [%s] not found,j[%s],k[%s],%s,%s\n",c1+c2,j,k,HLAnames.get(j),HLAnames.get(k));
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            

                if ( DEBUG ){
                    //Debugging: print updated likelihoods for 2 sets of HLA alleles, as well as normalized likelihoods for all 10 genotypes
                    out.printf("%s%s:%5.0f\t%s%s:%5.0f\t",r1,r2,Prob[j1][k1],s1,s2,Prob[j2][k2]);
                    for ( DiploidGenotype g : DiploidGenotype.values() ) {
                        out.printf("%s %5.0f\t",g.toString(),likelihood - homreflikelihood);
                    }
                    out.printf("\n");
                }
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
            ArrayList<Integer> TopAlleles = new ArrayList<Integer>();

            double maxA = 0;  int i_maxA =0;   int j_maxA = 0;
            double maxA2 = 0; int i_maxA_2 =0; int j_maxA_2 = 0;

            double maxB = 0;  int i_maxB =0;   int j_maxB = 0;
            double maxB2 = 0; int i_maxB_2 =0; int j_maxB_2 = 0;

            double maxC = 0;  int i_maxC =0;   int j_maxC = 0;
            double maxC2 = 0; int i_maxC_2 =0; int j_maxC_2 = 0;

            out.printf("\n");

            //Find max scores for each HLA gene.
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = i; j < numHLAlleles; j++){
                    //Print likelihoods for all alleles
                    if (!DEBUG){
                        out.printf("%s\t%s\t%5.0f\n",HLAnames.get(i),HLAnames.get(j),Prob[i][j]);
                    }

                    if (HLAnames.get(i).indexOf("HLA_A") > -1 && HLAnames.get(j).indexOf("HLA_A") > -1){
                        if (Prob[i][j] > maxA){
                            maxA2 = maxA; i_maxA_2 = i_maxA; j_maxA_2 = j_maxA;
                            maxA = Prob[i][j]; i_maxA = i; j_maxA = j;
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_B") > -1 && HLAnames.get(j).indexOf("HLA_B") > -1){
                        if (Prob[i][j] > maxB){
                            maxB2 = maxB; i_maxB_2 = i_maxB; j_maxB_2 = j_maxB;
                            maxB = Prob[i][j]; i_maxB = i; j_maxB = j;
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_C") > -1 && HLAnames.get(j).indexOf("HLA_C") > -1){
                        if (Prob[i][j] > maxC){
                            maxC2 = maxC; i_maxC_2 = i_maxC; j_maxC_2 = j_maxC;
                            maxC = Prob[i][j]; i_maxC = i; j_maxC = j;
                        }
                    }
		}
            }


            //Print alleles with the highest likelihoods
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = i; j < numHLAlleles; j++){
                    if (HLAnames.get(i).indexOf("HLA_A") > -1 && HLAnames.get(j).indexOf("HLA_A") > -1){
                        if (Prob[i][j] == maxA && maxA > 0){
                            out.printf("%s\t%s\t%5.0f\t%5.0f\tBEST\n",HLAnames.get(i),HLAnames.get(j),maxA,maxA-maxA2);
                            if (!TopAlleles.contains(i)){TopAlleles.add(i);}
                            if (!TopAlleles.contains(j)){TopAlleles.add(j);}
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_B") > -1 && HLAnames.get(j).indexOf("HLA_B") > -1){
                        if (Prob[i][j] == maxB && maxB > 0){
                            out.printf("%s\t%s\t%5.0f\t%5.0f\tBEST\n",HLAnames.get(i),HLAnames.get(j),maxB,maxB-maxB2);
                            if (!TopAlleles.contains(i)){TopAlleles.add(i);}
                            if (!TopAlleles.contains(j)){TopAlleles.add(j);}
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_C") > -1 && HLAnames.get(j).indexOf("HLA_C") > -1){
                        if (Prob[i][j] == maxC && maxC > 0){
                            out.printf("%s\t%s\t%5.0f\t%5.0f\tBEST\n",HLAnames.get(i),HLAnames.get(j),maxC,maxC-maxC2);
                            if (!TopAlleles.contains(i)){TopAlleles.add(i);}
                            if (!TopAlleles.contains(j)){TopAlleles.add(j);}
                        }
                    }

                }
            }

            //Check top alleles for phasing
            for (int i = 0; i< TopAlleles.size(); i++){
                if ( DEBUG ){
                //Debugging: print list of alleles to be checked for phasing
                    out.printf("%s\t%s",TopAlleles.get(i),HLAnames.get(TopAlleles.get(i)));
                    out.printf("\n");
                }
                
            }

            //2nd Highest likelihoods
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = i; j < numHLAlleles; j++){
                    if (Prob[i][j] == maxA2){
                        //out.printf("2nd Highest likelihood: %5.0f in %s and %s; i=%s, j=%s\n",maxA2,HLAnames.get(i),HLAnames.get(j),i,j);
                    }
                }
            }

            PrintedOutput = true;
            //out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n",((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
        }
    }
}
