/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
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

    Hashtable Scores = new Hashtable();

    //HLAreads.add("Italian Riviera");

    public Pair<Long, Long> reduceInit() {
        try{
            out.printf("Reading HLA database ...");
            FileInputStream fstream = new FileInputStream(HLAdatabaseFile);
            // Get the object of DataInputStream
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine; String [] s = null;
            //Read File Line By Line
            while ((strLine = br.readLine()) != null)   {
                s = strLine.split("\\t");
                if (s.length>=10){
                    HLAreads.add(CigarFormatted(s[5],s[9]));
                    HLAcigars.add(s[5]);
                    HLAnames.add(s[0]);
                    HLApositions.add(s[3]);
                    //System.out.println (s[0] + "\t" + s[3]);
                }
            }
            in.close();
            int n = HLApositions.size();
            numHLAlleles = n;
            HLAstartpos = new int[n]; HLAstoppos = new int[n];
            Prob = new double[n][n];
            for (int i = 0; i < n; i++){
                HLAstartpos[i]=Integer.parseInt(HLApositions.get(i));
                HLAstoppos[i]=HLAstartpos[i]+HLAreads.get(i).length()-1;
                for (int j = 0; j <n; j++){
                    Prob[i][j]=0;
                }
                //System.out.println (HLAstartpos[i] + " to " + HLAstoppos[i]);
            }
            out.printf("DONE!\n");
        }catch (Exception e){//Catch exception if any
          System.err.println("Error: " + e.getMessage());
        }

        //out.printf("%s,%s\n",HLAnames.get(22),HLAnames.get(68));
        //String read = HLAreads.get(22);
        //String cigar = HLAcigars.get(22);
        //out.printf("%s,%s\n",cigar,read);
        //HLAreads.set(22, CigarFormatted(cigar,read));
        //HLAreads.set(68, CigarFormatted(HLAcigars.get(68),HLAreads.get(68)));
        out.printf("Comparing reads to database ...\n");
        int j1 = 80, k1 = 96, j2 = 80, k2 = 114;
        out.printf("%s,%s\t%s,%s\n",HLAnames.get(j1),HLAnames.get(k1),HLAnames.get(j2),HLAnames.get(k2));
        return new Pair<Long,Long>(0l,0l);
    }

    private String CigarFormatted(String cigar, String read){
        String formattedRead = "";
        char c;
        String count;
        int cigarPlaceholder = 0; int subcigarLength = 0;
        int readPlaceholder = 0; int subreadLength = 0;
        //out.printf("%s\n",cigar);
        for (int i = 0; i < cigar.length(); i++){
            c = cigar.charAt(i);
            if (c == 'M'){
                subcigarLength = i-cigarPlaceholder;
                count = cigar.substring(cigarPlaceholder, i);
                
                //read.substring(readPlaceholder, subreadLength)
                
                subreadLength = Integer.parseInt(count);
                formattedRead = formattedRead + read.substring(readPlaceholder, readPlaceholder+subreadLength);

                //out.printf("%sM,%s,%s\n", count,readPlaceholder,readPlaceholder+subreadLength);
                cigarPlaceholder = i+1;
                readPlaceholder = readPlaceholder + subreadLength;
                
            } else if (c == 'I'){
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);
                //formattedRead = formattedRead + read.substring(readPlaceholder, subreadLength);
                //out.printf("%sI,%s,%s\n", count,readPlaceholder,readPlaceholder+subreadLength);
                cigarPlaceholder = i+1;
                readPlaceholder = readPlaceholder + subreadLength;
                
            } else if (c == 'D'){
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);
                String deletion = "";
                for (int j = 1; j <= subreadLength; j++){
                    deletion = deletion + "D";
                }
                //out.printf("%sD,%s,%s\n", count,readPlaceholder,readPlaceholder+subreadLength);
                formattedRead = formattedRead + deletion;
                cigarPlaceholder = i+1;
                //readPlaceholder = readPlaceholder + subreadLength;
                
            }

        }
        //out.printf("%s\n",formattedRead);
        
        return formattedRead;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        
        //Create new context with mapping quality >= 20
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        GenomeLoc Gloc = context.getLocation();
        //List<SAMRecord> newreads = Arrays.asList();
        //List<Integer> newoffsets = Arrays.asList();

        //remove reads with mapping quality < 20

        //for (int i = 0; i < reads.size(); i++){
            //out.printf("i=%s,s=%s",i,reads.size());
            //SAMRecord r = reads.get(i);
            //if (r.getMappingQuality() < 20){
                //out.printf(",%s removed",i);
                //reads.remove(i);
                //offsets.remove(i);
                //i--;
            //}
            //out.printf("\n");
        //}
        
        //reset context to include only high quality reads
        context = new AlignmentContext(Gloc,reads,offsets);

        //Create pileup of reads at this locus using filtered context
        ReadBackedPileup pileup = new ReadBackedPileup(ref.getBase(), context);

        String bases = pileup.getBases();

        long loc = context.getPosition();

        if( context.getReads().size() > 0 ) {
            //out.printf("RG for first read: %s%n",context.getReads().get(0).getReadName());
            int numAs = 0;
            int numCs = 0;
            int numGs = 0;
            int numTs = 0;
            String c1 = ""; String c2 = "";
            long pos_k = 0; long pos_j = 0;
            
            for (int i = 0;i < context.getReads().size();i++){
                //out.printf("%s%s ",i,bases.charAt(i));
                char base = bases.charAt(i);
                if (base == 'A'){numAs++;}
                if (base == 'C'){numCs++;}
                if (base == 'T'){numTs++;}
                if (base == 'G'){numGs++;}

            }
            //out.printf("%s\t", context.getLocation());
            //out.printf("ref=%s\t", ref);
            //out.printf("%sAs\t%sCs\t%sTs\t%sGs\t",numAs,numCs,numTs,numGs);
            int j1 = 80, k1 = 96, j2 = 80, k2 = 114;


            //Calculate posterior probabilities!
            NewHotnessGenotypeLikelihoods G = new NewHotnessGenotypeLikelihoods(GenotypeLikelihoods.HUMAN_HETEROZYGOSITY);
            //SSGGenotypeCall geno = G.callGenotypes(tracker, ref, pileup);

            for (int i = 0; i < pileup.getReads().size(); i++) {
                SAMRecord read = pileup.getReads().get(i);
                int offset = pileup.getOffsets().get(i);
                char base = read.getReadString().charAt(offset);
                byte qual = read.getBaseQualities()[offset];
                if (qual >= 20){
                    G.add(ref.getBase(), base, qual);
                }
            }

            double mLikelihoods[] = G.likelihoods;

            for (int i = 0; i< mLikelihoods.length; i++){
            //    out.printf("%5.0f\t",mLikelihoods[i]);
            }

            //Store confidence scores
            Scores = new Hashtable();
            Scores.put("AA", mLikelihoods[0]);
            Scores.put("AC", mLikelihoods[1]);
            Scores.put("AG", mLikelihoods[2]);
            Scores.put("AT", mLikelihoods[3]);
            Scores.put("CC", mLikelihoods[4]);
            Scores.put("CG", mLikelihoods[5]);
            Scores.put("CT", mLikelihoods[6]);
            Scores.put("GG", mLikelihoods[7]);
            Scores.put("GT", mLikelihoods[8]);
            Scores.put("TT", mLikelihoods[9]);

            

            //Update probabilities for combinations of alleles
            //For each HLA allele
            String r1 = "", r2 = "", s1 = "", s2 = "";
            j1 = 80; k1 = 96; j2 = 80; k2 = 114;
            for (int j = 0; j < numHLAlleles; j++){
                //out.print(j+","+loc + "," + HLAstartpos[j] + "," + HLAstoppos[j]);
                //check if allele overlaps current position
                if (loc >= HLAstartpos[j] && loc <= HLAstoppos[j]){
                    pos_j = loc - HLAstartpos[j];
                    c1 = Character.toString(Character.toUpperCase(HLAreads.get(j).charAt((int) pos_j)));

                    if (j == j1){ //C*010201
                        r1 = c1;
                    }
                    if (j == k1){ //C*070101
                        r2 = c1;
                    }
                    if(j == j2){ //C*010201
                        s1 = c1;
                    }
                    if (j == k2){ //C*150201
                        s2 = c1;
                    }
                    
                    for (int k = 0; k < numHLAlleles; k++){
                        //out.print(k+","+loc + "," + HLAstartpos[j] + "," + HLAstoppos[j]);
                        if (loc >= HLAstartpos[k] && loc <= HLAstoppos[k]){
                            pos_k = loc - HLAstartpos[k];
                            c2 = Character.toString(Character.toUpperCase(HLAreads.get(k).charAt((int) pos_k)));
                            //out.printf("j=%s,k=%s,C1=%s,C2=%s,",j,k,c1,c2);
                            if (Scores.containsKey(c1 + c2)){
                                String base=c1 + c2;
                                Prob[j][k]= Prob[j][k]+Double.parseDouble((String) Scores.get(base));
                                //out.printf("P(%s)=%s\t",base,Prob[j][k]);
                            }else if(Scores.containsKey(c2 + c1)){
                                String base=c2 + c1;
                                Prob[j][k]= Prob[j][k]+Double.parseDouble((String) Scores.get(base));
                                //out.printf("P(%s)=%s\t",base,Prob[j][k]);
                            }
                        
                        }
                    }

                }
            }

            //out.printf("%s%s:%5.0f\t%s%s:%5.0f\t",r1,r2,Prob[j1][k1],s1,s2,Prob[j2][k2]);
            //out.printf("%s %5.0f\t",base,score);

            if ( !suppressPrinting ){
            }
           
            //out.printf("\t");
            //mGenotypes.get(0).getConfidenceScore().
            //out.printf("RG for first read: %s%n",context.getReads().get(0).getReadName());
        }
        return context.getReads().size();

    }

    private double[] priorsArray(String priorsString) {
        String[] pstrs = priorsString.split(",");
        double[] pdbls = new double[pstrs.length];

        for (int i = 0; i < pstrs.length; i++) {
            pdbls[i] = Double.valueOf(pstrs[i]);
        }

        return pdbls;
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
        double max = 0; int max_i =0; int max_j = 0;
        double max_2 = 0; int max_i_2 =0; int max_j_2 = 0;
        double max_3 = 0; int max_i_3 =0; int max_j_3 = 0;

        out.printf("\n");
        for (int i = 0; i < numHLAlleles; i++){
            for (int j = 0; j < numHLAlleles; j++){
                if (i >= j ){
                    out.printf("%s\t%s\t%5.0f\n",HLAnames.get(i),HLAnames.get(j),Prob[i][j]);
                }
                if (Prob[i][j] > max){
                    max_3 = max_2; max_i_3 = max_i_2; max_j_3 = max_j_2;
                    max_2 = max;   max_i_2 = max_i;   max_j_2 = max_j;
                    max = Prob[i][j]; max_i = i; max_j = j;
                }
            }
        }
        out.printf("\n");
        out.printf("Highest score: %5.0f in %s and %s\n",max,HLAnames.get(max_i),HLAnames.get(max_j));
        out.printf("2nd Highest score: %5.0f in %s and %s\n",max_2,HLAnames.get(max_i_2),HLAnames.get(max_j_2));
        out.printf("3rd Highest score: %5.0f in %s and %s\n",max_3,HLAnames.get(max_i_3),HLAnames.get(max_j_3));
        //out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n",((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
    }
}
