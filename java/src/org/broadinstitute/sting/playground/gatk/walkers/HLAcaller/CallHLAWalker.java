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
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
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
    //String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/CEU_Founders_HLA.freq";

    String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";


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
    double[] SingleAlleleFrequencies;
    double[][] LOD; String[][] Alleles;
    double[][] ActualLikelihoods;
    int[][] PhasingScores;
    double[][] CombinedAlleleFrequencies;
    int j1, k1, j2, k2;
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;

    Hashtable AlleleFrequencies = new Hashtable();

    Hashtable Scores = new Hashtable();
    Hashtable SNPs = new Hashtable();
    ArrayList<Integer> SNPlocations = new ArrayList<Integer>();
    Integer SNPcount = 0;
    int[][] SNPcorrelation;

    boolean DatabaseLoaded = false;
    boolean PrintedOutput = false;
    boolean DEBUG = false;

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
                SingleAlleleFrequencies = new double[n];
                LOD = new double[n][n];
                ActualLikelihoods = new double[n][n];
                PhasingScores = new int[n][n];
                CombinedAlleleFrequencies = new double[n][n];

                for (i = 0; i < n; i++){
                    //Find start and stop positions for each allele
                    HLAstartpos[i]=Integer.parseInt(HLApositions.get(i));
                    HLAstoppos[i]=HLAstartpos[i]+HLAreads.get(i).length()-1;
                    SingleAlleleFrequencies[i]=0.0;
                    //Initialize matrix of probabilities / likelihoods
                    for (int j = 0; j <n; j++){
                        LOD[i][j]=0;
                        ActualLikelihoods[i][j]=0;
                        PhasingScores[i][j]=0;
                        CombinedAlleleFrequencies[i][j]=0.0;
                    }
                    //For debugging: get index for specific alleles
                    if (HLAnames.get(i).equals("HLA_B*0801"))
                        j1 = i;
                    if (HLAnames.get(i).equals("HLA_B*5601"))
                        k1 = i;
                    if (HLAnames.get(i).equals("HLA_B*0813"))
                        j2 = i;
                    if (HLAnames.get(i).equals("HLA_B*5613"))
                        k2 = i;
                }
                out.printf("DONE! Read %s alleles\n",HLAreads.size());
            }catch (Exception e){//Catch exception if any
              System.err.println("Error: " + e.getMessage());
            }

            try{
                out.printf("Reading allele frequences ...");
                FileInputStream fstream = new FileInputStream(CaucasianAlleleFrequencyFile);
                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine; String [] s = null;
                //Read File Line By Line
                int count = 0;
                while ((strLine = br.readLine()) != null)   {
                    s = strLine.split("\\t");
                    AlleleFrequencies.put(s[0], s[1]);
                    count++;
                }
                in.close();
                out.printf("Done! Read %s alleles\n",count);
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
            } else if (c == 'H'){
                //(Headers removed here)***

                //If reaches H for insertion, get number before 'H' and skip that many characters in sequence
                count = cigar.substring(cigarPlaceholder, i);
                subreadLength = Integer.parseInt(count);

                //increment cigar placeholder without adding inserted bases to sequence (effectively removes insertion).
                cigarPlaceholder = i+1;
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

            //Check for bad bases and ensure mapping quality myself. This works.
            for (int i = 0; i < context.getReads().size(); i++) {
                
                SAMRecord read = context.getReads().get(i);
                int offset = context.getOffsets().get(i);
                char base = read.getReadString().charAt(offset);
                byte qual = read.getBaseQualities()[offset];
                int mapquality = read.getMappingQuality();
                if (mapquality >= 5 && BaseUtils.simpleBaseToBaseIndex(base) != -1) {

                    String name = read.getReadName();
                    if (!AllReadNames.contains(name)){
                        AllReadNames.add(name);
                        AllReads.add(read);
                        if( DEBUG){
                            //out.print("\n" + read.getReadName() + "\t" + read.getCigarString() + "\t" + read.getReadLength() + "\n");
                        }
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

                Double likelihood = 0.0;
                
                for ( DiploidGenotype g : DiploidGenotype.values() ) {
                    likelihood = G.getLikelihood(g);

                    Scores.put(g.toString(), likelihood);
                    //also hash other combination not stored by DiploidGenotype
                    if (g.toString().equals("AC")) {
                        Scores.put("CA", likelihood);
                    } else if (g.toString().equals("AG")){
                        Scores.put("GA", likelihood);
                    } else if (g.toString().equals("AT")){
                        Scores.put("TA", likelihood);
                    } else if (g.toString().equals("CG")){
                        Scores.put("GC", likelihood);
                    } else if (g.toString().equals("CT")){
                        Scores.put("TC", likelihood);
                    } else if (g.toString().equals("GT")){
                        Scores.put("TG", likelihood);
                    }
                }

                //Get likelihood score for homozygous ref: used to normalize likelihoood scores at 0.
                String homref = String.valueOf(ref.getBase())+String.valueOf(ref.getBase());
                Double homreflikelihood = Double.parseDouble((String) Scores.get(homref).toString());

                //Add SNP if it is a SNP and hasn't been added before
                for ( DiploidGenotype g : DiploidGenotype.values() ) {
                    likelihood = G.getLikelihood(g);
                    if (likelihood > homreflikelihood && !SNPs.containsKey(Long.toString(loc))){
                        SNPcount++; SNPs.put(Long.toString(loc),SNPcount); SNPlocations.add(Integer.valueOf(Long.toString(loc)));
                        
                    }
                }

                //Update likelihood for each combinations of alleles
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

                        //Only check HLA A-A, B-B, C-C combinations
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
                                        ActualLikelihoods[j][k] = LOD[j][k]+ likelihood;
                                        LOD[j][k]= ActualLikelihoods[j][k] - homreflikelihood;

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
                    out.printf("%s%s:%5.0f\t%s%s:%5.0f\t",r1,r2,LOD[j1][k1],s1,s2,LOD[j2][k2]);
                    for ( DiploidGenotype g : DiploidGenotype.values() ) {
                        out.printf("%s %5.0f\t",g.toString(),Scores.get(g.toString()));
                    }
                    out.printf("\n");
                }
            }
        }
        
        return context.getReads().size();
    }

    private int GetPhaseScore(int alleleindex){
        int score = 0;
        ArrayList<Integer> SNPsInRead = new ArrayList<Integer>();
        ArrayList<Integer> readindex = new ArrayList<Integer>();
        Hashtable indexer = new Hashtable();
        indexer.put('A', (Integer) 0);
        indexer.put('C', (Integer) 1);
        indexer.put('G', (Integer) 2);
        indexer.put('T', (Integer) 3);
        indexer.put('D', (Integer) 4); // D for deletion

        //Get HLA allele sequence and position given index
        String allele = HLAreads.get(alleleindex);
        int allelestart = HLAstartpos[alleleindex], allelestop = HLAstoppos[alleleindex];

        //Finds SNPs in allele
        for (int i = allelestart; i <= allelestop; i++){
            if (SNPs.containsKey(String.valueOf(i))){
                //Stores matrix index
                SNPsInRead.add((Integer) SNPs.get(String.valueOf(i)));
                //stores position along read
                readindex.add((Integer) i - allelestart);
            }
        }

        //sum score for every pair of SNPs in the allele
        for (int i = 0; i < SNPsInRead.size(); i++){
            for (int j = i+1; j < SNPsInRead.size(); j ++){
                char c1 = allele.charAt((int) readindex.get(i));
                char c2 = allele.charAt((int) readindex.get(j));
                if (indexer.get(c1) != null && indexer.get(c2) != null){
                    int a = (SNPsInRead.get(i)-1)*5 + (Integer) indexer.get(c1);
                    int b = (SNPsInRead.get(j)-1)*5 + (Integer) indexer.get(c2);
                    score += SNPcorrelation[a][b];
                }
            }
        }
        return score;
    }

    private void UpdateCorrelation(SAMRecord read, boolean PRINT){
        //Updates correlation table with SNPs from specific read (for phasing)
        String s = CigarFormatted(read.getCigarString(), read.getReadString());
        ArrayList<Integer> SNPsInRead = new ArrayList<Integer>();
        ArrayList<Integer> readindex = new ArrayList<Integer>();
        Hashtable indexer = new Hashtable();

        indexer.put('A', (Integer) 0);
        indexer.put('C', (Integer) 1);
        indexer.put('G', (Integer) 2);
        indexer.put('T', (Integer) 3);
        indexer.put('D', (Integer) 4); // D for deletion


        int readstart = read.getAlignmentStart();

        //Find all SNPs in read
        for (int i = read.getAlignmentStart(); i <= read.getAlignmentEnd(); i++){
            if (SNPs.containsKey(String.valueOf(i))){
                //Stores matrix index
                SNPsInRead.add((Integer) SNPs.get(String.valueOf(i)));

                //stores position along read
                readindex.add((Integer) i - readstart);
            }
        }

        //Update correlation table; for each combination of SNP positions
        for (int i = 0; i < SNPsInRead.size(); i++){
            for (int j = i+1; j < SNPsInRead.size(); j ++){
                char c1 = s.charAt((int) readindex.get(i));
                char c2 = s.charAt((int) readindex.get(j));

                if (indexer.get(c1) != null && indexer.get(c2) != null){

                    int a = (SNPsInRead.get(i)-1)*5 + (Integer) indexer.get(c1);
                    int b = (SNPsInRead.get(j)-1)*5 + (Integer) indexer.get(c2);
                    if (PRINT){
                        out.printf("ReadIndex[%s,%s] of %s\t", readindex.get(i),readindex.get(j),read.getAlignmentEnd()-readstart);
                        out.printf("SNP#[%s,%s] of %s\tPOS:%s[%s]\t%s[%s]\t", SNPsInRead.get(i),SNPsInRead.get(j), SNPs.size(), readindex.get(i)+readstart,c1, readindex.get(j)+readstart,c2);
                        out.printf("MatrixIndex[%s,%s] of %s",a,b,SNPcorrelation.length);
                        out.printf("\tc1=%s,c2=%s",c1,c2);
                        out.printf("\ta=%s,b=%s",a,b);
                        out.printf("\tSNPcorrelation[a][b]=%s\n",SNPcorrelation[a][b]);
                    }
                    SNPcorrelation[a][b]+=1;
                }
            }
        }
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
            out.print("\nDone calculating likelihoods\n");

            ArrayList<Integer> TopAlleles = new ArrayList<Integer>();

            double maxA = 0;  int i_maxA =0;   int j_maxA = 0; int maxAphase = 0; double maxAfreq = 0.0;
            double maxA2 = 0; int i_maxA_2 =0; int j_maxA_2 = 0;

            double maxB = 0;  int i_maxB =0;   int j_maxB = 0; int maxBphase = 0; double maxBfreq = 0.0;
            double maxB2 = 0; int i_maxB_2 =0; int j_maxB_2 = 0;

            double maxC = 0;  int i_maxC =0;   int j_maxC = 0; int maxCphase = 0; double maxCfreq = 0.0;
            double maxC2 = 0; int i_maxC_2 =0; int j_maxC_2 = 0;

            
            out.print("Finding allele pair with highest likelihood\n");
            //Find max likelihood scores for each HLA gene.
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = i; j < numHLAlleles; j++){
                    //Print likelihoods for all alleles
                    if (DEBUG){
                        out.printf("%s\t%s\t%5.0f\n",HLAnames.get(i),HLAnames.get(j),LOD[i][j]);
                    }
                    if (HLAnames.get(i).indexOf("HLA_A") > -1 && HLAnames.get(j).indexOf("HLA_A") > -1){
                        if (LOD[i][j] > maxA){
                            maxA2 = maxA; i_maxA_2 = i_maxA; j_maxA_2 = j_maxA;
                            maxA = LOD[i][j]; i_maxA = i; j_maxA = j;
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_B") > -1 && HLAnames.get(j).indexOf("HLA_B") > -1){
                        if (LOD[i][j] > maxB){
                            maxB2 = maxB; i_maxB_2 = i_maxB; j_maxB_2 = j_maxB;
                            maxB = LOD[i][j]; i_maxB = i; j_maxB = j;
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_C") > -1 && HLAnames.get(j).indexOf("HLA_C") > -1){
                        if (LOD[i][j] > maxC){
                            maxC2 = maxC; i_maxC_2 = i_maxC; j_maxC_2 = j_maxC;
                            maxC = LOD[i][j]; i_maxC = i; j_maxC = j;
                        }
                    }
		}
            }

            //Record alleles in highest likelihood combinations
            for (Integer i = 0; i < numHLAlleles; i++){
                for (Integer j = i; j < numHLAlleles; j++){
                    if (HLAnames.get(i).indexOf("HLA_A") > -1 && HLAnames.get(j).indexOf("HLA_A") > -1){
                        if (LOD[i][j] == maxA && maxA > 0){
                            if (!TopAlleles.contains(i)){TopAlleles.add(i);}
                            if (!TopAlleles.contains(j)){TopAlleles.add(j);}
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_B") > -1 && HLAnames.get(j).indexOf("HLA_B") > -1){
                        if (LOD[i][j] == maxB && maxB > 0){
                            if (!TopAlleles.contains(i)){TopAlleles.add(i);}
                            if (!TopAlleles.contains(j)){TopAlleles.add(j);}
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_C") > -1 && HLAnames.get(j).indexOf("HLA_C") > -1){
                        if (LOD[i][j] == maxC && maxC > 0){
                            if (!TopAlleles.contains(i)){TopAlleles.add(i);}
                            if (!TopAlleles.contains(j)){TopAlleles.add(j);}
                        }
                    }

                }
            }

            out.printf("\nCalculating SNP correlation matrix for %s SNPs\n",SNPcount);
            SNPcorrelation = new int[SNPs.size()*5][SNPs.size()*5];

            
            //Create correlation matrix and update correlation scores for all reads
            for (int i = 0; i < AllReads.size(); i++){
                if (i == 2045){
                    UpdateCorrelation(AllReads.get(i), false);
                }else{
                    UpdateCorrelation(AllReads.get(i), false);
                }
                //out.printf("[%s,%s]\n", ((Integer) SNPs.get("31431982")) * 4,((Integer) SNPs.get("31432003")) * 4 );
                //out.printf("%s\t[%s,%s]\t31431982[A] 31432003[A]\t%s\n",i,((Integer) SNPs.get("31431982")) * 4,((Integer) SNPs.get("31432003")) * 4 , SNPcorrelation[((Integer) SNPs.get("31431982")) * 4][((Integer) SNPs.get("31432003")) * 4]);
            }

            Hashtable indexer = new Hashtable();
            indexer.put('A', (Integer) 0);
            indexer.put('C', (Integer) 1);
            indexer.put('G', (Integer) 2);
            indexer.put('T', (Integer) 3);
            indexer.put('D', (Integer) 4); // D for deletion
            char[] bases = {'A','C','G','T','D'};

            if ( false ){
                //prints entries in the correlation matrix that are > 0
                out.print("\n");
                for (int i = 0; i < SNPs.size(); i++){
                    int loc1 = SNPlocations.get(i);
                    for (int j = i ; j < SNPs.size(); j++){
                        int loc2 = SNPlocations.get(j);
                        for (char c1 : bases){
                            for (char c2 : bases){
                                int a = i*5 + (Integer) indexer.get(c1);
                                int b = j*5 + (Integer) indexer.get(c2);
                                if (SNPcorrelation[a][b] > 0){
                                    out.printf("[i,j]=[%s,%s]\t[a,b]=[%s,%s]\tPOS:%s[%s],%s[%s]\tCorr=%s\n",i,j,a,b,loc1,c1,loc2,c2,SNPcorrelation[a][b]);
                                }
                            }
                        }
                    }
                }
            }


            int k, readstart, readstop, allelestart, allelestop, pos_k;
            out.printf("Calculating phase scores for %s top alleles\n",TopAlleles.size());

            //Calculate Phase score for each allele
            int[] SinglePhaseScores = new int[TopAlleles.size()];
            for (int i = 0; i < TopAlleles.size(); i++){
                SinglePhaseScores[i] = GetPhaseScore(TopAlleles.get(i));
                if ( DEBUG ){
                    //Debugging: print list of alleles to be checked for phasing
                    out.printf("index=%s\t%s\tscore=%s\n",TopAlleles.get(i),HLAnames.get(TopAlleles.get(i)),SinglePhaseScores[i]);
                }
            }

            out.print("Calculating phasing score for pairs of alleles\n");
            //Calculate phasing score and population frequencies for pairs of alleles, and find pairs with the highest scores
            String alleleA, alleleB;
            Double freq1 = 0.0, freq2 = 0.0;
            for (Integer i = 0; i < numHLAlleles; i++){
                for (Integer j = i; j < numHLAlleles; j++){
                    if (HLAnames.get(i).indexOf("HLA_A") > -1 && HLAnames.get(j).indexOf("HLA_A") > -1){
                        if (LOD[i][j] == maxA && maxA > 0){
                            PhasingScores[i][j]= SinglePhaseScores[TopAlleles.indexOf(i)] + SinglePhaseScores[TopAlleles.indexOf(j)];
                            if (PhasingScores[i][j] > maxAphase){maxAphase = PhasingScores[i][j];}
                            alleleA=HLAnames.get(i).substring(4); if (AlleleFrequencies.containsKey(alleleA)){freq1 = Double.parseDouble((String) AlleleFrequencies.get(alleleA).toString());}else{freq1=0.0001;}
                            alleleB=HLAnames.get(j).substring(4); if (AlleleFrequencies.containsKey(alleleB)){freq2 = Double.parseDouble((String) AlleleFrequencies.get(alleleB).toString());}else{freq2=0.0001;}
                            SingleAlleleFrequencies[i]=freq1; SingleAlleleFrequencies[j]=freq2; CombinedAlleleFrequencies[i][j]=freq1*freq2;
                            if (CombinedAlleleFrequencies[i][j] > maxAfreq){maxAfreq = CombinedAlleleFrequencies[i][j];}
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_B") > -1 && HLAnames.get(j).indexOf("HLA_B") > -1){
                        if (LOD[i][j] == maxB && maxB > 0){
                            PhasingScores[i][j]= SinglePhaseScores[TopAlleles.indexOf(i)] + SinglePhaseScores[TopAlleles.indexOf(j)];
                            if (PhasingScores[i][j] > maxBphase){maxBphase = PhasingScores[i][j];}
                            alleleA=HLAnames.get(i).substring(4); if (AlleleFrequencies.containsKey(alleleA)){freq1 = Double.parseDouble((String) AlleleFrequencies.get(alleleA).toString());}else{freq1=0.0001;}
                            alleleB=HLAnames.get(j).substring(4); if (AlleleFrequencies.containsKey(alleleB)){freq2 = Double.parseDouble((String) AlleleFrequencies.get(alleleB).toString());}else{freq2=0.0001;}
                            SingleAlleleFrequencies[i]=freq1; SingleAlleleFrequencies[j]=freq2; CombinedAlleleFrequencies[i][j]=freq1*freq2;
                            if (freq1*freq2 > maxBfreq){maxBfreq = freq1*freq2;}
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_C") > -1 && HLAnames.get(j).indexOf("HLA_C") > -1){
                        if (LOD[i][j] == maxC && maxC > 0){
                            PhasingScores[i][j]= SinglePhaseScores[TopAlleles.indexOf(i)] + SinglePhaseScores[TopAlleles.indexOf(j)];
                            if (PhasingScores[i][j] > maxCphase){maxCphase = PhasingScores[i][j];}
                            alleleA=HLAnames.get(i).substring(4); if (AlleleFrequencies.containsKey(alleleA)){freq1 = Double.parseDouble((String) AlleleFrequencies.get(alleleA).toString());}else{freq1=0.0001;}
                            alleleB=HLAnames.get(j).substring(4); if (AlleleFrequencies.containsKey(alleleB)){freq2 = Double.parseDouble((String) AlleleFrequencies.get(alleleB).toString());}else{freq2=0.0001;}
                            SingleAlleleFrequencies[i]=freq1; SingleAlleleFrequencies[j]=freq2; CombinedAlleleFrequencies[i][j]=freq1*freq2;
                            if (freq1*freq2 > maxCfreq){maxCfreq = freq1*freq2;}
                        }
                    }
                }
            }

            //Print allele pairs with highest likelihood score and highest phasing score
            for (Integer i = 0; i < numHLAlleles; i++){
                for (Integer j = i; j < numHLAlleles; j++){
                    if (HLAnames.get(i).indexOf("HLA_A") > -1 && HLAnames.get(j).indexOf("HLA_A") > -1){
                        if (LOD[i][j] == maxA && maxA > 0 && PhasingScores[i][j] == maxAphase && CombinedAlleleFrequencies[i][j] == maxAfreq){
                            out.printf("%s\t%s\tLOD=%5.0f\tCONF=%5.0f\tPhase=%s\tfreq1=%s\tfreq2=%s\tCombined_Freq=%.6f\tBEST\n",HLAnames.get(i),HLAnames.get(j),maxA,maxA-maxA2,PhasingScores[i][j],SingleAlleleFrequencies[i],SingleAlleleFrequencies[j],CombinedAlleleFrequencies[i][j]);
                        }else if(LOD[i][j] == maxA && maxA > 0 && PhasingScores[i][j] == maxAphase){
                            out.printf("%s\t%s\tLOD=%5.0f\tCONF=%5.0f\tPhase=%s\tfreq1=%s\tfreq2=%s\tCombined_Freq=%.6f\n",HLAnames.get(i),HLAnames.get(j),maxA,maxA-maxA2,PhasingScores[i][j],SingleAlleleFrequencies[i],SingleAlleleFrequencies[j],CombinedAlleleFrequencies[i][j]);
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_B") > -1 && HLAnames.get(j).indexOf("HLA_B") > -1){
                        if (LOD[i][j] == maxB && maxB > 0 && PhasingScores[i][j] == maxBphase && CombinedAlleleFrequencies[i][j] == maxBfreq){
                            out.printf("%s\t%s\tLOD=%5.0f\tCONF=%5.0f\tPhase=%s\tfreq1=%s\tfreq2=%s\tCombined_Freq=%.6f\tBEST\n",HLAnames.get(i),HLAnames.get(j),maxB,maxB-maxB2,PhasingScores[i][j],SingleAlleleFrequencies[i],SingleAlleleFrequencies[j],CombinedAlleleFrequencies[i][j]);
                        }else if(LOD[i][j] == maxB && maxB > 0 && PhasingScores[i][j] == maxBphase){
                            out.printf("%s\t%s\tLOD=%5.0f\tCONF=%5.0f\tPhase=%s\tfreq1=%s\tfreq2=%s\tCombined_Freq=%.6f\n",HLAnames.get(i),HLAnames.get(j),maxB,maxB-maxB2,PhasingScores[i][j],SingleAlleleFrequencies[i],SingleAlleleFrequencies[j],CombinedAlleleFrequencies[i][j]);
                        }
                    } else if (HLAnames.get(i).indexOf("HLA_C") > -1 && HLAnames.get(j).indexOf("HLA_C") > -1){
                        if (LOD[i][j] == maxC && maxC > 0 && PhasingScores[i][j] == maxCphase && CombinedAlleleFrequencies[i][j] == maxCfreq){
                            out.printf("%s\t%s\tLOD=%5.0f\tCONF=%5.0f\tPhase=%s\tfreq1=%s\tfreq2=%s\tCombined_Freq=%.6f\tBEST\n",HLAnames.get(i),HLAnames.get(j),maxC,maxC-maxC2,PhasingScores[i][j],SingleAlleleFrequencies[i],SingleAlleleFrequencies[j],CombinedAlleleFrequencies[i][j]);
                        }else if (LOD[i][j] == maxC && maxC > 0 && PhasingScores[i][j] == maxCphase){
                            out.printf("%s\t%s\tLOD=%5.0f\tCONF=%5.0f\tPhase=%s\tfreq1=%s\tfreq2=%s\tCombined_Freq=%.6f\n",HLAnames.get(i),HLAnames.get(j),maxC,maxC-maxC2,PhasingScores[i][j],SingleAlleleFrequencies[i],SingleAlleleFrequencies[j],CombinedAlleleFrequencies[i][j]);
                        }
                    }
                }
            }

            //2nd Highest likelihoods
            for (int i = 0; i < numHLAlleles; i++){
                for (int j = i; j < numHLAlleles; j++){
                    if (LOD[i][j] == maxA2){
                        //out.printf("2nd Highest likelihood: %5.0f in %s and %s; i=%s, j=%s\n",maxA2,HLAnames.get(i),HLAnames.get(j),i,j);
                    }
                }
            }

            PrintedOutput = true;
            //out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n",((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
        }
    }
}
