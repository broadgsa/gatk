/*
 * FindClosestAlleleWalker finds the most similar HLA allele per read
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.*;

import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.lang.Math;
/**
 *
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class FindClosestAlleleWalker extends ReadWalker<Integer, Integer> {
    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.4digitUnique.sam";
    String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";
    
    boolean DatabaseLoaded = false;
    boolean DEBUG = true;

    ArrayList<String> HLAreads = new ArrayList<String>();
    ArrayList<String> HLAcigars = new ArrayList<String>();
    ArrayList<String> HLAnames = new ArrayList<String>();
    ArrayList<String> HLApositions = new ArrayList<String>();
    double[] SingleAlleleFrequencies;

    int numHLAlleles = 0;
    int[] HLAstartpos;
    int[] HLAstoppos;

    Hashtable AlleleFrequencies = new Hashtable();
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;


    public Integer reduceInit() { 
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


                for (i = 0; i < n; i++){
                    //Find start and stop positions for each allele
                    HLAstartpos[i]=Integer.parseInt(HLApositions.get(i));
                    HLAstoppos[i]=HLAstartpos[i]+HLAreads.get(i).length()-1;
                    SingleAlleleFrequencies[i]=0.0;
                    //Initialize matrix of probabilities / likelihoods

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

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
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
            } else if (c == 'H' || c == 'S'){
                //(H = Headers or S = Soft clipped removed here)***

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


    public Integer map(char[] ref, SAMRecord read) {
        int readstart = read.getAlignmentStart();
        int readstop = read.getAlignmentEnd();
        double[] nummatched = new double[HLAreads.size()];
        double[] concordance = new double[HLAreads.size()];
        double[] numcompared = new double[HLAreads.size()];
        double maxConcordance = 0;
        String s1 = CigarFormatted(read.getCigarString(), read.getReadString());
        char c1, c2;
        String s2 = "", name = "";

        //For every allele that overlaps with current read
        for (int i = 0; i < HLAreads.size(); i++){
            nummatched[i] = 0;
            //Get concordance between read and specific allele
            if (DEBUG){
                //out.printf("%s\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(), HLAnames.get(i),readstart,readstop,HLAstartpos[i],HLAstoppos[i]);
            }
            if (readstart <= HLAstoppos[i] && readstop >= HLAstartpos[i]){
                s2 = HLAreads.get(i);
                for (int j = read.getAlignmentStart(); j <= read.getAlignmentEnd(); j++){
                    c1 = s1.charAt(j-readstart);
                    c2 = s2.charAt(j-HLAstartpos[i]);
                    if (DEBUG){
                        //out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(),HLAnames.get(i),j,j-readstart,j-HLAstartpos[i],c1,c2);
                    }
                    if (c1 != 'D'){
                        numcompared[i]++;
                        if (c1 == c2){
                            nummatched[i]++;
                        }
                    }
                }
            }
            concordance[i]=nummatched[i]/numcompared[i];
            if (concordance[i] > maxConcordance){maxConcordance = concordance[i];}
        }
        double freq, maxFreq = 0.0;
        for (int i = 0; i < HLAreads.size(); i++){
            if (concordance[i] == maxConcordance && maxConcordance > 0){
                name=HLAnames.get(i).substring(4);
                if (AlleleFrequencies.containsKey(name)){
                    freq = Double.parseDouble((String) AlleleFrequencies.get(name).toString());
                }else{
                    freq=0.0001;
                }
                if (freq > maxFreq){maxFreq = freq;}
            }
        }

        out.printf("%s", read.getReadName());
        for (int i = 0; i < HLAreads.size(); i++){
            if (concordance[i] == maxConcordance && maxConcordance > 0){
                name=HLAnames.get(i).substring(4);
                if (AlleleFrequencies.containsKey(name)){
                    freq = Double.parseDouble((String) AlleleFrequencies.get(name).toString());
                }else{
                    freq=0.0001;
                }
                if (freq == maxFreq){
                    out.printf("\t%s\t%.3f\t%.0f\t%.0f\t%.3f",HLAnames.get(i),concordance[i],numcompared[i],numcompared[i]-nummatched[i],freq);
                }
            }
        }
        out.print("\n");
        return 1;
    }


    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}

