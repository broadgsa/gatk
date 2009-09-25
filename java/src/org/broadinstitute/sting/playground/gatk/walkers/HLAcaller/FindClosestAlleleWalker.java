/*
 * FindClosestAlleleWalker finds the most similar HLA allele per read
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.playground.gatk.walkers.HLAcaller.ReadCigarFormatter;
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
    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.sam";
    String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";
    
    boolean DatabaseLoaded = false;
    boolean DEBUG = false;

    ArrayList<String> HLAreads = new ArrayList<String>();
    ArrayList<String> HLAcigars = new ArrayList<String>();
    ArrayList<String> HLAnames = new ArrayList<String>();
    ArrayList<String> HLApositions = new ArrayList<String>();
    double[] SingleAlleleFrequencies;

    int numHLAlleles = 0;
    int[] HLAstartpos;
    int[] HLAstoppos;
    int minstartpos = 0;
    int maxstoppos = 0;

    int HLA_A_start = 30018310;
    int HLA_A_end = 30021211;
    
    ArrayList<String> PolymorphicSites = new ArrayList<String>();

    Hashtable AlleleFrequencies = new Hashtable();
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;
    ReadCigarFormatter formatter = new ReadCigarFormatter();

    public Integer reduceInit() { 
    if (!DatabaseLoaded){
            try{
                out.printf("Reading HLA database ...\n");
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
                        HLAreads.add(formatter.FormatRead(s[5],s[9]));
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
                    if (minstartpos == 0){minstartpos = HLAstartpos[i];}
                    minstartpos = Math.min(minstartpos, HLAstartpos[i]);
                    maxstoppos = Math.max(maxstoppos, HLAstoppos[i]);
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
                    AlleleFrequencies.put(s[0].substring(0, 6), s[1]);
                    count++;
                }
                in.close();
                out.printf("Done! Read %s alleles\n",count);
            }catch (Exception e){//Catch exception if any
              System.err.println("Error: " + e.getMessage());
            }

            char c;
            DatabaseLoaded = true;
            //Find polymorphic sites in dictionary
            for (int pos = minstartpos; pos <= maxstoppos; pos++){
                c = '0';
                for (int i = 0; i < HLAreads.size(); i++){
                    if (pos >= HLAstartpos[i] && pos <= HLAstoppos[i]){
                        if (c == '0'){c = HLAreads.get(i).charAt(pos-HLAstartpos[i]);}
                        if (HLAreads.get(i).charAt(pos-HLAstartpos[i]) != c){
                            PolymorphicSites.add(String.valueOf(pos));
                            break;
                        }
                    }
                }
            }
            out.printf("%s polymorphic sites found in HLA dictionary\n",PolymorphicSites.size());

            out.printf("Comparing reads to database ...\n");

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
    }


    public Integer map(char[] ref, SAMRecord read) {
        int readstart = read.getAlignmentStart();
        int readstop = read.getAlignmentEnd();
        double[] nummatched = new double[HLAreads.size()];
        double[] concordance = new double[HLAreads.size()];
        double[] numcompared = new double[HLAreads.size()];
        double maxConcordance = 0;
        String s1 = formatter.FormatRead(read.getCigarString(), read.getReadString());
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
                        if (PolymorphicSites.contains(String.valueOf(j))){
                            numcompared[i]++;
                            if (c1 == c2){
                                nummatched[i]++;
                            }
                        }else if (c1 != c2){
                            numcompared[i]++;
                        }
                    }
                }
            }
            concordance[i]=nummatched[i]/numcompared[i];
            
            if (concordance[i] > maxConcordance){maxConcordance = concordance[i];}
        }
        double freq, freq2, maxFreq = 0.0;
        for (int i = 0; i < HLAreads.size(); i++){
            if (concordance[i] == maxConcordance && maxConcordance > 0){
                if (HLAnames.get(i).length() >= 10){
                    name=HLAnames.get(i).substring(4,10);
                }else{
                    name=HLAnames.get(i).substring(4);
                }
                if (AlleleFrequencies.containsKey(name)){
                    freq = Double.parseDouble((String) AlleleFrequencies.get(name).toString());
                    if (DEBUG){
                        out.printf("%s\t%s\t%s\t%s\t%s\t%.4f\n",read.getReadName(),name,nummatched[i],numcompared[i],concordance[i],freq);
                    }
                }else{
                    freq=0.0001;
                }
                if (freq > maxFreq){maxFreq = freq;}
            }
        }

        if (read.getReadName().length() >= 10){
            name=read.getReadName().substring(4,10);
        }else{
            name=read.getReadName().substring(4);
        }
        
        if (AlleleFrequencies.containsKey(name)){
            freq = Double.parseDouble((String) AlleleFrequencies.get(name).toString());
        }else{
            freq=0.0001;
        }

        out.printf("%s\t%.4f", read.getReadName(),freq);
        for (int i = 0; i < HLAreads.size(); i++){
            if (concordance[i] == maxConcordance && maxConcordance > 0){
                if (HLAnames.get(i).length() >= 10){
                    name=HLAnames.get(i).substring(4,10);
                }else{
                    name=HLAnames.get(i).substring(4);
                }
                if (AlleleFrequencies.containsKey(name)){
                    freq = Double.parseDouble((String) AlleleFrequencies.get(name).toString());
                }else{
                    freq=0.0001;
                }
                if (freq == maxFreq){

                    out.printf("\t%s\t%.4f\t%.3f\t%.0f\t%.0f",HLAnames.get(i),freq,concordance[i],numcompared[i],numcompared[i]-nummatched[i]);
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

