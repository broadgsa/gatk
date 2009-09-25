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
public class ImputeAllelesWalker extends ReadWalker<Integer, Integer> {
    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.sam";
    String ClosestAllelesFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.closest";

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
    int HLA_B_start = 31430239;
    int HLA_B_end = 31432914;
    int HLA_C_start = 31344925;
    int HLA_C_end = 31347827;


    ArrayList<String> PolymorphicSites = new ArrayList<String>();

    Hashtable ClosestAllele = new Hashtable();
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
                out.printf("Reading closest allele file ...");
                FileInputStream fstream = new FileInputStream(ClosestAllelesFile);
                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine; String [] s = null;
                //Read File Line By Line
                int count = 0;
                while ((strLine = br.readLine()) != null)   {
                    s = strLine.split("\\t");
                    ClosestAllele.put(s[0], s[2]);
//                    out.printf("loading: %s\t%s\n",s[0],s[2]);
                    count++;
                }
                in.close();
                out.printf("Done! Read %s alleles\n",count);
            }catch (Exception e){//Catch exception if any
              System.err.println("Error: " + e.getMessage());
            }

            char c;
            DatabaseLoaded = true;
            
            out.printf("Imputing alleles ...\n");

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
    }


    public Integer map(char[] ref, SAMRecord read) {
        int readstart = read.getAlignmentStart();
        int readstop = read.getAlignmentEnd();
        int startimputation = 0, stopimputation = 0;

        String s1 = formatter.FormatRead(read.getCigarString(), read.getReadString());
        char c;
        String readstring = "", name = "", cigar = "", qualitystring = "";
        int numM = 0, numI = 0, numD = 0;

        name = read.getReadName();
        String matchedAllele = (String) ClosestAllele.get(name);

        //out.printf("%s\t%s\n",name,matchedAllele);
        int index = HLAnames.indexOf(matchedAllele);
        String matchedRead = HLAreads.get(index);

        if (name.indexOf("HLA_A") > -1){
            startimputation = HLA_A_start;
            stopimputation = HLA_A_end;
        } else if (name.indexOf("HLA_B") > -1){
            startimputation = HLA_B_start;
            stopimputation = HLA_B_end;
        } else if (name.indexOf("HLA_C") > -1){
            startimputation = HLA_C_start;
            stopimputation = HLA_C_end;
        }

        for (int i = startimputation; i <= stopimputation; i++){
            //if position is within read
            if (i >= readstart && i <= readstop){
                c = s1.charAt(i-readstart);
                //if position is not missing
                if (c != 'D'){
                    readstring = readstring + c;
                    qualitystring = qualitystring + 'I';
                    numM++;
                    if (numD > 0){
                        cigar = cigar + String.valueOf(numD) + "D";
                        numD = 0;
                    } else if (numI > 0){
                        cigar = cigar + String.valueOf(numI) + "I";
                        numI = 0;
                    }
                //if position is missing, get base from matched allele
                }else{
                    c = matchedRead.charAt(i-HLAstartpos[index]);
                    //if matched allele is also missing / deleted at position
                    if (c == 'D'){
                        numD++;
                        if (numM > 0){
                            cigar = cigar + String.valueOf(numM) + "M";
                            numM = 0;
                        }
                    //if matched allele is not missing / deleted at position
                    }else{
                        readstring = readstring + c;
                        qualitystring = qualitystring + 'I';
                        numM++;
                        if (numD > 0){
                            cigar = cigar + String.valueOf(numD) + "D";
                            numD = 0;
                        } else if (numI > 0){
                            cigar = cigar + String.valueOf(numI) + "I";
                            numI = 0;
                        }
                    }
                }
            //if position is outside of range of read, look at matched allele
            }else{
                //if within range of matched allele
                if (i >= HLAstartpos[index] && i <= HLAstoppos[index]){
                    c = matchedRead.charAt(i-HLAstartpos[index]);
                    //if matched allele is also missing / deleted at position
                    if (c == 'D'){
                        numD++;
                        if (numM > 0){
                            cigar = cigar + String.valueOf(numM) + "M";
                            numM = 0;
                        }
                    //if matched allele is not missing / deleted at position
                    }else{
                        readstring = readstring + c;
                        qualitystring = qualitystring + 'I';
                        numM++;
                        if (numD > 0){
                            cigar = cigar + String.valueOf(numD) + "D";
                            numD = 0;
                        } else if (numI > 0){
                            cigar = cigar + String.valueOf(numI) + "I";
                            numI = 0;
                        }
                    }
                }else{
                    numD++;
                    if (numM > 0){
                        cigar = cigar + String.valueOf(numM) + "M";
                        numM = 0;
                    }
                }
            }
        }

        if (numM > 0){
            cigar = cigar + String.valueOf(numM) + "M";
        }else if(numD > 0){
            cigar = cigar + String.valueOf(numD) + "D";
        }else if(numI > 0){
            cigar = cigar + String.valueOf(numI) + "I";
        }
        
        out.printf("%s\t0\t6\t%s\t99\t%s\t*\t0\t0\t%s\t%s\n",name,startimputation,cigar,readstring,qualitystring);
        
        
        return 1;
    }




    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}

