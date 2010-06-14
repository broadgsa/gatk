/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import java.io.*;
import java.util.ArrayList;
import java.util.Hashtable;
/**
 *
 * @author shermanjia
 */
public class SimilarityFileReader {
    ArrayList<String> ReadsToDiscard = new ArrayList<String>();
    ArrayList<String> AllelesToSearch = new ArrayList<String>();
    Hashtable AlleleCount = new Hashtable();
    Hashtable LocusCount = new Hashtable();
    Hashtable Concordance = new Hashtable();
    Hashtable NumMatches = new Hashtable();
    Hashtable NumMismatches = new Hashtable();

    public ArrayList<String> GetReadsToDiscard(){
        return ReadsToDiscard;
    }

    public ArrayList<String> GetAllelesToSearch(){
        return AllelesToSearch;
    }

    public String[] GetReadsToDiscardArray(){
        return ReadsToDiscard.toArray(new String[ReadsToDiscard.size()]);
    }

    public Hashtable GetAlleleCount(){
        return AlleleCount;
    }

    public Hashtable GetLocusCount(){
        return LocusCount;
    }

    public Hashtable GetConcordance(){
        return Concordance;
    }

    public Hashtable GetNumMatches(){
        return NumMatches;
    }

    public Hashtable GetNumMismatches(){
        return NumMismatches;
    }

    public void ReadFile(String filename, int minAllowedMismatches){
        try{
            FileInputStream fstream = new FileInputStream(filename);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine; String [] s = null, alleles = null, a = null; String allele;
            //Read File Line By Line
            int i = 0;
            while ((strLine = br.readLine()) != null)   {
                s = strLine.split("\\t");
                if (s.length >= 6){
                    Double matchFraction = Double.valueOf(s[3]);
                    int numMismatches = Integer.valueOf(s[5]);
                    int numMatches = Integer.valueOf(s[4]);
                    Concordance.put(s[0],matchFraction);
                    NumMatches.put(s[0], s[4]);
                    NumMismatches.put(s[0], numMismatches);
                    if ((matchFraction < 0.8 && numMismatches > 3) || (numMismatches > minAllowedMismatches) || numMatches < 10){
                        ReadsToDiscard.add(s[0]);
                    }else{
                        Hashtable fourDigitAlleles = new Hashtable();
                        alleles = s[6].split("\\,");
                        if (alleles.length > 0){
                            a = alleles[0].split("\\_");
                            s = a[1].split("\\*");
                            if (!LocusCount.containsKey(s[0])){
                                LocusCount.put(s[0], 1);
                            }else{
                                LocusCount.put(s[0], (Integer) LocusCount.get(s[0]) + 1);
                            }
                        }
                        for (int j = 0; j < alleles.length; j++){
                            a = alleles[j].split("\\_");
                            s = a[1].split("\\*");
                            allele = s[0] + "*" + s[1].substring(0,4);
                            
                            if (!fourDigitAlleles.containsKey(allele)){
                                fourDigitAlleles.put(allele, allele);
                                if (!AlleleCount.containsKey(allele)){
                                    AlleleCount.put(allele, 1);
                                }else{
                                    AlleleCount.put(allele, (Integer) AlleleCount.get(allele) + 1);
                                }

                                if ((Integer) AlleleCount.get(allele) > 1 && !AllelesToSearch.contains(allele)){
                                    AllelesToSearch.add(allele);
                                }
                            }
                        }
                    }
                }
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            //System.err.println("SimilarityFile Error: " + e.getMessage());
        }
    }
}

