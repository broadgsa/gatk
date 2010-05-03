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
    Hashtable Concordance = new Hashtable();
    Hashtable NumMatches = new Hashtable();
    Hashtable NumMismatches = new Hashtable();

    public ArrayList<String> GetReadsToDiscard(){
        return ReadsToDiscard;
    }

    public String[] GetReadsToDiscardArray(){
        return ReadsToDiscard.toArray(new String[ReadsToDiscard.size()]);
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
            String strLine; String [] s = null;
            //Read File Line By Line
            int i = 0;
            while ((strLine = br.readLine()) != null)   {
                s = strLine.split("\\t");
                if (s.length >= 6){
                    Double matchFraction = Double.valueOf(s[4]);
                    int numMismatches = Integer.valueOf(s[6]);

                    Concordance.put(s[0],matchFraction);
                    NumMatches.put(s[0], s[5]);
                    NumMismatches.put(s[0], numMismatches);
                    if ((matchFraction < 0.9 && numMismatches > 3) || (numMismatches > minAllowedMismatches)){
                        ReadsToDiscard.add(s[0]);
                    }
                }
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("SimilarityFile Error: " + e.getMessage());
        }
    }
}

