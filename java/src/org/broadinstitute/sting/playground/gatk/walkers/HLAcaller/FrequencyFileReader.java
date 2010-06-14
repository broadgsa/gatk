package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import java.io.*;
import java.util.Hashtable;
/**
 * File reader used by other Walkers to read HLA allele frequencies.
 * @author shermanjia
 */
public class FrequencyFileReader {
    Hashtable MaxFrequencies = new Hashtable();
    Hashtable CommonAlleles = new Hashtable();
    Hashtable [] AlleleFrequencies = null;
    String [] Populations = null;

    public Hashtable [] GetAlleleFrequencies(){
        //return allele frequencies for all populations
        return AlleleFrequencies;
    }

    public Hashtable GetCommonAlleles(){
        //return list of common alleles
        return CommonAlleles;
    }

    public Hashtable GetMaxFrequencies(){
        //return list of common alleles
        return MaxFrequencies;
    }

    public String[] GetPopulations(){
        //Return name of populations
        return Populations;
    }

    public void ReadFile(String filename, String ethnicity){
        try{
            int linenum = 0;
            FileInputStream fstream = new FileInputStream(filename);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine; String [] s = null;
            //Read File Line By Line
            while ((strLine = br.readLine()) != null)   {
                linenum++;
                s = strLine.split("\\t");
                if (linenum == 1){
                    //Determine number of populations, create a hash table for each population
                    AlleleFrequencies = new Hashtable[s.length-1];
                    Populations = new String[s.length-1];
                    for (int i = 1; i < s.length; i++){
                        Populations[i-1]=s[i];
                        AlleleFrequencies[i-1] = new Hashtable();
                    }
                }else{
                    //assign allele frequencies for each population
                    for (int i = 1; i < s.length; i++){
                        if (Double.valueOf(s[i]) > 0.0001){
                            CommonAlleles.put(s[0], s[0]);
                        }
                        AlleleFrequencies[i-1].put(s[0],s[i]);
                        if (!MaxFrequencies.containsKey(s[0])){
                            MaxFrequencies.put(s[0], s[i]);
                        }else if (Double.valueOf(MaxFrequencies.get(s[0]).toString()) < Double.valueOf(s[i])){
                            MaxFrequencies.put(s[0], s[i]);
                        }
                    }
                }
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Exception in FrequencyFileReader (" + e.getMessage() + ").");
        }
    }
}

