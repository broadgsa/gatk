/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import java.io.*;
import java.util.Hashtable;
/**
 *
 * @author shermanjia
 */
public class FrequencyFileReader {
    Hashtable AlleleFrequencies = new Hashtable();

    public Hashtable GetAlleleFrequencies(){
        return AlleleFrequencies;
    }

    public void ReadFile(String filename){
        try{
            FileInputStream fstream = new FileInputStream(filename);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine; String [] s = null;
            //Read File Line By Line
            while ((strLine = br.readLine()) != null)   {
                s = strLine.split("\\t");
                AlleleFrequencies.put(s[0], s[1]);
                //System.out.printf("Loaded: %s\t%s\n",s[0],AlleleFrequencies.get(s[0]).toString());
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("FrequencyFileReader Error: " + e.getMessage());
        }
    }
}

