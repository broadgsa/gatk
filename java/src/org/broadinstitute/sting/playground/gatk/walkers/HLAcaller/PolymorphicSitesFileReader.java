/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import java.io.*;
import java.util.ArrayList;
/**
 *
 * @author shermanjia
 */
public class PolymorphicSitesFileReader {
    ArrayList<Integer> PolymorphicSites = new ArrayList<Integer>();
    ArrayList<Integer> NonPolymorphicSites = new ArrayList<Integer>();


    public Integer[] GetPolymorphicSites(){
        return PolymorphicSites.toArray(new Integer[PolymorphicSites.size()]);
    }

    public Integer[] GetNonPolymorphicSites(){
        return NonPolymorphicSites.toArray(new Integer[NonPolymorphicSites.size()]);
    }

    public void AddSites(Integer [] sites){
        for (int i = 0; i < sites.length; i++){
            if (!PolymorphicSites.contains(sites[i])){
                PolymorphicSites.add(sites[i]);
            }
        }
    }

    public void ReadFile(String filename){
        try{
            FileInputStream fstream = new FileInputStream(filename);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine; String [] s = null;
            //Read File Line By Line
            int i = 0;
            while ((strLine = br.readLine()) != null)   {
                s = strLine.split("\\t");
                if (Double.valueOf(s[8]) > 0.1){
                    PolymorphicSites.add(Integer.valueOf(s[0]));
                }else{
                    NonPolymorphicSites.add(Integer.valueOf(s[0]));
                }
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("PolymorphicSitesFileReader Error: " + e.getMessage());
        }
    }
}

