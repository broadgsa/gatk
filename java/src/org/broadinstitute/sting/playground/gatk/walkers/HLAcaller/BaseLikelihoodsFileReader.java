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
public class BaseLikelihoodsFileReader {
    double[][] baseLikelihoods;
    int[] positions;
    ArrayList<Integer> polymorphicSites = new ArrayList<Integer>();

    public Integer[] GetPolymorphicSites(){
        return polymorphicSites.toArray(new Integer[polymorphicSites.size()]);
    }

    public double[][] GetBaseLikelihoods(){
        return baseLikelihoods;
    }

    public int[] GetPositions(){
        return positions;
    }

    public void ReadFile(String filename, boolean findPolymorphicSites){
        try{
            //System.out.printf("INFO  Reading base likelihoods file ... ");
            FileInputStream fstream = new FileInputStream(filename);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine; String [] s = null, pos = null;
            //Determine size of file
            int n = 0; char ref;
            while ((strLine = br.readLine()) != null)   {
                if (strLine.indexOf("INFO") == -1){n++;}
            }

            baseLikelihoods = new double[n][10];
            positions = new int[n];
            double[] localLikelihoods = new double[10];
            

            //System.out.printf("%s lines of data found ... ",n);
            in.close();
            //Read and store data

            fstream = new FileInputStream(filename);
            in = new DataInputStream(fstream);
            br = new BufferedReader(new InputStreamReader(in));
            n = 0;
            while ((strLine = br.readLine()) != null)   {
                if (strLine.indexOf("INFO") == -1){
                    s = strLine.split("\\t");
                    pos = s[0].split(":");
                    ref = s[1].charAt(0);
                    positions[n] = Integer.valueOf(pos[1]);
                    for (int i = 3; i <= 12; i++){
                        baseLikelihoods[n][i-3] = Double.valueOf(s[i]);
                        localLikelihoods[i-3] = baseLikelihoods[n][i-3];
                    }
                    if (findPolymorphicSites){
                        if (IsPolymorphicSite(localLikelihoods,ref)){
                            polymorphicSites.add(positions[n]);
                        }
                    }
                    n++;
                }
            }

        }catch (Exception e){//Catch exception if any
            System.err.println("BaseLikelihoodsFileReader Error: " + e.getMessage());
        }
    }

    private int IndexOf(char c){
        switch(c){
            case 'A': return 0;
            case 'C': return 4;
            case 'G': return 7;
            case 'T': return 9;
            default: return -1;
        }
    }

    private boolean IsPolymorphicSite(double[] likelihoods, char ref){
        boolean isPolymorphicSite = false;
        double homreflikelihood = likelihoods[IndexOf(ref)];
        for (int i = 0; i < 10; i++){
            if (likelihoods[i] > homreflikelihood){
                isPolymorphicSite = true;
            }
        }
        return isPolymorphicSite;
    }
}

