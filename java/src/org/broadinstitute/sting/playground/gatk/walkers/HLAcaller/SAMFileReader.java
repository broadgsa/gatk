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
public class SAMFileReader {
    ArrayList<String> ReadStrings = new ArrayList<String>();
    ArrayList<String> CigarStrings = new ArrayList<String>();
    ArrayList<String> ReadNames = new ArrayList<String>();
    ArrayList<Integer> ReadStartPositions = new ArrayList<Integer>();
    ArrayList<Integer> ReadStopPositions = new ArrayList<Integer>();
    int minstartpos;
    int maxstoppos;

    CigarParser formatter = new CigarParser();

    public String[] GetReads(){
        return ReadStrings.toArray(new String[ReadStrings.size()]);
    }

    public String[] GetReadNames(){
        return ReadNames.toArray(new String[ReadNames.size()]);
    }

    public String[] GetCigarStrings(){
        return CigarStrings.toArray(new String[CigarStrings.size()]);
    }

    public Integer[] GetStartPositions(){
        return ReadStartPositions.toArray(new Integer[ReadStartPositions.size()]);
    }

    public Integer[] GetStopPositions(){
        return ReadStopPositions.toArray(new Integer[ReadStopPositions.size()]);
    }

    public Integer GetMinStartPos(){
        return minstartpos;
    }

    public Integer GetMaxStopPos(){
        return maxstoppos;
    }

    public int GetReadIndex(String readname){
        if (ReadNames.contains(readname)){
            return ReadNames.indexOf(readname);
        }else{
            return -1;
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
                if (s.length>=10){
                    //Parse the reads with cigar parser
                    String read = formatter.FormatRead(s[5],s[9]);
                    ReadStrings.add(read);
                    CigarStrings.add(s[5]);
                    ReadNames.add(s[0]);
                    ReadStartPositions.add(Integer.valueOf(s[3]));
                    ReadStopPositions.add(Integer.valueOf(s[3]) + read.length() - 1);
                    if (i == 0){
                        minstartpos = Integer.valueOf(s[3]);
                        maxstoppos = Integer.valueOf(Integer.valueOf(s[3]) + read.length() - 1);
                    }
                    minstartpos = Math.min(minstartpos, Integer.valueOf(s[3]));
                    maxstoppos = Math.max(maxstoppos, Integer.valueOf(Integer.valueOf(s[3]) + read.length() - 1));
                    i++;
                }
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("SAMFileReader Error: " + e.getMessage());
        }
    }
}

