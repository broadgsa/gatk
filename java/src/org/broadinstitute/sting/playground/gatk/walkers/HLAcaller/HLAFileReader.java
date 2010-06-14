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
public class HLAFileReader {
    ArrayList<String> Sequences = new ArrayList<String>();
    ArrayList<String> Names = new ArrayList<String>();
    ArrayList<Integer> StartPositions = new ArrayList<Integer>();
    ArrayList<Integer> StopPositions = new ArrayList<Integer>();
    int minstartpos;
    int maxstoppos;

    CigarParser formatter = new CigarParser();

    public String[] GetNames(){
        return Names.toArray(new String[Names.size()]);
    }

    public String[] GetSequences(){
        return Sequences.toArray(new String[Sequences.size()]);
    }

    public Integer[] GetStartPositions(){
        return StartPositions.toArray(new Integer[StartPositions.size()]);
    }

    public Integer[] GetStopPositions(){
        return StopPositions.toArray(new Integer[StopPositions.size()]);
    }


    public Integer GetMinStartPos(){
        return minstartpos;
    }

    public Integer GetMaxStopPos(){
        return maxstoppos;
    }
    
    public int GetIndex(String readname){
        if (Names.contains(readname)){
            return Names.indexOf(readname);
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
            while ((strLine = br.readLine()) != null)   {
                s = strLine.split("\\t");
                Sequences.add(s[3]);
                Names.add(s[0]);
                StartPositions.add(Integer.valueOf(s[1]));
                StopPositions.add(Integer.valueOf(s[2]));
                minstartpos = Math.min(minstartpos, Integer.valueOf(s[1]));
                maxstoppos = Math.max(maxstoppos, Integer.valueOf(s[2]));
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("HLAFileReader Error: " + e.getMessage());
        }
    }
}

