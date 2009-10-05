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
public class TextFileReader {
    ArrayList<String> lines = new ArrayList<String>();
    int numLines = 0;

    public String[] GetLines(){
        return lines.toArray(new String[lines.size()]);
    }

    public int GetNumLines(){
        return numLines;
    }

    public void ReadFile(String filename){
        try{
            FileInputStream fstream = new FileInputStream(filename);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine;
            //Read File Line By Line
            while ((strLine = br.readLine()) != null)   {
                lines.add(strLine);
                numLines++;
            }
            in.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("TextFileReader Error: " + e.getMessage());
        }
    }
}

