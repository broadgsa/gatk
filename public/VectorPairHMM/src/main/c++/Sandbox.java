/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.vectorpairhmm;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.io.File;
import java.util.Scanner;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;

public class Sandbox {

    private long setupTime = 0;
    private long computeTime = 0;
    //Used to copy references to byteArrays to JNI from reads
    protected class JNIReadDataHolderClass {
        public byte[] readBases = null;
        public byte[] readQuals = null;
        public byte[] insertionGOP = null;
        public byte[] deletionGOP = null;
        public byte[] overallGCP = null;
    }

    //Used to copy references to byteArrays to JNI from haplotypes
    protected class JNIHaplotypeDataHolderClass {
        public byte[] haplotypeBases = null;
    }

    /**
     * Return 64-bit mask representing machine capabilities
     * Bit 0 is LSB, bit 63 MSB
     * Bit 0 represents sse4.2 availability
     * Bit 1 represents AVX availability
     */
    public native long jniGetMachineType();
    public static final long enableAll = 0xFFFFFFFFFFFFFFFFl;

    
    /**
     * Function to initialize the fields of JNIReadDataHolderClass and JNIHaplotypeDataHolderClass from JVM.
     * C++ codegets FieldIDs for these classes once and re-uses these IDs for the remainder of the program. Field IDs do not
     * change per JVM session
     * @param readDataHolderClass class type of JNIReadDataHolderClass
     * @param haplotypeDataHolderClass class type of JNIHaplotypeDataHolderClass
     * @param mask 64 bit integer identical to the one received from jniGetMachineType(). Users can disable usage of some hardware features by zeroing bits in the mask
     * */
    private native void jniInitializeClassFieldsAndMachineMask(Class<?> readDataHolderClass, Class<?> haplotypeDataHolderClass, long mask);

    private static Boolean isVectorLoglessPairHMMLibraryLoaded = false;
    //The constructor is called only once inside PairHMMLikelihoodCalculationEngine
    public Sandbox() {
        synchronized(isVectorLoglessPairHMMLibraryLoaded) {
            //Load the library and initialize the FieldIDs
            if(!isVectorLoglessPairHMMLibraryLoaded) {
                System.loadLibrary("VectorLoglessPairHMM");
                isVectorLoglessPairHMMLibraryLoaded = true;
                jniInitializeClassFieldsAndMachineMask(JNIReadDataHolderClass.class, JNIHaplotypeDataHolderClass.class, enableAll);        //need to do this only once
            }
        }
    }

    private native void jniInitializeHaplotypes(final int numHaplotypes,  JNIHaplotypeDataHolderClass[] haplotypeDataArray);
    private JNIHaplotypeDataHolderClass[] mHaplotypeDataArray = null;
    
    //Used to transfer data to JNI
    //Since the haplotypes are the same for all calls to computeLikelihoods within a region, transfer the haplotypes only once to the JNI per region
    public void initialize(final List<JNIHaplotypeDataHolderClass> haplotypes) {
        int numHaplotypes = haplotypes.size();
        mHaplotypeDataArray = new JNIHaplotypeDataHolderClass[numHaplotypes];
        int idx = 0;
        for(final JNIHaplotypeDataHolderClass currHaplotype : haplotypes)
        {
            mHaplotypeDataArray[idx] = new JNIHaplotypeDataHolderClass();
            mHaplotypeDataArray[idx].haplotypeBases = currHaplotype.haplotypeBases;
            ++idx;
        }
        jniInitializeHaplotypes(numHaplotypes, mHaplotypeDataArray);
    }
    /**
     * Tell JNI to release arrays - really important if native code is directly accessing Java memory, if not
     * accessing Java memory directly, still important to release memory from C++
     */
    private native void jniFinalizeRegion();

    
    public void finalizeRegion()
    {
        jniFinalizeRegion();
    }

    /**
     * Real compute kernel
     */
    private native void jniComputeLikelihoods(int numReads, int numHaplotypes, JNIReadDataHolderClass[] readDataArray,
            JNIHaplotypeDataHolderClass[] haplotypeDataArray, double[] likelihoodArray, int maxNumThreadsToUse);
    
    public void computeLikelihoods(final List<JNIReadDataHolderClass> reads, final List<JNIHaplotypeDataHolderClass> haplotypes) {
      //System.out.println("Region : "+reads.size()+" x "+haplotypes.size());
      long startTime = System.nanoTime();
      int readListSize = reads.size();
      int numHaplotypes = haplotypes.size();
      int numTestcases = readListSize*numHaplotypes;
      JNIReadDataHolderClass[] readDataArray = new JNIReadDataHolderClass[readListSize];
      int idx = 0;
      for(JNIReadDataHolderClass read : reads)
      {
        readDataArray[idx] = new JNIReadDataHolderClass();
        readDataArray[idx].readBases    = read.readBases; 
        readDataArray[idx].readQuals    = read.readQuals;
        readDataArray[idx].insertionGOP = read.insertionGOP;
        readDataArray[idx].deletionGOP  = read.deletionGOP;
        readDataArray[idx].overallGCP   = read.overallGCP;
        ++idx;
      }

      double[] mLikelihoodArray = new double[readListSize*numHaplotypes];      //to store results
      setupTime += (System.nanoTime() - startTime);
      //for(reads)
      //   for(haplotypes)
      //       compute_full_prob()
      jniComputeLikelihoods(readListSize, numHaplotypes, readDataArray, mHaplotypeDataArray, mLikelihoodArray, 12);

      computeTime += (System.nanoTime() - startTime);
    }
    
    /**
     * Print final profiling information from native code
     */
    public native void jniClose();
    public void close()
    {
        System.err.println("Time spent in setup for JNI call : " + (setupTime * 1e-9) + " compute time : " + (computeTime * 1e-9));
        jniClose();
    }

    public void parseSandboxFile(String filename)
    {
      File file = new File(filename);
      Scanner input = null;
      try
      {
        input = new Scanner(file); 
      }
      catch(FileNotFoundException e)
      {
          System.err.println("File "+filename + " cannot be found/read");
          return;
      }
      int idx = 0;
      int numReads = 0;
      int numHaplotypes = 0;
      int readIdx = 0, testCaseIdx = 0, haplotypeIdx = 0;
      LinkedList<JNIHaplotypeDataHolderClass> haplotypeList = new LinkedList<JNIHaplotypeDataHolderClass>();
      LinkedList<JNIReadDataHolderClass> readList = new LinkedList<JNIReadDataHolderClass>();
      
      byte[][] byteArray = new byte[6][];
      boolean firstLine = true;
      String[] currTokens = new String[8];
      while(input.hasNextLine())
      {
        String line = input.nextLine();
        Scanner lineScanner = new Scanner(line);
        idx = 0;
        while(lineScanner.hasNext())
          currTokens[idx++] = lineScanner.next();
        if(idx == 0)
          break;
        assert(idx >= 6);
        //start of new region
        if(idx == 8)
        {
          if(!firstLine)
          {
            initialize(haplotypeList);
            computeLikelihoods(readList, haplotypeList);
            finalizeRegion();
          }
          try
          {
            numReads = Integer.parseInt(currTokens[6]);
          }
          catch(NumberFormatException e)
          {
            numReads = 1;
          }
          try
          {
            numHaplotypes = Integer.parseInt(currTokens[7]);
          }
          catch(NumberFormatException e)
          {
            numHaplotypes = 1;
          }
          haplotypeIdx = readIdx = testCaseIdx = 0;
          readList.clear();
          haplotypeList.clear();
        }
        if(haplotypeIdx < numHaplotypes)
        {
          JNIHaplotypeDataHolderClass X = new JNIHaplotypeDataHolderClass();
          X.haplotypeBases = currTokens[0].getBytes();
          haplotypeList.add(X);
        }
        if(testCaseIdx%numHaplotypes == 0)
        {
          JNIReadDataHolderClass X = new JNIReadDataHolderClass();
          X.readBases = currTokens[1].getBytes();
          for(int i=2;i<6;++i)
          {
            byteArray[i] = currTokens[i].getBytes();
            for(int j=0;j<byteArray[i].length;++j)
              byteArray[i][j] -= 33;       //normalize
          }
          X.readQuals = byteArray[2];
          X.insertionGOP = byteArray[3];
          X.deletionGOP = byteArray[4];
          X.overallGCP = byteArray[5];
          readList.add(X);
        }
        ++testCaseIdx;
        ++haplotypeIdx;

        lineScanner.close();
        firstLine = false;
      }
      if(haplotypeList.size() > 0 && readList.size() > 0)
      {
        initialize(haplotypeList);
        computeLikelihoods(readList, haplotypeList);
        finalizeRegion();
      }

      close();
      input.close();
    }

    private native void doEverythingNative(String filename);

    public static void main(String[] args)
    {
      if(args.length <= 0)
      {
        System.err.println("Needs 1 argument - <filename>");
        System.exit(-1);
      }
      //// Get runtime
      //java.lang.Runtime rt = java.lang.Runtime.getRuntime();
      //// Start a new process: UNIX command ls
      //String cmd = "/home/karthikg/broad/gsa-unstable/public/c++/VectorPairHMM/checker "+args[0];
      //try
      //{
        //System.out.println(cmd);
        //java.lang.Process p = rt.exec(cmd);
        //try
        //{
          //p.waitFor();
          //java.io.InputStream is = p.getInputStream();
          //java.io.BufferedReader reader = new java.io.BufferedReader(new InputStreamReader(is));
          //// And print each line
          //String s = null;
          //while ((s = reader.readLine()) != null) {
            //System.out.println(s);
          //}
          //is.close();
        //}
        //catch(InterruptedException e)
        //{
          //System.err.println(e);
        //}
      //}
      //catch(IOException e)
      //{
        //System.err.println(e);
      //}
      Sandbox t = new Sandbox();
      //t.doEverythingNative(args[0]);
      t.parseSandboxFile(args[0]);
    }
}
