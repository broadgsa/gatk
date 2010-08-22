/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.Vector;
import java.util.Collections;
import java.io.PrintStream;

/**
 * Calculates likelihood of observing the data given pairs of HLA alleles. NOTE: run CalculateBaseLikelihoods first! Usage: java -jar GenomeAnalysisTK.jar -T CalculateAlleleLikelihoods -I /humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.imputed.4digit.bam -R /broad/1KG/reference/human_b36_both.fasta -L /humgen/gsa-scr1/GSA/sjia/454_HLA/HAPMAP270/HLA_exons.interval -bl INPUT.baselikelihoods -eth\
nicity Caucasian | grep -v "INFO"  | grep -v "DONE!" > OUTPUT.allelelikelihoods
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CalculateAlleleLikelihoodsWalker extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(fullName = "baseLikelihoods", shortName = "bl", doc = "Base likelihoods file", required = true)
    public String baseLikelihoodsFile = "";

    @Argument(fullName = "debugHLA", shortName = "debugHLA", doc = "Print debug", required = false)
    public boolean DEBUG = false;

    @Argument(fullName = "debugAlleles", shortName = "debugAlleles", doc = "Print likelihood scores for these alleles", required = false)
    public String debugAlleles = "";

    @Argument(fullName = "onlyfrequent", shortName = "onlyfrequent", doc = "Only consider alleles with frequency > 0.0001", required = false)
    public boolean FREQUENT = false;

    @Argument(fullName = "ethnicity", shortName = "ethnicity", doc = "Use allele frequencies for this ethnic group", required = false)
    public String ethnicity = "CaucasianUSA";

    String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_Caucasians.freq";
    String BlackAlleleFrequencyFile     = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_BlackUSA.freq";
    String AlleleFrequencyFile          = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";
    String UniqueAllelesFile            = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/UniqueAlleles4Digit";
    String HLAdatabaseFile              = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_DICTIONARY.txt";
    String HLA2DigitFile              = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_DICTIONARY_2DIGIT.txt";

    Hashtable AlleleFrequencies,UniqueAlleles,Alleles2Digit;
    
    CigarParser formatter = new CigarParser();
    double[][] baseLikelihoods;
    int[] positions;
    boolean loaded = false;

    String[] HLAnames, HLAreads, HLAnames2, HLAreads2;
    Integer[] HLAstartpos, HLAstoppos, HLAstartpos2, HLAstoppos2;
    ArrayList<String> HLAnamesAL, HLAreadsAL, Loci, AllelesToSearch;
    ArrayList<Integer> HLAstartposAL, HLAstopposAL;

    public Integer reduceInit() {
        if (!loaded){
            loaded = true;
            BaseLikelihoodsFileReader baseLikelihoodsReader = new BaseLikelihoodsFileReader();
            baseLikelihoodsReader.ReadFile(baseLikelihoodsFile, false);
            baseLikelihoods = baseLikelihoodsReader.GetBaseLikelihoods();
            positions = baseLikelihoodsReader.GetPositions();

            HLAnamesAL = new ArrayList<String>();
            HLAreadsAL = new ArrayList<String>();
            HLAstartposAL = new ArrayList<Integer>();
            HLAstopposAL = new ArrayList<Integer>();

            out.printf("INFO  Reading HLA alleles ... ");
            HLAFileReader HLADictionaryReader = new HLAFileReader();
            HLADictionaryReader.ReadFile(HLAdatabaseFile);
            HLAreads = HLADictionaryReader.GetSequences();
            HLAnames = HLADictionaryReader.GetNames();
            HLAstartpos = HLADictionaryReader.GetStartPositions();
            HLAstoppos = HLADictionaryReader.GetStopPositions();

            HLADictionaryReader = new HLAFileReader();
            HLADictionaryReader.ReadFile(HLA2DigitFile);
            HLAreads2 = HLADictionaryReader.GetSequences();
            HLAnames2 = HLADictionaryReader.GetNames();
            HLAstartpos2 = HLADictionaryReader.GetStartPositions();
            HLAstoppos2 = HLADictionaryReader.GetStopPositions();
            out.printf("Done! %s HLA alleles loaded.\n",HLAreads.length);

            //out.printf("INFO Common alleles:\n");
            for (int i = 1; i < UniqueAlleles.size(); i++){
                //out.printf("INFO %s\n",UniqueAlleles.values().toArray()[i]);
            }
            //out.printf("INFO  Reading HLA dictionary ...");


        }
        return 0;
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        //HLAnamesAL.add(read.getReadName());
        //HLAreadsAL.add(formatter.FormatRead(read.getCigarString(), read.getReadString()));
        //HLAstartposAL.add(read.getAlignmentStart());
        //HLAstopposAL.add(read.getAlignmentEnd());
        //out.printf("INFO\t%s\t%s\t%s\t%s\n",read.getReadName(),read.getAlignmentStart(),read.getAlignmentEnd(),formatter.FormatRead(read.getCigarString(), read.getReadString()));
        return 1;
    }

    private int GenotypeIndex(char a, char b){
        switch(a){
            case 'A':
                switch(b){
                    case 'A': return 0;
                    case 'C': return 1;
                    case 'G': return 2;
                    case 'T': return 3;
                };
            case 'C':
                switch(b){
                    case 'A': return 1;
                    case 'C': return 4;
                    case 'G': return 5;
                    case 'T': return 6;
                };
            case 'G':
                switch(b){
                    case 'A': return 2;
                    case 'C': return 5;
                    case 'G': return 7;
                    case 'T': return 8;
                };
            case 'T':
                switch(b){
                    case 'A': return 3;
                    case 'C': return 6;
                    case 'G': return 8;
                    case 'T': return 9;
                };
            default: return -1;
        }
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer numreads) {
        //out.printf("Done! %s alleles found\n", numreads);
        //HLAnames = HLAnamesAL.toArray(new String[numreads]);
        //HLAreads = HLAreadsAL.toArray(new String[numreads]);
        //HLAstartpos = HLAstartposAL.toArray(new Integer[numreads]);
        //HLAstoppos = HLAstopposAL.toArray(new Integer[numreads]);

        double[][] AlleleLikelihoods = new double[numreads][numreads];

        String name1, name2;
        double frq1, frq2;


        double minfrq = 0;
        if (FREQUENT){
            minfrq = 0.0001;
        }
        int numcombinations = 0;
        out.printf("NUM\tAllele1\tAllele2\tSSG\n");

        //debugging specific alleles
        int index1 = -1, index2 = -1;
        if (!debugAlleles.equals("")){
            String s[] = debugAlleles.split(",");
            for (int i = 0; i < numreads; i++){
                if (HLAnames[i].equals(s[0])){
                    index1 = i;
                }
                if (HLAnames[i].equals(s[1])){
                    index2 = i;
                }
                if (index1 > -1 && index2 > -1){
                    out.printf("INFO: debugging %s\t%s\t%s\t%s\n",s[0],s[1],index1,index2);
                    double dl = CalculateLikelihood(index1,index2,HLAreads2,true);
                    break;
                }
            }
        }

        //Pre-process homozygous combinations to determine top possible alleles (for efficiency)
        int numreads2 = HLAnames2.length;
        Alleles2Digit = new Hashtable();
        Loci = new ArrayList<String>();
        double[] AlleleLikelihoods2 = new double[numreads];
        for (int i = 0; i < numreads; i++){
            name1 = HLAnames[i].substring(4);
            String [] n1 = name1.split("\\*");
            numcombinations++;
            AlleleLikelihoods2[i] = CalculateLikelihood(i,i,HLAreads,false);
            if (AlleleLikelihoods2[i] < 0){
                name2 = n1[0] + "*" + n1[1].substring(0, 2);
                if (!Loci.contains(n1[0])){Loci.add(n1[0]);}
                if (!Alleles2Digit.containsKey(name2)){
                    Alleles2Digit.put(name2, AlleleLikelihoods2[i]);
                }else if ((Double) Alleles2Digit.get(name2) < AlleleLikelihoods2[i]){
                    Alleles2Digit.put(name2, AlleleLikelihoods2[i]);
                }
            }
        }

        //Sort alleles at 2 digit resolution for each locus
        AllelesToSearch = new ArrayList<String>();
        for (int i = 0; i < Loci.size(); i++){
            Enumeration k = Alleles2Digit.keys();
            Hashtable AllelesAtLoci = new Hashtable();

            //find alleles at the locus
            while( k.hasMoreElements() ){
                name1 = k.nextElement().toString();
                String [] n1 = name1.split("\\*");
                if (Loci.get(i).equals(n1[0])){
                    AllelesAtLoci.put(-1 * (Double) Alleles2Digit.get(name1), name1);
                }
            }

            //Sort alleles at locus, mark top six 2-digit classes for deep search
            int num = 1;
            Vector v = new Vector(AllelesAtLoci.keySet());
            Collections.sort(v);
            for (Enumeration e = v.elements(); e.hasMoreElements();) {
                Double key = Double.valueOf(e.nextElement().toString());
                String allele = AllelesAtLoci.get(key).toString();
                if (num <= 10){
                    AllelesToSearch.add(allele);
                    
                    num++;
                }
                //out.printf("%s\t%s\n",allele,key);
            }
        }

        //Iterate through allele pairs to calculate likelihoods
        if (true){
            numcombinations = 0;
            for (int i = 0; i < numreads; i++){
                name1 = HLAnames[i].substring(4);
                String [] n1 = name1.split("\\*");
                if (AllelesToSearch.contains(n1[0] + "*" + n1[1].substring(0, 2))){
                    //out.printf("1: %s\n",name1);
                    //frq1 = Double.parseDouble((String) AlleleFrequencies.get(name1).toString());
                    //if (frq1 > minfrq){
                    for (int j = i; j < numreads; j++){
                        name2 = HLAnames[j].substring(4);
                        String [] n2 = name2.split("\\*");
                        if (AllelesToSearch.contains(n2[0] + "*" + n2[1].substring(0, 2))){
                            if ((HLAstartpos[i] < HLAstoppos[j]) && (HLAstartpos[j] < HLAstoppos[i])){
                                numcombinations++;
                                AlleleLikelihoods[i][j] = CalculateLikelihood(i,j,HLAreads,false);
                                if (AlleleLikelihoods[i][j] < 0){
                                    out.printf("%s\t%s\t%s\t%.2f\n",numcombinations,name1,name2,AlleleLikelihoods[i][j]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private double CalculateLikelihood(int a1, int a2, String[] HLAalleles, boolean debug){
        //Calculates likelihood for specific allele pair
        String read1 = HLAalleles[a1];
        String read2 = HLAalleles[a2];
        int start1 = HLAstartpos[a1];
        int start2 = HLAstartpos[a2];
        int stop1 = HLAstoppos[a1];
        int stop2 = HLAstoppos[a2];
        double likelihood = 0;
        int pos, index;
        char c1, c2;
        

        for (int i = 0; i < positions.length; i++){
            pos = positions[i];
            if (pos < stop1 && pos > start1 && pos < stop2 && pos > start2){
                index = GenotypeIndex(read1.charAt(pos-start1),read2.charAt(pos-start2));
                if (index > -1){
                    likelihood = likelihood + baseLikelihoods[i][index];
                    if (debug){
                        c1 = read1.charAt(pos-start1);
                        c2 = read2.charAt(pos-start2);
                        out.printf("INFO: DEBUG %s\t%s\t%s\t%s\t%s\t%s\t%.2f\n",HLAnames[a1],HLAnames[a2],pos,c1,c2,index,likelihood);
                    }
                }
            }
        }
        return likelihood;
    }
}
