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
import org.broadinstitute.sting.commandline.Argument;

import java.util.ArrayList;
import java.util.Hashtable;


/**
 * Calculates likelihood of observing the data given pairs of HLA alleles. NOTE: run CalculateBaseLikelihoods first! Usage: java -jar GenomeAnalysisTK.jar -T CalculateAlleleLikelihoods -I /humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.imputed.4digit.bam -R /broad/1KG/reference/human_b36_both.fasta -L /humgen/gsa-scr1/GSA/sjia/454_HLA/HAPMAP270/HLA_exons.interval -bl INPUT.baselikelihoods -eth\
nicity Caucasian | grep -v "INFO"  | grep -v "DONE!" > OUTPUT.allelelikelihoods
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CalculateAlleleLikelihoodsWalker extends ReadWalker<Integer, Integer> {
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
    String UniqueAllelesFile            = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/UniqueAllelesCommon";
    Hashtable AlleleFrequencies,UniqueAlleles;
    
    CigarParser formatter = new CigarParser();
    double[][] baseLikelihoods;
    int[] positions;
    boolean loaded = false;

    String[] HLAnames, HLAreads;
    Integer[] HLAstartpos, HLAstoppos;
    ArrayList<String> HLAnamesAL, HLAreadsAL;
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

            if (!ethnicity.equals("CaucasianUSA")){
                AlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_" + ethnicity + ".freq";
            }
            out.printf("INFO  Reading HLA allele frequencies ... ");
            FrequencyFileReader HLAfreqReader = new FrequencyFileReader();
            HLAfreqReader.ReadFile(AlleleFrequencyFile,UniqueAllelesFile);
            AlleleFrequencies = HLAfreqReader.GetAlleleFrequencies();
            UniqueAlleles = HLAfreqReader.GetUniqueAlleles();
            out.printf("Done! Frequencies for %s HLA alleles loaded.\n",AlleleFrequencies.size());

            //out.printf("INFO Common alleles:\n");
            for (int i = 1; i < UniqueAlleles.size(); i++){
                //out.printf("INFO %s\n",UniqueAlleles.values().toArray()[i]);
            }
            out.printf("INFO  Reading HLA dictionary ...");


        }
        return 0;
    }

    public Integer map(char[] ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        HLAnamesAL.add(read.getReadName());
        HLAreadsAL.add(formatter.FormatRead(read.getCigarString(), read.getReadString()));
        HLAstartposAL.add(read.getAlignmentStart());
        HLAstopposAL.add(read.getAlignmentEnd());
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
        out.printf("Done! %s alleles found\n", numreads);
        HLAnames = HLAnamesAL.toArray(new String[numreads]);
        HLAreads = HLAreadsAL.toArray(new String[numreads]);
        HLAstartpos = HLAstartposAL.toArray(new Integer[numreads]);
        HLAstoppos = HLAstopposAL.toArray(new Integer[numreads]);

        double[][] AlleleLikelihoods = new double[numreads][numreads];

        String name1, name2;
        double frq1, frq2;


        double minfrq = 0;
        if (FREQUENT){
            minfrq = 0.0001;
        }
        int numcombinations = 0;
        out.printf("NUM\tAllele1\tAllele2\tSSG\n");

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
                    double dl = CalculateLikelihood(index1,index2,true);
                    break;
                }
            }
        }
        
        for (int i = 0; i < numreads; i++){
            name1 = HLAnames[i].substring(4);
            String [] n1 = name1.split("\\*");
//            out.printf("1: %s\n",n1[0] + "*" + n1[1].substring(0, 3));
            if (UniqueAlleles.containsKey(n1[0] + "*" + n1[1].substring(0, 4))){
                //out.printf("1: %s\n",name1);
                //frq1 = Double.parseDouble((String) AlleleFrequencies.get(name1).toString());
                //if (frq1 > minfrq){
                    for (int j = i; j < numreads; j++){
                        name2 = HLAnames[j].substring(4);
                        String [] n2 = name2.split("\\*");
//                        out.printf("2: %s\n",n2[0] + "*" + n2[1].substring(0, 3));
                        if (UniqueAlleles.containsKey(n2[0] + "*" + n2[1].substring(0, 4))){
                        
                //            frq2 = Double.parseDouble((String) AlleleFrequencies.get(name2).toString());
                //            if (frq2 > minfrq){
                                if ((HLAstartpos[i] < HLAstoppos[j]) && (HLAstartpos[j] < HLAstoppos[i])){
                                    numcombinations++;
                                    AlleleLikelihoods[i][j] = CalculateLikelihood(i,j,false);
                                    out.printf("%s\t%s\t%s\t%.2f\n",numcombinations,name1,name2,AlleleLikelihoods[i][j]);
                                }
                //            }else{
                //                if (DEBUG){out.printf("%s has allele frequency%.5f\n",name2,frq2);}
                //            }
                //        }else{
                //            if (DEBUG){out.printf("%s not found in allele frequency file\n",name2);}
                        }
                    }
                //}else{
                //    if (DEBUG){out.printf("%s has allele frequency%.5f\n",name1,frq1);}
                //}
            //}else{
            //    if (DEBUG){out.printf("%s not found in allele frequency file\n",name1);}
            }
        }
    }

    private double CalculateLikelihood(int a1, int a2, boolean debug){
        String read1 = HLAreads[a1];
        String read2 = HLAreads[a2];
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
