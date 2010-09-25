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
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.*;
import java.util.Map.Entry;
import java.io.PrintStream;

/**
 * Calculates the likelihood of observing data given phase info from pairs of HLA alleles. Note: Run FindClosestAlleleWalker first! Usage: java -jar $GATK -T HLACaller -I INPUT.bam -R /broad/1KG/reference/human_b36_both.fasta -L /humgen/gsa-scr1/GSA/sjia/454_HLA/HAPMAP270/HLA_exons.interval -phaseInterval /humgen/gsa-scr1/GSA/sjia/454_HLA/HAPMAP270/HLA_exons.interval -bl IMPUT.baselikelihoods [-filter $ID.filter -minAllowe\
dMismatches 7] -ethnicity Caucasian | grep -v "INFO"  | grep -v "DEBUG" | grep -v "DONE!" > OUTPUT.phaselikelihoods
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class HLACallerWalker extends ReadWalker<Integer, Integer> {
    @Output
    private PrintStream out;

    @Argument(fullName = "baseLikelihoods", shortName = "bl", doc = "Base likelihoods file", required = true)
    public String baseLikelihoodsFile = "";

    @Argument(fullName = "debugHLA", shortName = "debugHLA", doc = "Print debug", required = false)
    public boolean DEBUG = false;

    @Argument(fullName = "filter", shortName = "filter", doc = "file containing reads to exclude", required = false)
    public String filterFile = "";

    @Argument(fullName = "ethnicity", shortName = "ethnicity", doc = "Use allele frequencies for this ethnic group", required = false)
    public String ethnicity = "CaucasiansUSA";

    @Argument(fullName = "debugAlleles", shortName = "debugAlleles", doc = "Print likelihood scores for these alleles", required = false)
    public String debugAlleles = "";

    @Argument(fullName = "useInterval", shortName = "useInterval", doc = "Use only these intervals in phase calculation", required = false)
    public String IntervalsFile = "";

    @Argument(fullName = "minFreq", shortName = "minFreq", doc = "only consider alleles greater than this frequency", required = false)
    public double minFrequency = 0.0;

    @Argument(fullName = "maxAllowedMismatches", shortName = "maxAllowedMismatches", doc = "Max number of mismatches tolerated per read (default 7)", required = false)
    public int MAXALLOWEDMISMATCHES = 6;

    @Argument(fullName = "minRequiredMatches", shortName = "minRequiredMatches", doc = "Min number of matches required per read (default 7)", required = false)
    public int MINREQUIREDMATCHES = 5;

    @Argument(fullName = "HLAfrequencies", shortName = "HLAfrequencies", doc = "HLA allele frequencies file", required = true)
    public String AlleleFrequencyFile = "HLA_FREQUENCIES.txt";

    @Argument(fullName = "HLAdictionary", shortName = "HLAdictionary", doc = "HLA dictionary file", required = true)
    public String HLAdatabaseFile = "HLA_DICTIONARY.txt";

    @Argument(fullName = "turnOffVerboseOutput", shortName = "noVerbose", doc = "Do not output verbose probability descriptions (INFO lines) ", required = false)
    protected boolean NO_VERBOSE = false;

    // Initializing variables
    
    HLAFileReader HLADictionaryReader = new HLAFileReader();
    boolean HLAdataLoaded = false;
    String[] HLAnames, HLAreads, Populations;
    ArrayList<String> ReadsToDiscard;
    Integer[] HLAstartpos, HLAstoppos, PolymorphicSites;


    int[][] numObservations, totalObservations, intervals;
    int[] SNPnumInRead, SNPposInRead, positions;
    CigarParser cigarparser = new CigarParser();
    Hashtable MaxLikelihoods = new Hashtable();
    Hashtable MaxFrequencies, CommonAlleles, AlleleCount, LocusCount;
    Hashtable[] AlleleFrequencies;
    int numIntervals;
    double[][] baseLikelihoods;

    ArrayList AllelesToSearch = new ArrayList<String>();

    // setting error rates for phasing algorithm (1% expected error rate for any genotype)
    
    double P_err = 0.01;
    double P_correct = 1 - P_err;
    double L_err = Math.log10(P_err);
    double L_correct = Math.log10(P_correct);

    public Integer reduceInit() {

        if (!HLAdataLoaded){
            HLAdataLoaded = true;

            //Load HLA dictionary

            HLADictionaryReader.ReadFile(HLAdatabaseFile);
            HLAreads = HLADictionaryReader.GetSequences();
            HLAnames = HLADictionaryReader.GetNames();
            HLAstartpos = HLADictionaryReader.GetStartPositions();
            HLAstoppos = HLADictionaryReader.GetStopPositions();

            //Load pre-processing file for misaligned reads and list of alleles to search

            if (!filterFile.equals("")){
                //If pre-processing file exists, load contents
                SimilarityFileReader similarityReader = new SimilarityFileReader();
                similarityReader.ReadFile(filterFile,MAXALLOWEDMISMATCHES,MINREQUIREDMATCHES);
                ReadsToDiscard = similarityReader.GetReadsToDiscard();
                AllelesToSearch = similarityReader.GetAllelesToSearch();
                AlleleCount = similarityReader.GetAlleleCount();
                LocusCount = similarityReader.GetLocusCount();
                if (!NO_VERBOSE) {
                    for (int i = 0; i < AllelesToSearch.size(); i++){
                        out.printf("INFO\tAllelesToSearch\t%s\t%s\n",AllelesToSearch.get(i),AlleleCount.get(AllelesToSearch.get(i)));
                    }
                }
            }else{
                ReadsToDiscard = new ArrayList<String>();
                AlleleCount = new Hashtable();
                String name, d4_name; String [] n;
                for (int i = 0; i < HLAnames.length; i++){
                    name = HLAnames[i].substring(4);
                    n = name.split("\\*");
                    d4_name = n[0] + "*" + n[1].substring(0, 4);
                    if (!AllelesToSearch.contains(d4_name)){
                        AllelesToSearch.add(d4_name);
                        AlleleCount.put(d4_name, 0);
                    }
                    if (!LocusCount.containsKey(n[0])){
                        LocusCount.put(n[0], 0);
                    }
                }
            }

            //Load genotypes and find polymorphic sites (sites that differ from reference)
            
            BaseLikelihoodsFileReader baseLikelihoodsReader = new BaseLikelihoodsFileReader();
            baseLikelihoodsReader.ReadFile(baseLikelihoodsFile, true);
            baseLikelihoods = baseLikelihoodsReader.GetBaseLikelihoods();
            positions = baseLikelihoodsReader.GetPositions();
            PolymorphicSites = baseLikelihoodsReader.GetPolymorphicSites();
            if (!NO_VERBOSE) {
                out.printf("INFO\t%s polymorphic sites found\n",PolymorphicSites.length);
            }

            int l = PolymorphicSites.length;
            SNPnumInRead = new int[l];
            SNPposInRead = new int[l];
            numObservations = new int[l*5][l*5];
            totalObservations = new int[l][l];

            //Load allele frequencies for different populations

            FrequencyFileReader HLAfreqReader = new FrequencyFileReader();
            HLAfreqReader.ReadFile(AlleleFrequencyFile,ethnicity);
            AlleleFrequencies = HLAfreqReader.GetAlleleFrequencies();
            MaxFrequencies = HLAfreqReader.GetMaxFrequencies();
            CommonAlleles = HLAfreqReader.GetCommonAlleles();
            Populations = HLAfreqReader.GetPopulations();

            //Load genomic intervals for bam file
            
            if (!IntervalsFile.equals("")){
                TextFileReader fileReader = new TextFileReader();
                fileReader.ReadFile(IntervalsFile);
                String[] lines = fileReader.GetLines();
                intervals = new int[lines.length][2];
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split(":");
                    String[] intervalPieces = s[1].split("-");
                    intervals[i][0] = Integer.valueOf(intervalPieces[0]);
                    intervals[i][1] = Integer.valueOf(intervalPieces[1]);
                }
                numIntervals = intervals.length;
            }

            
        }
        return 0;
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        if (!ReadsToDiscard.contains(read.getReadName())){
            UpdateCorrelation(read);
        }else{
            
        }
        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer numreads) {
        String name1, name2, d4_name1, d4_name2, d2_name1, d2_name2;
        Double frq1 = 0.0, frq2 = 0.0, log1 = 0.0, log2 = 0.0,alleleLikelihood= 0.0, phaseLikelihood=0.0, likelihood = 0.0;
        int numCombinations = 0;

        //For debugging specific alleles
        if (!debugAlleles.equals("")){
            String s[] = debugAlleles.split(",");
            int index1 = HLADictionaryReader.GetIndex(s[0]);
            int index2 = HLADictionaryReader.GetIndex(s[1]);
            out.printf("INFO: debugging %s\t%s\t%s\t%s\n",s[0],s[1],index1,index2);
            if (index1 > -1 && index2 > -1){
                alleleLikelihood = CalculateAlleleLikelihood(index1,index2,HLAreads,true);
                phaseLikelihood = CalculatePhaseLikelihood(index1,index2,true,false);
            }
        }

        double max;
        ArrayList Output = new ArrayList<String>();
        ArrayList Likelihoods = new ArrayList<Double>();
        Hashtable TotalProb = new Hashtable();
        //Search pairs of alleles that satisfy initial search criteria

        // Allele 1

        for (int i = 0; i < HLAnames.length; i++){
            name1 = HLAnames[i].substring(4);
            String [] n1 = name1.split("\\*");
            d4_name1 = n1[0] + "*" + n1[1].substring(0, 4);
            d2_name1 = n1[0] + "*" + n1[1].substring(0, 2);
            if (AllelesToSearch.contains(d4_name1)){
                
                if (MaxFrequencies.containsKey(d4_name1)){
                    frq1 = Double.parseDouble(MaxFrequencies.get(d4_name1).toString());
                }else{
                    if (n1[1].length() > 4){if (n1[1].substring(4, 5).equals("N")){frq1 = .00000005;}else{frq1 = .000001;}}else{frq1 = .000001;}
                }

                if (frq1 > minFrequency){
                    
                    // Allele 2

                    for (int j = i; j < HLAnames.length; j++){
                        name2 = HLAnames[j].substring(4);
                        String [] n2 = name2.split("\\*");
                        d4_name2 = n2[0] + "*" + n2[1].substring(0, 4);
                        d2_name2 = n2[0] + "*" + n2[1].substring(0, 2);
                        if (n1[0].equals(n2[0]) && (AllelesToSearch.contains(d4_name2))){
                            if (MaxFrequencies.containsKey(d4_name2)){
                                frq2 = Double.parseDouble(MaxFrequencies.get(d4_name2).toString());
                            }else{
                                if (n2[1].length() > 4){if (n2[1].substring(4, 5).equals("N")){frq2 = .00000005;}else{frq2 = .000001;}}else{frq2 = .000001;}
                            }

                            if (frq2 > minFrequency){
                                
                            //Calculate allele and phase likelihoods for each allele pair
                                alleleLikelihood = CalculateAlleleLikelihood(i,j,HLAreads,false);
                                numCombinations++;

                                //If there is data at the allele pair, continue with other calculations

                                if (alleleLikelihood < 0){
                                    phaseLikelihood = CalculatePhaseLikelihood(i,j,false,false);
                                    log1=Math.log10(frq1);
                                    log2=Math.log10(frq2);

                                    //sum likelihoods

                                    likelihood = alleleLikelihood+phaseLikelihood+log1+log2;
                                    if (!MaxLikelihoods.containsKey(n1[0])){MaxLikelihoods.put(n1[0], likelihood);}

                                    if (likelihood > (Double) MaxLikelihoods.get(n1[0])) {
                                        MaxLikelihoods.put(n1[0], likelihood);
                                    }
                                    Likelihoods.add(likelihood);
                                    String data = String.format("%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",n1[0],name1,name2,alleleLikelihood,phaseLikelihood,log1,log2,likelihood);
                                    Output.add(data);
                                    if (!NO_VERBOSE) {
                                        out.printf("INFO\t%s\n",data);
                                    }
                                    if (DEBUG){

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //Print output
        out.printf("Locus\tA1\tA2\tGeno\tPhase\tFrq1\tFrq2\tL\tProb\tReads1\tReads2\tLocus\tEXP");
        for (int i = 0; i < Populations.length; i++){
            out.printf("\t%s",Populations[i]);
        }
        out.printf("\n");

        //Calculate probabilities for each locus
        Double probSum = 0.0, prob = 0.0, f1 = 0.0, f2 = 0.0, aLikelihood4 = 0.0, pLikelihood4 = 0.0;
        Integer count = 0;
        Hashtable HLA4DigitProbs = new Hashtable();
        Hashtable HLA4DigitLs = new Hashtable();
        Hashtable HLA4DigitCount = new Hashtable();
        Hashtable HLA4DigitF1 = new Hashtable();
        Hashtable HLA4DigitF2 = new Hashtable();
        Hashtable HLA4DigitA = new Hashtable();
        Hashtable HLA4DigitP = new Hashtable();

        String key;
        Enumeration keys = LocusCount.keys();
        while (keys.hasMoreElements()){
            String locus = keys.nextElement().toString();
            
            probSum = 0.0;
            ArrayList localOutput = new ArrayList<String>();
            ArrayList localLikelihoods = new ArrayList<Double>();

            //Sum probabilities for each locus

            for (int j = 0; j < Output.size(); j++){
                String data = Output.get(j).toString();
                String [] d = data.split("\\t");
                if (d[0].equals(locus)){
                    localOutput.add(data);
                    likelihood = (Double)Likelihoods.get(j)-(Double)MaxLikelihoods.get(locus);
                    localLikelihoods.add(likelihood);
                    probSum = probSum + Math.pow(10, likelihood);
                    //out.printf("INFO\t%s\t%s\t%.2f\t%.2f\t%.2f\t%s\n",locus,data,likelihood,(Double)MaxLikelihoods.get(locus),(Double)Likelihoods.get(j),probSum);
                }
            }
            
            //aggregate statistics for 4-digit types
            
            String A1 = "", A2 = "", a1 = "", a2 = "";
            String [] s, s1, s2;
            Double prob4digit = 0.0;
            int n = 0;

            for (int j = 0; j < localOutput.size(); j++){
                String data = localOutput.get(j).toString();
                prob = Math.pow(10, (Double)localLikelihoods.get(j))/probSum;
                
                if (prob > 0.005){
                    s = data.split("\\t");
                    s1 = s[1].split("\\*");
                    s2 = s[2].split("\\*");
                    a1 = s1[0] + "*" + s1[1].substring(0,4);
                    a2 = s2[0] + "*" + s2[1].substring(0,4);
                    key = a1 + "," + a2;
                    aLikelihood4 = Double.valueOf(s[3]);
                    pLikelihood4 = Double.valueOf(s[4]);
                    likelihood = aLikelihood4 + pLikelihood4 + f1 + f2;
                    f1 = Double.valueOf(s[5]);
                    f2 = Double.valueOf(s[6]);
                    if (!HLA4DigitProbs.containsKey(key)){
                        HLA4DigitProbs.put(key, prob);
                        HLA4DigitLs.put(key, likelihood);
                        HLA4DigitCount.put(key, 1);
                        HLA4DigitF1.put(key,f1);
                        HLA4DigitF2.put(key,f2);
                        HLA4DigitA.put(key,aLikelihood4);
                        HLA4DigitP.put(key,pLikelihood4);
                    }else{
                        prob = prob + Double.valueOf(HLA4DigitProbs.get(key).toString());
                        HLA4DigitProbs.put(key, prob);
                        likelihood = likelihood + Double.valueOf(HLA4DigitLs.get(key).toString());
                        HLA4DigitLs.put(key, likelihood);
                        n = Integer.valueOf(HLA4DigitCount.get(key).toString()) + 1;
                        HLA4DigitCount.put(key, n);
                        aLikelihood4 = aLikelihood4 + Double.valueOf(HLA4DigitA.get(key).toString());
                        HLA4DigitA.put(key, aLikelihood4);
                        pLikelihood4 = pLikelihood4 + Double.valueOf(HLA4DigitP.get(key).toString());
                        HLA4DigitP.put(key, pLikelihood4);
                    }
                    
                }
            }
        }

        //Print results
        Enumeration P = HLA4DigitProbs.keys();
        String K = ""; String [] s, s1, s2;
        double count1, count2, locusCount, accountedFor;

        // Sort hashtable.
        Vector v = new Vector(HLA4DigitProbs.keySet());
        Collections.sort(v);

        // Display (sorted) hashtable.
        for (Enumeration e = v.elements(); e.hasMoreElements();) {
            K = (String)e.nextElement();
            prob = (Double) HLA4DigitProbs.get(K);

            likelihood = (Double) HLA4DigitLs.get(K);
            count = (Integer) HLA4DigitCount.get(K);
            s = K.split("\\,");
            s1 = s[0].split("\\*"); name1 = s1[1];
            s2 = s[1].split("\\*"); name2 = s2[1];
            aLikelihood4 = (Double) HLA4DigitA.get(K);
            pLikelihood4 = (Double) HLA4DigitP.get(K);
            f1 = (Double) HLA4DigitF1.get(K);
            f2 = (Double) HLA4DigitF2.get(K);
            count1 = Double.valueOf(AlleleCount.get(s[0]).toString());
            count2 = Double.valueOf(AlleleCount.get(s[1]).toString());
            locusCount = Double.valueOf(LocusCount.get(s1[0]).toString());
            if (s[0].equals(s[1])){
                accountedFor = count1 / locusCount;
            }else{
                accountedFor = (count1 + count2) / locusCount;
            }
            if (prob > 0.1){
                out.printf("%s\t%s\t%s\t%.1f\t%.1f\t%.2f\t%.2f\t%.1f\t%.2f\t%.0f\t%.0f\t%.0f\t%.2f",s1[0],name1,name2,aLikelihood4/count,pLikelihood4/count,f1,f2,likelihood/count,prob,count1,count2,locusCount,accountedFor);
                for (int i = 0; i < Populations.length; i++){
                    if (AlleleFrequencies[i].containsKey(s[0])){f1 = Double.valueOf(AlleleFrequencies[i].get(s[0]).toString());}else{f1=.000001;}
                    if (AlleleFrequencies[i].containsKey(s[1])){f2 = Double.valueOf(AlleleFrequencies[i].get(s[1]).toString());}else{f2=.000001;}
                    if (!Double.isInfinite(-1*Math.log10(f1*f2))){out.printf("\t%.2f",Math.log10(f1*f2));}else{out.printf("\t-INF");}
                }
                out.print("\n");
            }
            if (!NO_VERBOSE) {
                out.printf("INFO\t%s\t%s\t%s\t%.1f\t%.1f\t%.2f\t%.2f\t%.1f\t%.2f\t%.0f\t%.0f\t%.0f\t%.2f",s1[0],name1,name2,aLikelihood4/count,pLikelihood4/count,f1,f2,likelihood/count,prob,count1,count2,locusCount,accountedFor);
            }
            for (int i = 0; i < Populations.length; i++){
                if (AlleleFrequencies[i].containsKey(s[0])){f1 = Double.valueOf(AlleleFrequencies[i].get(s[0]).toString());}else{f1=.000001;}
                if (AlleleFrequencies[i].containsKey(s[1])){f2 = Double.valueOf(AlleleFrequencies[i].get(s[1]).toString());}else{f2=.000001;}
                if (!Double.isInfinite(-1*Math.log10(f1*f2))){out.printf("\t%.2f",Math.log10(f1*f2));}else{out.printf("\t-INF");}
            }
            out.print("\n");
        }
    }

    Comparator<Double> valueComparator = new Comparator<Double>() {
        @Override public int compare(Double val1, Double val2) {
            return val1.compareTo(val2);
        }
    };
    
    private Integer[] InitializePolymorphicSites(){
        int HLA_A_start = 30018310, HLA_A_end = 30021211, num_A_positions = HLA_A_end - HLA_A_start + 1;
        int HLA_B_start = 31430239, HLA_B_end = 31432914, num_B_positions = HLA_B_end - HLA_B_start + 1;
        int HLA_C_start = 31344925, HLA_C_end = 31347827, num_C_positions = HLA_C_end - HLA_C_start + 1;
        Integer[] polymorphicSites = new Integer[num_A_positions+num_B_positions+num_C_positions];
        for (int i = 0; i < num_A_positions; i++){
            polymorphicSites[i]=HLA_A_start + i;
        }
        for (int i = 0; i < num_C_positions; i++){
            polymorphicSites[i+num_A_positions]=HLA_C_start + i;
        }
        for (int i = 0; i < num_B_positions; i++){
            polymorphicSites[i+num_A_positions+num_C_positions]=HLA_B_start + i;
        }
        return polymorphicSites;
    }

    private int IndexOf(char c){
        switch(c){
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            //case 'D': return 4;
            default: return -1;
        }
    }


    private boolean IsWithinInterval(int pos){
        boolean isWithinInterval = false;
        for (int i = 0; i < numIntervals; i++){
            if (pos >= intervals[i][0] && pos <= intervals[i][1]){
                isWithinInterval = true;
                break;
            }
        }
        return isWithinInterval;
    }
    
    private void UpdateCorrelation(SAMRecord read){
        //Updates correlation table with SNPs from specific read (for phasing)
        String s = cigarparser.FormatRead(read.getCigarString(), read.getReadString());
        ArrayList<Integer> SNPsInRead = new ArrayList<Integer>();
        ArrayList<Integer> readindex = new ArrayList<Integer>();

        int readstart = read.getAlignmentStart();
        int readend = read.getAlignmentEnd();
        int numPositions = PolymorphicSites.length;
        char c1, c2;
        int a, b, i, j, SNPcount = 0;
        
        //Find all SNPs in read
        for (i = 0; i < numPositions; i++){
            if (PolymorphicSites[i] > readstart && PolymorphicSites[i] < readend){
                SNPnumInRead[i] = SNPcount;
                SNPposInRead[i] = PolymorphicSites[i]-readstart;
                SNPcount++;
            }else{
                SNPnumInRead[i] = -1;
                SNPposInRead[i] = -1;
            }
        }

        //Update correlation table; for each combination of SNP positions
        for (i = 0; i < numPositions; i++){
            if (SNPnumInRead[i] > -1){
                c1 = s.charAt(SNPposInRead[i]);
                if (IndexOf(c1) > -1){
                    for (j = i+1; j < numPositions; j ++){
                        if (SNPnumInRead[j] > -1){
                            c2 = s.charAt(SNPposInRead[j]);
                            if (IndexOf(c2) > -1){
                                a = i*5 + IndexOf(c1);
                                b = j*5 + IndexOf(c2);

                                numObservations[a][b]++;
                                totalObservations[i][j]++;
                                if (DEBUG){
                                    //out.printf("INFO  %s\t%s %s\t[i=%s,j=%s]\t[%s,%s]\t[%s,%s]\n",read.getReadName(),PolymorphicSites[i],PolymorphicSites[j],i,j,c1,c2,a,b);
                                }
                            }
                        }
                    }

                }
            }
        }
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

    private double CalculateAlleleLikelihood(int a1, int a2, String[] HLAalleles, boolean debug){
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
                    if (!NO_VERBOSE || debug){
                        c1 = read1.charAt(pos-start1);
                        c2 = read2.charAt(pos-start2);
                        out.printf("INFO: DEBUG %s\t%s\t%s\t%s\t%s\t%s\t%.2f\n",HLAnames[a1],HLAnames[a2],pos,c1,c2,index,likelihood);
                    }
                }
            }
        }
        return likelihood;
    }

    private double CalculatePhaseLikelihood(int alleleIndex1, int alleleIndex2, boolean PRINTDEBUG, boolean SINGLEALLELE){
        //calculate the likelihood that the particular combination of alleles satisfies the phase count data
        double likelihood = 0, prob = 0;
        int readstart1 = HLAstartpos[alleleIndex1]; int readend1 = HLAstoppos[alleleIndex1];
        int readstart2 = HLAstartpos[alleleIndex2]; int readend2 = HLAstoppos[alleleIndex2];
        int combinedstart = Math.max(readstart1,readstart2);
        int combinedstop  = Math.min(readend1,readend2);
        
        int numPositions = PolymorphicSites.length, SNPcount = 0;
        int i, j, a1, a2, b1, b2;
        char c11, c12, c21, c22;
        int numInPhase = 0, numOutOfPhase = 0;
        double sumInPhase = 0.0, sumObservations = 0.0;


        //Find all SNPs in read
        for (i = 0; i < numPositions; i++){

            if (PolymorphicSites[i] > combinedstart && PolymorphicSites[i] < combinedstop  ){ // && IsWithinInterval(PolymorphicSites[i])
                if (PRINTDEBUG){
                    out.printf("DEBUG\t%s\t%s\n",PolymorphicSites[i],SNPcount);
                }
                SNPnumInRead[i] = SNPcount;
                SNPposInRead[i] = PolymorphicSites[i]-combinedstart;
                SNPcount++;
            }else{
                SNPnumInRead[i] = -1;
                SNPposInRead[i] = -1;
            }
        }
        String s1 = HLAreads[alleleIndex1];
        String s2 = HLAreads[alleleIndex2];
        if (PRINTDEBUG){
            out.printf("DEBUG  %s SNPs found in %s and %s between %s and %s\n",SNPcount,HLAnames[alleleIndex1], HLAnames[alleleIndex2],combinedstart,combinedstop);
        }
        //Iterate through every pairwise combination of SNPs, and update likelihood for the allele combination
        for (i = 0; i < numPositions; i++){
            if (SNPnumInRead[i] > -1){
                c11 = s1.charAt(SNPposInRead[i]);
                c21 = s2.charAt(SNPposInRead[i]);
                if (IndexOf(c11) > -1 && IndexOf(c21) > -1){
                    for (j = i+1; j < numPositions; j ++){
                        if (SNPnumInRead[j] > -1 && totalObservations[i][j] > 0){
                            c12 = s1.charAt(SNPposInRead[j]);
                            c22 = s2.charAt(SNPposInRead[j]);
                            if (IndexOf(c12) > -1 && IndexOf(c22) > -1){
                                a1 = i*5 + IndexOf(c11);
                                b1 = j*5 + IndexOf(c12);
                                a2 = i*5 + IndexOf(c21);
                                b2 = j*5 + IndexOf(c22);
                                //check if the two alleles are identical at the chosen 2 locations
                                if ((c11 == c21) && (c12 == c22)){
                                    numInPhase = numObservations[a1][b1];
                                }else{
                                    numInPhase = numObservations[a1][b1] + numObservations[a2][b2];
                                }
                                numOutOfPhase = totalObservations[i][j] - numInPhase;
                                sumInPhase += (double) numInPhase;
                                sumObservations += (double) totalObservations[i][j];
                                if (SINGLEALLELE){
                                    likelihood = sumInPhase / sumObservations;
                                }else{
                                    likelihood += numInPhase * L_correct + numOutOfPhase * L_err;
                                }
                                //prob = Math.max((double) numInPhase / (double) totalObservations[i][j], 0.0001);
                                //likelihood += Math.log10(prob);
                                //likelihood = Math.max(Math.log10(sumInPhase / sumObservations),-10);
                            
                                if (PRINTDEBUG){
                                    out.printf("DEBUG  %s %s %s[%s%s] %s[%s%s]\t[%s,%s]\t[%s,%s] [%s,%s]\t%s / %s\t%s / %s\t    %.2f\n",HLAnames[alleleIndex1],HLAnames[alleleIndex2],PolymorphicSites[i],c11,c21,PolymorphicSites[j],c12,c22, i,j,a1,b1,a2,b2,numInPhase,totalObservations[i][j],sumInPhase,sumObservations,likelihood);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
        return likelihood;
    }

    private void ExtraCode(){
        String name1, name2;
        //Pre-process homozygous combinations to determine top possible alleles (for efficiency)
        Hashtable Alleles2Digit = new Hashtable();
        Hashtable Phase2Digit = new Hashtable();
        Hashtable Count2Digit = new Hashtable();
        
        Hashtable AllelesAtLocus = new Hashtable();
        ArrayList Loci = new ArrayList<String>();
        double[] AlleleLikelihoods2 = new double[HLAnames.length];
        double[] PhaseLikelihoods2 = new double[HLAnames.length];
        for (int i = 0; i < HLAnames.length; i++){
            name1 = HLAnames[i].substring(4);
            String [] n1 = name1.split("\\*");
            AlleleLikelihoods2[i] = CalculateAlleleLikelihood(i,i,HLAreads,false);
            PhaseLikelihoods2[i] = CalculatePhaseLikelihood(i,i,false,true);
            if (AlleleLikelihoods2[i] < 0){
                name2 = n1[0] + "*" + n1[1].substring(0, 4);
                if (!Loci.contains(n1[0])){
                    Loci.add(n1[0]);
                    MaxLikelihoods.put(n1[0], 0.0);
                    AllelesAtLocus.put(n1[0], 1);
                }else{
                    AllelesAtLocus.put(n1[0], 1+(Integer)AllelesAtLocus.get(n1[0]));
                }
                if (!Alleles2Digit.containsKey(name2)){
                    Alleles2Digit.put(name2, AlleleLikelihoods2[i]);
                    Phase2Digit.put(name2, PhaseLikelihoods2[i]);
                    Count2Digit.put(name2, 1.0);
                }else {
                    if (AlleleLikelihoods2[i] > (Double) Alleles2Digit.get(name2)){
                        Alleles2Digit.put(name2, AlleleLikelihoods2[i]);
                    }
                    if (PhaseLikelihoods2[i] > (Double) Phase2Digit.get(name2)){
                        Phase2Digit.put(name2, PhaseLikelihoods2[i]);
                    }
                    Count2Digit.put(name2,1.0+(Double)Count2Digit.get(name2));
                }
            }
        }

        //Sort alleles at 2 digit resolution for each locus

        for (int i = 0; i < Loci.size(); i++){
            Enumeration k = Alleles2Digit.keys();
            Hashtable AllelesAtLoci = new Hashtable();
            HashMap map = new HashMap();
            int numalleles = 0;
            //find alleles at the locus
            while( k.hasMoreElements() ){
                name1 = k.nextElement().toString();
                String [] n1 = name1.split("\\*");
                if (Loci.get(i).equals(n1[0])){
                    numalleles++;
                    map.put(name1,-1 * (Double) Alleles2Digit.get(name1));
                    AllelesAtLoci.put(-1 * (Double) Alleles2Digit.get(name1), name1);
                    //out.printf("%s\t%.2f\n",name1,-1 * (Double) Alleles2Digit.get(name1));
                }

            }

            //Sort alleles at locus, mark top six 2-digit classes for deep search
            List<Map.Entry<String, Double>> entries = new ArrayList<Entry<String, Double>>(map.entrySet());
            Collections.sort(entries, new Comparator<Entry<String, Double>>() {
                public int compare(Entry<String, Double> e1, Entry<String, Double> e2) {
                    return e1.getValue().compareTo(e2.getValue());
                }
            });
            int num = 1;
            for (Map.Entry<String, Double> entry : entries) {
                if (num <= Math.max(5,entries.size()/8)){
                    AllelesToSearch.add(entry.getKey());
                    if (!NO_VERBOSE) {
                        out.printf("INFO\t%s\t%.2f\t%.2f\n",entry.getKey(),entry.getValue(),Phase2Digit.get(entry.getKey()));
                    }
                    num++;
                }else if (!NO_VERBOSE) {
                    if (!AllelesToSearch.contains(entry.getKey())){
                        out.printf("INFO\t%s\t%.2f\t%.2f\tNotSearched\n",entry.getKey(),entry.getValue(),Phase2Digit.get(entry.getKey()));
                    }else{
                        out.printf("INFO\t%s\t%.2f\t%.2f\n",entry.getKey(),entry.getValue(),Phase2Digit.get(entry.getKey()));
                    }
                }
            }

            if (!NO_VERBOSE) {
                out.printf("INFO\n");
            }
        }
    }
}
