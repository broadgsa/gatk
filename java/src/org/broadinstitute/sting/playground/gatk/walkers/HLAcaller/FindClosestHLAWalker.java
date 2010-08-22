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
import java.io.PrintStream;

/**
 * Finds the most similar HLA allele for each read (helps detect misalignments). Usage: java -jar GenomeAnalysisTK.jar -T FindClosestHLA -I INPUT.bam -R /broad/1KG/reference/human_b36_both.fasta -L INPUT.interval | grep -v INFO | sort -k1 > OUTPUT
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class FindClosestHLAWalker extends ReadWalker<Integer, Integer> {
    @Output
    protected PrintStream out;

    @Argument(fullName = "debugRead", shortName = "debugRead", doc = "Print match score for read", required = false)
    public String debugRead = "";

    @Argument(fullName = "findFirst", shortName = "findFirst", doc = "For each read, stop when first HLA allele is found with concordance = 1", required = false)
    public boolean findFirst = false;

    @Argument(fullName = "DEBUG", shortName = "DEBUG", doc = "Debug walker", required = false)
    public boolean DEBUG = false;
    
    @Argument(fullName = "debugAllele", shortName = "debugAllele", doc = "Print match score for allele", required = false)
    public String debugAllele = "";
    
    @Argument(fullName = "useInterval", shortName = "useInterval", doc = "Use only these intervals", required = false)
    public String intervalFile = "";

    @Argument(fullName = "dictionary", shortName = "dictionary", doc = "bam file of HLA ditionary", required = false)
    public String HLAdictionaryFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.sam";

    @Argument(fullName = "onlyfrequent", shortName = "onlyfrequent", doc = "Only consider alleles with frequency > 0.0001", required = false)
    public boolean ONLYFREQUENT = false;
    
    @Argument(fullName = "HLAdictionary", shortName = "HLAdictionary", doc = "HLA dictionary file", required = true)
    public String HLAdatabaseFile = "HLA_DICTIONARY.txt";

    @Argument(fullName = "PolymorphicSites", shortName = "PolymorphicSites", doc = "file containing polymorphic sites within the HLA", required = true)
    public String PolymorphicSitesFile = "HLA_POLYMORPHIC_SITES.txt";
    
    HLAFileReader HLADictionaryReader = new HLAFileReader();    

    boolean DatabaseLoaded = false;
    ArrayList<String> ClosestAlleles = new ArrayList<String>();

    String[] HLAnames, HLAreads;
    Integer[] HLAstartpos, HLAstoppos, PolymorphicSites,NonPolymorphicSites;
    double[] SingleAlleleFrequencies;

    double[] nummatched, concordance, numcompared;
    int numHLAlleles = 0;
    int minstartpos = 0;
    int maxstoppos = 0;
    int numpolymorphicsites = 0, numnonpolymorphicsites = 0, pos =0;

    Hashtable AlleleFrequencies = new Hashtable();
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;
    CigarParser formatter = new CigarParser();
    int [][] intervals; int numIntervals;

    public Integer reduceInit() { 
        if (!DatabaseLoaded){
            DatabaseLoaded = true;

            //Load HLA dictionary
            out.printf("INFO  Loading HLA dictionary ... ");

            HLADictionaryReader.ReadFile(HLAdatabaseFile);
            HLAreads = HLADictionaryReader.GetSequences();
            HLAnames = HLADictionaryReader.GetNames();
            HLAstartpos = HLADictionaryReader.GetStartPositions();
            HLAstoppos = HLADictionaryReader.GetStopPositions();
            minstartpos = HLADictionaryReader.GetMinStartPos();
            maxstoppos = HLADictionaryReader.GetMaxStopPos();
            
            out.printf("Done! %s HLA alleles loaded.\n",HLAreads.length);

            nummatched = new double[HLAreads.length];
            concordance = new double[HLAreads.length];
            numcompared = new double[HLAreads.length];

            //Load list of polymorphic sites
            PolymorphicSitesFileReader siteFileReader = new PolymorphicSitesFileReader();
            siteFileReader.ReadFile(PolymorphicSitesFile);
            PolymorphicSites = siteFileReader.GetPolymorphicSites();
            NonPolymorphicSites = siteFileReader.GetNonPolymorphicSites();
            numpolymorphicsites = PolymorphicSites.length;
            numnonpolymorphicsites = NonPolymorphicSites.length;

            if (!intervalFile.equals("")){
                TextFileReader fileReader = new TextFileReader();
                fileReader.ReadFile(intervalFile);
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
            
            out.printf("INFO  %s polymorphic and %s non-polymorphic sites found in HLA dictionary\n",numpolymorphicsites,numnonpolymorphicsites);
            out.printf("INFO  Comparing reads to database ...\n");

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
    }

    private double CalculateConcordance(SAMRecord read){
        int readstart = read.getAlignmentStart();
        int readstop = read.getAlignmentEnd();
        char c1, c2;
        double maxConcordance = 0.0, freq = 0.0, minFreq = 0.0;
        String s1 = formatter.FormatRead(read.getCigarString(), read.getReadString());
        String s2;
        int allelestart, allelestop;

        if (ONLYFREQUENT){
            minFreq = 0.0001;
        }

        for (int i = 0; i < HLAreads.length; i++){
            nummatched[i] = 0; concordance[i] = 0; numcompared[i] = 0;
            freq = GetAlleleFrequency(HLAnames[i]);
            //Get concordance between read and specific allele
            if (readstart <= HLAstoppos[i] && readstop >= HLAstartpos[i] && freq > minFreq){
                s2 = HLAreads[i];
                
                allelestart = HLAstartpos[i];
                allelestop = HLAstoppos[i];

                //Polymorphic sites: always increment denominator, increment numerator when bases are concordant
                for (int j = 0; j < numpolymorphicsites; j++){
                    pos = PolymorphicSites[j];
                    if (DEBUG == true){
                        out.printf("DEBUG\tPOS\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(),HLAnames[i],pos,allelestart,allelestop,IsWithin(pos,readstart,readstop), IsWithin(pos,allelestart,allelestop),IsWithinInterval(pos));
                    }
                    if (pos >= readstart && pos <= readstop && pos >= allelestart && pos <= allelestop && IsWithinInterval(pos)){
                        c1 = s1.charAt(pos-readstart);
                        c2 = s2.charAt(pos-allelestart);
                        if (c1 != 'D' && c2 != 'D'){//allow for deletions (sequencing errors)
                            numcompared[i]++;
                            if (c1 == c2){
                                nummatched[i]++;
                            }else{
                                if (debugRead.equals(read.getReadName()) && debugAllele.equals(HLAnames[i])){
                                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(), HLAnames[i], j, pos,c1,c2);
                                }
                            }
                        }
                    }
                }

                //Non-polymorphic sites: increment denominator only when bases are discordant
                if (numcompared[i] > 0){
                    for (int j = 0; j < numnonpolymorphicsites; j++){
                        pos = NonPolymorphicSites[j];
                        if (pos >= readstart && pos <= readstop && pos >= allelestart && pos <= allelestop && IsWithinInterval(pos)){
                            c1 = s1.charAt(pos-readstart);
                            c2 = s2.charAt(pos-allelestart);
                            if (c1 != c2 && c1 != 'D' && c2 != 'D'){//allow for deletions (sequencing errors)
                                numcompared[i]++;
                                if (debugRead.equals(read.getReadName()) && debugAllele.equals(HLAnames[i])){
                                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(), HLAnames[i], j, c1,c2);
                                }
                            }
                        }
                    }
                }
            
                //Update concordance array
                concordance[i]=nummatched[i]/numcompared[i];
                if (concordance[i] > maxConcordance){maxConcordance = concordance[i];}
                if (DEBUG == true){
                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(),HLAnames[i],concordance[i],numcompared[i],numcompared[i]-nummatched[i]);
                }
                if (debugRead.equals(read.getReadName()) && debugAllele.equals(HLAnames[i])){
                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(),HLAnames[i],concordance[i],numcompared[i],numcompared[i]-nummatched[i]);
                }
                if (findFirst && (concordance[i] == 1)){
                    break;
                }
            }
        
        }

        return maxConcordance;
    }

    private double FindMaxAlleleFrequency(double maxConcordance){
        //finds the max frequency of the alleles that share the maximum concordance with the read of interest
        double freq, maxFreq = 0.0;
        for (int i = 0; i < HLAreads.length; i++){
            if (concordance[i] == maxConcordance && maxConcordance > 0){
                freq = GetAlleleFrequency(HLAnames[i]);
                if (freq > maxFreq){maxFreq = freq;}
            }
        }
        return maxFreq;
    }

    private boolean IsWithin(int pos, int start, int stop){
        return pos >= start && pos <= stop;
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
    
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        //Calculate concordance for this read and all overlapping reads
        if (DEBUG == true){
            out.printf("%s\t%s\n",read.getReadName(),read.getMappingQuality());
        }
        if (read.getMappingQuality() > 0 || DEBUG == true){
            double maxConcordance = CalculateConcordance(read);
            String stats = "", topAlleles = "";
            if (maxConcordance > 0 || DEBUG == true){
                String readname = read.getReadName(), allelename = ""; double freq;
                //For input bam files that contain HLA alleles, find and print allele frequency
                out.printf("%s\t%s-%s", readname,read.getAlignmentStart(),read.getAlignmentEnd());

                //Print concordance statistics between this read and the most similar HLA allele(s)

                for (int i = 0; i < HLAreads.length; i++){
                    if (concordance[i] == maxConcordance){
                        //freq = GetAlleleFrequency(HLAnames[i]);
                        if (topAlleles.equals("")){
                            topAlleles = HLAnames[i];
                        }else{
                            topAlleles = topAlleles + "," + HLAnames[i];
                        }
                        stats = String.format("%.1f\t%.3f\t%.0f\t%.0f",1.0,concordance[i],numcompared[i],numcompared[i]-nummatched[i]);
                        
                    }
                }
                out.printf("\t%s\t%s\t%s\n",stats,topAlleles,maxConcordance);
            }
        }
        return 1;
    }

    private double GetAlleleFrequency(String allelename){
        double frequency = 0.0;
        //Truncate names to 4-digit "A*0101" format
        if (allelename.length() >= 10){
            allelename=allelename.substring(4,10);
        }else{
            allelename=allelename.substring(4);
        }
        if (AlleleFrequencies.containsKey(allelename)){
            frequency = Double.parseDouble((String) AlleleFrequencies.get(allelename).toString());
        }else{
            frequency=0.0001;
        }
        return frequency;
    }


    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
        // Double check traversal result to make count is the same.
        // TODO: Is this check necessary?
        out.println("[REDUCE RESULT] Traversal result is: " + result);
    }    
}

