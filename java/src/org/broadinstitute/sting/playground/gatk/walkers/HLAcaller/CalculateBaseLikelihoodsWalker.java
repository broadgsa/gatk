/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
/**
 *
 * @author shermanjia
 */
public class CalculateBaseLikelihoodsWalker extends LocusWalker<Integer, Pair<Long, Long>>{
    @Argument(fullName = "debugHLA", shortName = "debugHLA", doc = "Print debug", required = false)
    public boolean DEBUG = false;
    @Argument(fullName = "debugAlleles", shortName = "debugAlleles", doc = "Print likelihood scores for these alleles", required = false)
    public String inputAlleles = "";
    @Argument(fullName = "ethnicity", shortName = "ethnicity", doc = "Use allele frequencies for this ethnic group", required = false)
    public String ethnicity = "CaucasianUSA";
    @Argument(fullName = "filter", shortName = "filter", doc = "file containing reads to exclude", required = false)
    public String filterFile = "";
    @Argument(fullName = "minAllowedMismatches", shortName = "minAllowedMismatches", doc = "Min number of mismatches tolerated per read (default 7)", required = false)
    public int MINALLOWEDMISMATCHES = 7;

    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.imputed.4digit.sam";
    String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";
    String BlackAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_BlackUSA.freq";
    String AlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";

    SAMFileReader HLADictionaryReader = new SAMFileReader();
    String[] HLAreads, HLAnames;
    Integer[] HLAstartpos, HLAstoppos;
    Hashtable AlleleFrequencies;

    int[][] LOD, LikelihoodScores;
    ArrayList<String> ReadsToDiscard = new ArrayList<String>();

    ArrayList<SAMRecord> AllReads = new ArrayList<SAMRecord>();
    ArrayList<String> AllReadNames = new ArrayList<String>();

    boolean HLAdataLoaded = false;


    //Loads HLA dictionary, allele frequencies, and reads to filter
    public Pair<Long, Long> reduceInit() {
        if (!HLAdataLoaded){
            HLAdataLoaded = true;
            
            out.printf("INFO  Reading HLA database ... ");
            HLADictionaryReader.ReadFile(HLAdatabaseFile);
            HLAreads = HLADictionaryReader.GetReads();
            HLAnames = HLADictionaryReader.GetReadNames();
            HLAstartpos = HLADictionaryReader.GetStartPositions();
            HLAstoppos = HLADictionaryReader.GetStopPositions();
            InitializeVariables(HLAreads.length);
            out.printf("Done! %s HLA alleles loaded.\n",HLAreads.length);



            if (!ethnicity.equals("CaucasianUSA")){
                AlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_" + ethnicity + ".freq";
            }
            out.printf("INFO  Reading HLA allele frequencies ... ");
            FrequencyFileReader HLAfreqReader = new FrequencyFileReader();
            HLAfreqReader.ReadFile(AlleleFrequencyFile);
            AlleleFrequencies = HLAfreqReader.GetAlleleFrequencies();
            out.printf("Done! Frequencies for %s HLA alleles loaded.\n",AlleleFrequencies.size());

            
            if (!filterFile.equals("")){
                out.printf("INFO  Reading properties file ... ");
                SimilarityFileReader similarityReader = new SimilarityFileReader();
                similarityReader.ReadFile(filterFile,MINALLOWEDMISMATCHES);
                ReadsToDiscard = similarityReader.GetReadsToDiscard();
                out.printf("Done! Found %s misaligned reads to discard.\n",ReadsToDiscard.size());
                for (int i = 0; i < ReadsToDiscard.size(); i++){
                    out.printf("MISALIGNED %s\n", ReadsToDiscard.get(i).toString());
                }
            }
        }
        return new Pair<Long,Long>(0l,0l);
    }



    private void InitializeVariables(int n){
        LOD = new int[n][n];
        LikelihoodScores = new int[n][n];
        for (int i = 0; i < n; i++){

            for (int j = 0; j <n; j++){
                LOD[i][j]=0;
                LikelihoodScores[i][j]=0;
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        List<SAMRecord> reads = context.getReads();
        if(reads.size() > 0 ) {
            List<Integer> offsets = context.getOffsets();

            int numAs = 0, numCs = 0, numGs = 0, numTs = 0;
            //if (DEBUG){
                out.printf("%s\t%s\t", context.getLocation(),ref.getBase());
            //}

            //Calculate posterior probabilities
            GenotypeLikelihoods G = new GenotypeLikelihoods(BaseMismatchModel.THREE_STATE);
            SAMRecord read; int offset; char base; byte qual; int mapquality; String readname;

            //Check for bad bases and ensure mapping quality
            for (int i = 0; i < reads.size(); i++) {
                read = reads.get(i);
                offset = offsets.get(i);
                base = read.getReadString().charAt(offset);
                qual = read.getBaseQualities()[offset];
                //mapquality = read.getMappingQuality();
                if (!ReadsToDiscard.contains(read.getReadName()) && BaseUtils.simpleBaseToBaseIndex(base) != -1) {
                    
                    //consider base in likelihood calculations if it looks good and has high mapping score
                    G.add(base, qual, read, offset);
                    //if (DEBUG){
                        if (base == 'A'){numAs++;}
                        else if (base == 'C'){numCs++;}
                        else if (base == 'T'){numTs++;}
                        else if (base == 'G'){numGs++;}
                    //}

                }
            }
            //if (DEBUG) {
                out.printf("A[%s]C[%s]T[%s]G[%s]",numAs,numCs,numTs,numGs);
                for ( DiploidGenotype g : DiploidGenotype.values() ) {
                    out.printf("\t%.2f",G.getLikelihood(g));
                }
                out.printf("\n");
            //}
            
        }
        return context.getReads().size();
    }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = value.longValue() + sum.getFirst();
        long right = sum.getSecond() + 1l;
        return new Pair<Long,Long>(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {
        
    }
}
