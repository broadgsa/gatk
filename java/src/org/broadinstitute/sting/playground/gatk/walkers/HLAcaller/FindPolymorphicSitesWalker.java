package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;

import java.util.ArrayList;
import java.util.Hashtable;
/**
 * Finds polymorphic sites in the HLA dictionary. Usage: java -jar GenomeAnalysisTK.jar -T FindPolymorphicSites -I HLA_DICTIONARY.bam -R /broad/1KG/reference/human_b36_both.fasta -L INPUT.interval -findFirst | grep -v INFO | sort -k1 > OUTPUT
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class FindPolymorphicSitesWalker extends ReadWalker<Integer, Integer> {
    @Argument(fullName = "debugRead", shortName = "debugRead", doc = "Print match score for read", required = false)
    public String debugRead = "";

    @Argument(fullName = "findFirst", shortName = "findFirst", doc = "For each read, stop when first HLA allele is found with concordance = 1", required = false)
    public boolean findFirst = false;

    @Argument(fullName = "debugAllele", shortName = "debugAllele", doc = "Print match score for allele", required = false)
    public String debugAllele = "";

    @Argument(fullName = "ethnicity", shortName = "ethnicity", doc = "Use allele frequencies for this ethnic group", required = false)
    public String ethnicity = "Caucasian";

    @Argument(fullName = "onlyfrequent", shortName = "onlyfrequent", doc = "Only consider alleles with frequency > 0.0001", required = false)
    public boolean ONLYFREQUENT = false;

    //String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.imputed.4digit.sam";
    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_DICTIONARY.sam";
    
    SAMFileReader HLADictionaryReader = new SAMFileReader();

    //String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";
    String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_Caucasians.freq";
    String BlackAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_BlackUSA.freq";
    String AlleleFrequencyFile;
    String UniqueAllelesFile               = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/UniqueAlleles";

    String PolymorphicSitesFile = "/humgen/gsa-scr1/GSA/sjia/Sting/HLA.polymorphic.sites";

    boolean DatabaseLoaded = false;
    boolean DEBUG = false;

    String[] HLAnames, HLAreads;
    Integer[] HLAstartpos, HLAstoppos, PolymorphicSites,NonPolymorphicSites;
    double[] SingleAlleleFrequencies;

    double[] nummatched, concordance, numcompared;
    int numHLAlleles = 0;
    int minstartpos = 0;
    int maxstoppos = 0;

    int HLA_A_start = 30018310;
    int HLA_A_end = 30021211;

    Hashtable AlleleFrequencies = new Hashtable();
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;
    CigarParser formatter = new CigarParser();

    public Integer reduceInit() {
        if (!DatabaseLoaded){
            DatabaseLoaded = true;

            //Load HLA dictionary
            out.printf("INFO  Loading HLA dictionary ... ");

            HLADictionaryReader.ReadFile(HLAdatabaseFile);
            HLAreads = HLADictionaryReader.GetReads();
            HLAnames = HLADictionaryReader.GetReadNames();
            HLAstartpos = HLADictionaryReader.GetStartPositions();
            HLAstoppos = HLADictionaryReader.GetStopPositions();
            minstartpos = HLADictionaryReader.GetMinStartPos();
            maxstoppos = HLADictionaryReader.GetMaxStopPos();

            out.printf("Done! %s HLA alleles loaded.\n",HLAreads.length);

            nummatched = new double[HLAreads.length];
            concordance = new double[HLAreads.length];
            numcompared = new double[HLAreads.length];

            //Read allele frequencies
            if (ethnicity.equals("Black")){
                AlleleFrequencyFile = BlackAlleleFrequencyFile;
            }else{
                AlleleFrequencyFile = CaucasianAlleleFrequencyFile;
            }
            out.printf("INFO  Reading HLA allele frequencies ... ");
            FrequencyFileReader HLAfreqReader = new FrequencyFileReader();
            HLAfreqReader.ReadFile(AlleleFrequencyFile,UniqueAllelesFile);
            AlleleFrequencies = HLAfreqReader.GetAlleleFrequencies();
            out.printf("Done! Frequencies for %s HLA alleles loaded.\n",AlleleFrequencies.size());

            FindPolymorphicSites(minstartpos,maxstoppos);

            out.printf("INFO  %s polymorphic and %s non-polymorphic sites found in HLA dictionary\n",PolymorphicSites.length,NonPolymorphicSites.length);
            out.printf("INFO  Comparing reads to database ...\n");

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
    }

    private void FindPolymorphicSites(int start, int stop){
        boolean initialized, polymorphic, examined;
        char c = ' ';
        ArrayList<Integer> polymorphicsites = new ArrayList<Integer>();
        ArrayList<Integer> nonpolymorphicsites = new ArrayList<Integer>();
        //Find polymorphic sites in dictionary
        for (int pos = start; pos <= stop; pos++){
            initialized = false; polymorphic = false; examined = false;
            //look across all alleles at specific position to see if it is polymorphic
            for (int i = 0; i < HLAreads.length; i++){
                if (pos >= HLAstartpos[i] && pos <= HLAstoppos[i]){
                    if (!initialized){
                        c = HLAreads[i].charAt(pos-HLAstartpos[i]);
                        initialized = true;
                        examined = true;
                    }
                    if (HLAreads[i].charAt(pos-HLAstartpos[i]) != c){
                        polymorphicsites.add(pos);
                        out.printf("POLYMORPHIC\t6\t%s\n", pos);
                        polymorphic = true;
                        break;
                    }
                }

            }
            if (!polymorphic && examined){
                nonpolymorphicsites.add(pos);
                out.printf("CONSERVED\t6\t%s\n", pos);
            }
        }
        PolymorphicSites = polymorphicsites.toArray(new Integer[polymorphicsites.size()]);
        NonPolymorphicSites = nonpolymorphicsites.toArray(new Integer[nonpolymorphicsites.size()]);
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        //Calculate concordance for this read and all overlapping reads
        return 1;
    }


    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}

