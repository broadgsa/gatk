/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ÓSoftwareÓ), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ÓAS ISÓ, WITHOUT WARRANTY OF ANY KIND,
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
import java.util.Enumeration;
/**
 * Creates a ped file of SNPs and amino acids coded as SNPs given an input ped file with 4-digit HLA alleles. Usage: java -jar GenomeAnalysisTK.jar -T CreatePedFile --allelesFile INPUT.ped -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-sc\
r1/GSA/sjia/454_HLA/HLA/HLA.combined.4digitUnique.bam > OUTPUT.log
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CreatePedFileWalker extends ReadWalker<Integer, Integer> {
    @Argument(fullName = "allelesFile", shortName = "allelesFile", doc = "Create ped file for HLA alleles named in this file", required = true)
    public String alleleNamesFile = "";

    @Argument(fullName = "pedIntervals", shortName = "pedIntervals", doc = "Create genotypes in these intervals", required = false)
    public String pedIntervalsFile = "";

    @Argument(fullName = "HLAexonIntervals", shortName = "HLAexonIntervals", doc = "HLA exonic intervals", required = false)
    public String exonIntervalsFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_EXON_POSITIONS.txt";

    @Argument(fullName = "DNAcode", shortName = "DNAcode", doc = "Amino acid codes", required = false)
    public String dnaCodesFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/DNA_CODE.txt";

    String[] HLAnames, HLAreads, inputFileContents;
    Integer[] HLAstartpos, HLAstoppos;
    ArrayList<String> HLAnamesAL, HLAreadsAL;
    ArrayList<Integer> HLAstartposAL, HLAstopposAL;
    int[][] intervals; String[][] exonIntervals;
    int numIntervals;

    ReadCigarFormatter formatter = new ReadCigarFormatter();
    char c;
    boolean DEBUG = false;
    boolean FilesLoaded = false;
    int HLA_A_start    = 30018310, HLA_A_end    = 30021211;
    int HLA_C_start    = 31344925, HLA_C_end    = 31347827;
    int HLA_B_start    = 31430239, HLA_B_end    = 31432914;
    int HLA_DRB1_start = 32654846, HLA_DRB1_end = 32665497;
    int HLA_DQA1_start = 32713214, HLA_DQA1_end = 32718519;
    int HLA_DQB1_start = 32735991, HLA_DQB1_end = 32742362;
    int HLA_DPA1_start = 33144405, HLA_DPA1_end = 33149325;
    int HLA_DPB1_start = 33151797, HLA_DPB1_end = 33161993;


    String[] SNPnames;
    String SNPname;
    int start, end;
    Integer I;

    Hashtable indexer = new Hashtable();
    Hashtable DNAcode = new Hashtable();

    public Integer reduceInit() {
        if (!FilesLoaded){
            FilesLoaded = true;
            HLAnamesAL = new ArrayList<String>();
            HLAreadsAL = new ArrayList<String>();
            HLAstartposAL = new ArrayList<Integer>();
            HLAstopposAL = new ArrayList<Integer>();

            TextFileReader fileReader = new TextFileReader();
            fileReader.ReadFile(alleleNamesFile);
            inputFileContents = fileReader.GetLines();


            //Determine intervals
            if (!pedIntervalsFile.equals("")){
                fileReader = new TextFileReader();
                fileReader.ReadFile(pedIntervalsFile);
                String[] lines = fileReader.GetLines();
                intervals = new int[lines.length][2];
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split(":");
                    String[] intervalPieces = s[0].split("-");
                    intervals[i][0] = Integer.valueOf(intervalPieces[0]);
                    intervals[i][1] = Integer.valueOf(intervalPieces[1]);
                }
                numIntervals = intervals.length;
                for (int i = 0; i < numIntervals; i++){
                    out.printf("INFO  Interval %s: %s-%s\n",i+1,intervals[i][0],intervals[i][1]);
                }
            }

            //load HLA exonic intervals
            if (!exonIntervalsFile.equals("")){
                fileReader = new TextFileReader();
                fileReader.ReadFile(exonIntervalsFile);
                String[] lines = fileReader.GetLines();
                exonIntervals = new String[lines.length][5];
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split("\t");
                    String[] intervalPieces = s[1].split("-");
                    exonIntervals[i][1] = intervalPieces[0];
                    exonIntervals[i][2] = intervalPieces[1];
                    exonIntervals[i][0] = s[0]; // Locus
                    exonIntervals[i][3] = s[2]; // Exon number
                    exonIntervals[i][4] = s[3]; // +/- strand
                }
                numIntervals = exonIntervals.length;
                for (int i = 0; i < numIntervals; i++){
                    out.printf("INFO  HLA-%s %s (%s): %s-%s\n",exonIntervals[i][0],exonIntervals[i][3],exonIntervals[i][4],exonIntervals[i][1],exonIntervals[i][2]);
                }
            }

            //load amino-acid coding DNA triplets
            if (!dnaCodesFile.equals("")){
                fileReader = new TextFileReader();
                fileReader.ReadFile(dnaCodesFile);
                String[] lines = fileReader.GetLines();
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split("\t");
                    DNAcode.put(s[0],s[1]);
                }

                Enumeration e = DNAcode.keys();
                while( e.hasMoreElements() ){
                    String key = e.nextElement().toString();
                    out.printf("INFO %s encodes %s\n",key,DNAcode.get(key));
                }
            }
        }
        return 0;
    }

    private String[][] GetExonIntervals(String locus, boolean isForwardStrand){
        int numExons = 0; int exonNum;
        for (int i = 0; i < exonIntervals.length; i++){
            if (exonIntervals[i][0].equals(locus)){
                numExons++;
            }
        }
        String[][] ExonIntervals = new String[numExons][5];
        if (isForwardStrand){exonNum = 1;}else{exonNum = ExonIntervals.length;}
        for (int i = 0; i < exonIntervals.length; i++){
            if (exonIntervals[i][0].equals(locus)){
                ExonIntervals[exonNum-1]=exonIntervals[i];
                if (isForwardStrand){
                    exonNum++;
                }else{
                    exonNum--;
                }
            }
        }
        return ExonIntervals;
    }

    private int BaseCharToInt(char c){
        switch(c){
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 3;
            case 'T': return 4;
            default: return -1;
        }
    }

    private char Complement(char c){
        switch(c){
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            default: return '0';
        }
    }

    private char GetAminoAcid(String codon){
        if (DNAcode.containsKey(codon)){
            return DNAcode.get(codon).toString().charAt(0);
        }else{
            return '0';
        }
    }

    public Integer map(char[] ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        HLAnamesAL.add(read.getReadName());
        HLAreadsAL.add(formatter.FormatRead(read.getCigarString(), read.getReadString()));
        HLAstartposAL.add(read.getAlignmentStart());
        HLAstopposAL.add(read.getAlignmentEnd());
        return 1;
    }

    private String PrintGenotypes(String ID, String alleleName1, String alleleName2, int startpos, int stoppos){

        String error = "";
        //prints genotypes for allele1 and allele2 at given interval
        int i1 = GetAlleleIndex(alleleName1);
        int i2 = GetAlleleIndex(alleleName2);
        String s1, s2;
        int start1, start2, stop1, stop2;
        char c1, c2;

        if (i1 > -1){
            s1 = HLAreads[i1];
            start1 = HLAstartpos[i1];
            stop1 = HLAstoppos[i1];
        }else{
            error = error + "INFO  " + alleleName1 + " for " + ID + " not found in HLA dictionary\n";
            s1 = "";
            start1 = -1;
            stop1 = -1;
        }

        if (i2 > -1){
            s2 = HLAreads[i2];
            start2 = HLAstartpos[i2];
            stop2 = HLAstoppos[i2];
        }else{
            error = error + "INFO  " + alleleName2 + " for " + ID + " not found in HLA dictionary\n";
            s2 = "";
            start2 = -1;
            stop2 = -1;
        }
        
        for (int pos = startpos; pos <= stoppos; pos++){
            c1 = GetBase(pos,s1,start1,stop1);
            c2 = GetBase(pos,s2,start2,stop2);
            out.printf("\t%s %s",c1,c2);
        }
        return error;
    }

private String PrintAminoAcids(String ID, String alleleName1, String alleleName2, String[][] ExonIntervals){

        String error = "";
        //prints genotypes for allele1 and allele2 at given interval
        int i1 = GetAlleleIndex(alleleName1);
        int i2 = GetAlleleIndex(alleleName2);
        String s1, s2;
        int start1, start2, stop1, stop2;
        char c1, c2;
        boolean isForwardStrand = false;
        if (ExonIntervals[0][4].equals("+")){isForwardStrand=true;}

        int AAcount=0;
        int baseCount=0;
        String codon1 = ""; String codon2 = "";

        if (i1 > -1){
            s1 = HLAreads[i1];
            start1 = HLAstartpos[i1];
            stop1 = HLAstoppos[i1];
        }else{
            s1 = "";
            start1 = -1;
            stop1 = -1;
            error = error + "INFO  " + alleleName1 + " for " + ID + " not found in HLA dictionary\n";
        }

        if (i2 > -1){
            s2 = HLAreads[i2];
            start2 = HLAstartpos[i2];
            stop2 = HLAstoppos[i2];
        }else{
            s2 = "";
            start2 = -1;
            stop2 = -1;
            error = error + "INFO  " + alleleName2 + " for " + ID + " not found in HLA dictionary\n";
        }

        int i;
        for (int exonNum = 1; exonNum <= ExonIntervals.length; exonNum++){
            if (isForwardStrand){i=exonNum-1;}else{i=ExonIntervals.length-exonNum;}
            int exonStart = Integer.parseInt(ExonIntervals[i][1]);
            int exonStop = Integer.parseInt(ExonIntervals[i][2]);
            for (int pos = exonStart; pos <= exonStop; pos++){
                c1 = GetBase(pos,s1,start1,stop1);
                c2 = GetBase(pos,s2,start2,stop2);
                if (!isForwardStrand){
                    c1 = Complement(c1);
                    c2 = Complement(c2);
                }
                if (baseCount < 3){
                    if (isForwardStrand){
                        codon1 = codon1 + c1;
                        codon2 = codon2 + c2;
                    }else{
                        codon1 = c1 + codon1;
                        codon2 = c2 + codon2;
                    }
                    baseCount++;
                }

                if (baseCount == 3){
                    out.printf("\t%s %s",GetAminoAcid(codon1),GetAminoAcid(codon2));
                    baseCount = 0;
                    AAcount++;
                    codon1 = "";
                    codon2 = "";
                }
            }
        }
        if (baseCount > 0){
            //Print stop or start codon depending on strandedness
            if (isForwardStrand){out.printf("\tO O");}else{out.printf("\tM M");}
        }
        
        return error;
    }

    private char GetBase(int pos, String str, int start, int stop){
        char base;
        if (pos >= start && pos <= stop){
            base = str.charAt(pos-start);
            if (base == 'D'){base = '0';}
        }else{
            base = '0';
        }
        return base;
    }

    private int GetAlleleIndex(String alleleName){
        //Find first allele that matches name, or matches part of name for 2-digit allele
        int i;
        for (i = 0; i < HLAnames.length; i++){
            if (HLAnames[i].indexOf(alleleName) > -1){
                return i;
            }
        }
        return -1;
        
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    private String GetAlleleName(String locus, String sep, String allele){
        if (allele.length() > 1){
            return locus + sep + allele;
        }else{
            return locus + sep + "0000";
        }
    }

    public void onTraversalDone(Integer numreads) {
        HLAnames = HLAnamesAL.toArray(new String[numreads]);
        HLAreads = HLAreadsAL.toArray(new String[numreads]);
        HLAstartpos = HLAstartposAL.toArray(new Integer[numreads]);
        HLAstoppos = HLAstopposAL.toArray(new Integer[numreads]);
        String star = "*";
        String error = "";

        //out.printf("INFO %s alleles in dictionary\n",HLAnames.length);
        String[][] A_exons = GetExonIntervals("A",true);
        String[][] B_exons = GetExonIntervals("B",false);
        String[][] C_exons = GetExonIntervals("C",false);
        String[][] DRB1_exons = GetExonIntervals("DRB1",false);
        String[][] DQB1_exons = GetExonIntervals("DQB1",false);
        String[][] DQA1_exons = GetExonIntervals("DQA1",true);
        String[][] DPB1_exons = GetExonIntervals("DPB1",true);
        String[][] DPA1_exons = GetExonIntervals("DPA1",false);
        //Print individual info and genotypes
        for (int i = 0; i < inputFileContents.length; i++){
            String[] s = inputFileContents[i].split(" ");
            //out.printf("%s\t%s\n",inputFileContents[i],s.length);
            if (s.length > 10){
                error = "";
                out.printf("%s\t%s\t%s\t%s\t%s\t%s",s[0],s[1],s[2],s[3],s[4],s[5]);
                String HLA_A_1 = GetAlleleName("HLA_A",star,s[6]);
                String HLA_A_2 = GetAlleleName("HLA_A",star,s[7]);
                String HLA_B_1 = GetAlleleName("HLA_B",star,s[8]);
                String HLA_B_2 = GetAlleleName("HLA_B",star,s[9]);
                String HLA_C_1 = GetAlleleName("HLA_C",star,s[10]);
                String HLA_C_2 = GetAlleleName("HLA_C",star,s[11]);
                String HLA_DPA1_1 = GetAlleleName("HLA_DPA1",star,s[12]);
                String HLA_DPA1_2 = GetAlleleName("HLA_DPA1",star,s[13]);
                String HLA_DPB1_1 = GetAlleleName("HLA_DPB1",star,s[14]);
                String HLA_DPB1_2 = GetAlleleName("HLA_DPB1",star,s[15]);
                String HLA_DQA1_1 = GetAlleleName("HLA_DQA1",star,s[16]);
                String HLA_DQA1_2 = GetAlleleName("HLA_DQA1",star,s[17]);
                String HLA_DQB1_1 = GetAlleleName("HLA_DQB1",star,s[18]);
                String HLA_DQB1_2 = GetAlleleName("HLA_DQB1",star,s[19]);
                String HLA_DRB1_1 = GetAlleleName("HLA_DRB1",star,s[20]);
                String HLA_DRB1_2 = GetAlleleName("HLA_DRB1",star,s[21]);

                

                if (true) {
                    error = error + PrintGenotypes(s[1], HLA_A_1,HLA_A_2, HLA_A_start,HLA_A_end);
                    error = error + PrintGenotypes(s[1], HLA_C_1,HLA_C_2, HLA_C_start,HLA_C_end);
                    error = error + PrintGenotypes(s[1], HLA_B_1,HLA_B_2, HLA_B_start,HLA_B_end);
                    error = error + PrintGenotypes(s[1], HLA_DRB1_1,HLA_DRB1_2, HLA_DRB1_start,HLA_DRB1_end);
                    error = error + PrintGenotypes(s[1], HLA_DQA1_1,HLA_DQA1_2, HLA_DQA1_start,HLA_DQA1_end);
                    error = error + PrintGenotypes(s[1], HLA_DQB1_1,HLA_DQB1_2, HLA_DQB1_start,HLA_DQB1_end);
                    error = error + PrintGenotypes(s[1], HLA_DPA1_1,HLA_DPA1_2, HLA_DPA1_start,HLA_DPA1_end);
                    error = error + PrintGenotypes(s[1], HLA_DPB1_1,HLA_DPB1_2, HLA_DPB1_start,HLA_DPB1_end);

                    error = error + PrintAminoAcids(s[1], HLA_A_1,HLA_A_2, A_exons);
                    error = error + PrintAminoAcids(s[1], HLA_C_1,HLA_C_2, C_exons);
                    error = error + PrintAminoAcids(s[1], HLA_B_1,HLA_B_2, B_exons);
                    error = error + PrintAminoAcids(s[1], HLA_DRB1_1,HLA_DRB1_2, DRB1_exons);
                    error = error + PrintAminoAcids(s[1], HLA_DQA1_1,HLA_DQA1_2, DQA1_exons);
                    error = error + PrintAminoAcids(s[1], HLA_DQB1_1,HLA_DQB1_2, DQB1_exons);
                    error = error + PrintAminoAcids(s[1], HLA_DPA1_1,HLA_DPA1_2, DPA1_exons);
                    error = error + PrintAminoAcids(s[1], HLA_DPB1_1,HLA_DPB1_2, DPB1_exons);
                    out.printf("\n");
                    out.printf("%s",error);
                }
            }
        }

        //Prints SNP names for each site
        if (true){
            PrintSNPS(HLA_A_start,HLA_A_end);
            PrintSNPS(HLA_C_start,HLA_C_end);
            PrintSNPS(HLA_B_start,HLA_B_end);
            PrintSNPS(HLA_DRB1_start,HLA_DRB1_end);
            PrintSNPS(HLA_DQA1_start,HLA_DQA1_end);
            PrintSNPS(HLA_DQB1_start,HLA_DQB1_end);
            PrintSNPS(HLA_DPA1_start,HLA_DPA1_end);
            PrintSNPS(HLA_DPB1_start,HLA_DPB1_end);

            PrintAminoAcidSites(A_exons,"A",true);
            PrintAminoAcidSites(C_exons,"C",false);
            PrintAminoAcidSites(B_exons,"B",false);
            PrintAminoAcidSites(DRB1_exons,"DRB1",false);
            PrintAminoAcidSites(DQA1_exons,"DQA1",true);
            PrintAminoAcidSites(DQB1_exons,"DQB1",false);
            PrintAminoAcidSites(DPA1_exons,"DPA1",false);
            PrintAminoAcidSites(DPB1_exons,"DPB1",true);
        }

    }

    private void PrintSNPS(int startpos, int stoppos){
        for (int pos = startpos; pos <= stoppos; pos++){
            SNPname = "CHR6_POS" + String.valueOf(pos);
            out.printf("6\t%s\t0\t%s\n",SNPname,pos);
        }
    }

    private void PrintAminoAcidSites(String[][] ExonIntervals, String locus, boolean isForwardStrand){
        int AAcount=1; int baseCount = 1; int exonNum;

        if (!isForwardStrand){
            for (int i = 1; i <= ExonIntervals.length; i++){
                int exonStart = Integer.parseInt(ExonIntervals[i-1][1]);
                int exonStop = Integer.parseInt(ExonIntervals[i-1][2]);
                for (int pos = exonStart; pos <= exonStop; pos++){
                    if (baseCount == 3){
                        AAcount++;
                        baseCount = 1;
                    }else{
                        baseCount++;
                    }
                }
            }
        }

        for (int i = 1; i <= ExonIntervals.length; i++){
            if (isForwardStrand){exonNum = i;}else{exonNum = ExonIntervals.length - i + 1;}
            int exonStart = Integer.parseInt(ExonIntervals[exonNum-1][1]);
            int exonStop = Integer.parseInt(ExonIntervals[exonNum-1][2]);
            for (int pos = exonStart; pos <= exonStop; pos++){
                if (baseCount == 2){
                    SNPname = locus + "_AA" + String.valueOf(AAcount) + "_E" + exonNum + "_" + String.valueOf(pos);
                    out.printf("6\t%s\t0\t%s\n",SNPname,pos);
                }
                if (baseCount == 3){
                    if (isForwardStrand){AAcount++;}else{AAcount--;}
                    baseCount = 1;
                }else{
                    baseCount++;
                }
            }
        }
    }
}
