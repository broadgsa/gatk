/*
 * Creates a ped file given reads (for SNP analysis, imputation, etc)
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.ArrayList;
import java.util.Hashtable;
/**
 *
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CreatePedFileWalker extends ReadWalker<Integer, Integer> {
    @Argument(fullName = "allelesFile", shortName = "allelesFile", doc = "Create ped file for HLA alleles named in this file", required = true)
    public String alleleNamesFile = "";

    @Argument(fullName = "pedIntervals", shortName = "pedIntervals", doc = "Create genotypes in these intervals", required = false)
    public String pedIntervalsFile = "";

    String[] HLAnames, HLAreads, inputFileContents;
    Integer[] HLAstartpos, HLAstoppos;
    ArrayList<String> HLAnamesAL, HLAreadsAL;
    ArrayList<Integer> HLAstartposAL, HLAstopposAL;
    int[][] intervals;
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
        }
        return 0;
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


    public Integer map(char[] ref, SAMRecord read) {
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
            if (pos >= start1 && pos <= stop1){
                c1 = s1.charAt(pos-start1);
                if (c1 == 'D'){c1 = '0';}
            }else{
                c1 = '0';
            }

            if (pos >= start2 && pos <= stop2){
                c2 = s2.charAt(pos-start2);
                if (c2 == 'D'){c2 = '0';}
            }else{
                c2 = '0';
            }

            out.printf("\t%s %s",c1,c2);
        }
        return error;
    }

    private int GetAlleleIndex(String alleleName){
        int i;
        for (i = 0; i < HLAnames.length; i++){
            if (HLAnames[i].equals(alleleName)){
                return i;
            }
        }
        return -1;
        
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer numreads) {
        HLAnames = HLAnamesAL.toArray(new String[numreads]);
        HLAreads = HLAreadsAL.toArray(new String[numreads]);
        HLAstartpos = HLAstartposAL.toArray(new Integer[numreads]);
        HLAstoppos = HLAstopposAL.toArray(new Integer[numreads]);
        String star = "*";
        String error = "";

        //Print individual info and genotypes
        for (int i = 0; i < inputFileContents.length; i++){
            String[] s = inputFileContents[i].split(" ");
            //out.printf("%s\t%s\n",inputFileContents[i],s.length);
            if (s.length > 10){
                error = "";
                out.printf("%s\t%s\t%s\t%s\t%s\t%s",s[0],s[1],s[2],s[3],s[4],s[5]);
                String HLA_A_1 = "HLA_A" + star + s[6];
                String HLA_A_2 = "HLA_A" + star + s[7];
                String HLA_B_1 = "HLA_B" + star + s[8];
                String HLA_B_2 = "HLA_B" + star + s[9];
                String HLA_C_1 = "HLA_C" + star + s[10];
                String HLA_C_2 = "HLA_C" + star + s[11];
                String HLA_DPA1_1 = "HLA_DPA1" + star + s[12];
                String HLA_DPA1_2 = "HLA_DPA1" + star + s[13];
                String HLA_DPB1_1 = "HLA_DPB1" + star + s[14];
                String HLA_DPB1_2 = "HLA_DPB1" + star + s[15];
                String HLA_DQA1_1 = "HLA_DQA1" + star + s[16];
                String HLA_DQA1_2 = "HLA_DQA1" + star + s[17];
                String HLA_DQB1_1 = "HLA_DQB1" + star + s[18];
                String HLA_DQB1_2 = "HLA_DQB1" + star + s[19];
                String HLA_DRB1_1 = "HLA_DRB1" + star + s[20];
                String HLA_DRB1_2 = "HLA_DRB1" + star + s[21];

                error = error + PrintGenotypes(s[1], HLA_A_1,HLA_A_2, HLA_A_start,HLA_A_end);
                error = error + PrintGenotypes(s[1], HLA_C_1,HLA_C_2, HLA_C_start,HLA_C_end);
                error = error + PrintGenotypes(s[1], HLA_B_1,HLA_B_2, HLA_B_start,HLA_B_end);
                error = error + PrintGenotypes(s[1], HLA_DRB1_1,HLA_DRB1_2, HLA_DRB1_start,HLA_DRB1_end);
                error = error + PrintGenotypes(s[1], HLA_DQA1_1,HLA_DQA1_2, HLA_DQA1_start,HLA_DQA1_end);
                error = error + PrintGenotypes(s[1], HLA_DQB1_1,HLA_DQB1_2, HLA_DQB1_start,HLA_DQB1_end);
                error = error + PrintGenotypes(s[1], HLA_DPA1_1,HLA_DPA1_2, HLA_DPA1_start,HLA_DPA1_end);
                error = error + PrintGenotypes(s[1], HLA_DPB1_1,HLA_DPB1_2, HLA_DPB1_start,HLA_DPB1_end);
                
                out.printf("\n");
                out.printf("%s",error);
            }
        }

        //Prints SNP names for each site
        PrintSNPS(HLA_A_start,HLA_A_end);
        PrintSNPS(HLA_C_start,HLA_C_end);
        PrintSNPS(HLA_B_start,HLA_B_end);
        PrintSNPS(HLA_DRB1_start,HLA_DRB1_end);
        PrintSNPS(HLA_DQA1_start,HLA_DQA1_end);
        PrintSNPS(HLA_DQB1_start,HLA_DQB1_end);
        PrintSNPS(HLA_DPA1_start,HLA_DPA1_end);
        PrintSNPS(HLA_DPB1_start,HLA_DPB1_end);
    }

    private void PrintSNPS(int startpos, int stoppos){
        for (int pos = startpos; pos <= stoppos; pos++){
            SNPname = "CHR6_POS" + String.valueOf(pos);
            out.printf("6\t%s\t0\t%s\n",SNPname,pos);
        }
    }
}
