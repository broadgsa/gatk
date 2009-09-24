/*
 * Creates a ped file given reads (for SNP analysis, imputation, etc)
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.playground.gatk.walkers.HLAcaller.ReadCigarFormatter;
import org.broadinstitute.sting.gatk.walkers.*;

import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.lang.Math;
/**
 *
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CreatePedFileWalker extends ReadWalker<Integer, Integer> {
    ReadCigarFormatter formatter = new ReadCigarFormatter();
    char c;
    boolean DEBUG = false;
    int HLA_A_start = 30018310;
    int HLA_A_end = 30021211;
    String[] SNPnames;
    String SNPname;
    int start, end;
    Integer I;

    Hashtable indexer = new Hashtable();

    public Integer reduceInit() {
        SNPnames = new String[HLA_A_end-HLA_A_start+1];
        start = HLA_A_start;
        end = HLA_A_end;

        indexer.put('A', (Integer) 1);
        indexer.put('C', (Integer) 2);
        indexer.put('G', (Integer) 3);
        indexer.put('T', (Integer) 4);
        indexer.put('D', (Integer) 0); // D for deletion
        out.print("Reads:\n");
        return 0;
    }



    public Integer map(char[] ref, SAMRecord read) {

        int readstart = read.getAlignmentStart();
        int readstop = read.getAlignmentEnd();

        if(readstart <= HLA_A_end && readstop >= HLA_A_start){
            String s = formatter.FormatRead(read.getCigarString(),read.getReadString());
            String name = read.getReadName();

            out.printf("%s\t%s\t0\t0\tM",name,name);

            for (int i = start; i <= end; i++){

                if (i - readstart < s.length() && i - readstart >= 0){
                    c = s.charAt(i-readstart);
                    I = (Integer) indexer.get(c);
                    out.printf("\t%s %s",I,I);
                }else{
                    out.print("\t0 0");
                }
            }
            out.printf("\n");
        }
        return 1;
    }



    public Integer reduce(Integer value, Integer sum) {

        return value + sum;
    }

    public void onTraversalDone(Integer value) {
        out.print("\nSNP names:\n");
        for (int pos = start; pos <= end; pos++){
            SNPname = "M CHR6_POS" + String.valueOf(pos);
            SNPnames[pos-start]=SNPname;
            out.printf("%s\n",SNPname);
        }
    }
}
