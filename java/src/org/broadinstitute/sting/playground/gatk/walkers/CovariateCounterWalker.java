package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.ArrayList;
import java.util.List;

@WalkerName("CountCovariates")
public class CovariateCounterWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="MAX_READ_LENGTH", shortName="mrl", required=false,defaultValue="101")
    public int MAX_READ_LENGTH;

    @Argument(fullName="MAX_QUAL_SCORE", shortName="mqs", required=false,defaultValue="34")
    public int MAX_QUAL_SCORE;
    //public double MAX_REPORTED_QUAL_SCORE = 100.0;

    int NDINUCS = 16;
    RecalData[][][] data = new RecalData[MAX_READ_LENGTH+1][MAX_QUAL_SCORE+1][NDINUCS];
    //RecalData[][][] data = new RecalData;
    ArrayList<RecalData> flattenData = new ArrayList();

    int nuc2num[];
    char num2nuc[];

    private class RecalData {
        long N;
        long B;
        int pos;
        int qual;
        String dinuc;

        public RecalData(int pos, int qual, String dinuc ) {
            this.pos = pos;
            this.qual = qual;
            this.dinuc = dinuc;
        }

        public void inc(char curBase, char ref) {
            N++;
            B += nuc2num[curBase] == nuc2num[ref] ? 0 : 1;
            //out.printf("%s %s\n", curBase, ref);
        }

        public String headerString() {
            return ("pos, dinuc, qual, emp_qual, qual_diff, n, b");
        }

        public String toString() {
            double empiricalQual = -10 * Math.log10((double)B / N);

            if (empiricalQual > MAX_QUAL_SCORE) empiricalQual = MAX_QUAL_SCORE;
            return String.format("%3d,%s,%3d,%5.1f,%5.1f,%6d,%6d", pos, dinuc, qual, empiricalQual, qual-empiricalQual, N, B);
            //return String.format("%d\t%s\t%d\t%.1f\t%d\t%6d", pos, dinuc, qual, empiricalQual, N, B);
        }
    }

    public void initialize() {
        for ( int i = 0; i < MAX_READ_LENGTH+1; i++) {
            for ( int j = 0; j < MAX_QUAL_SCORE+1; j++) {
                for ( int k = 0; k < NDINUCS; k++) {
                    String dinuc = dinucIndex2bases(k);
                    RecalData datum = new RecalData(i, j, dinuc);
                    data[i][j][k] = datum;
                    flattenData.add(datum);
                }
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp == null ) {
            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();
            for (int i =0; i < reads.size(); i++ ) {
                SAMRecord read = reads.get(i);
                int offset = offsets.get(i);
                if ( offset > 0 ) {
                    char base = (char)read.getReadBases()[offset];
                    char prevBase = (char)read.getReadBases()[offset-1];
                    int qual = (int)read.getBaseQualities()[offset];
                    //out.printf("%d %d %d\n", offset, qual, bases2dinucIndex(prevBase,base));
                    RecalData datum = data[offset][qual][bases2dinucIndex(prevBase,base)];
                    datum.inc(base,ref);
                }
            }
        }
        return 1;
    }

    public void onTraversalDone(Integer result) {
        if (flattenData.size() > 0)
            out.println(flattenData.get(0).headerString());
        for ( RecalData datum : flattenData ) {
            out.println(datum);
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer a, Integer b) {
        return 0;
    }

    int bases2dinucIndex(char prevBase, char base) {
        return nuc2num[prevBase] * 4 + nuc2num[base];
    };

    String dinucIndex2bases(int index) {
        char data[] = {num2nuc[index / 4], num2nuc[index % 4]};
        return new String( data );
    }

    public CovariateCounterWalker() {
        nuc2num = new int[128];
        nuc2num['A'] = 0;
        nuc2num['C'] = 1;
        nuc2num['G'] = 2;
        nuc2num['T'] = 3;
        nuc2num['a'] = 0;
        nuc2num['c'] = 1;
        nuc2num['g'] = 2;
        nuc2num['t'] = 3;

        num2nuc = new char[4];
        num2nuc[0] = 'A';
        num2nuc[1] = 'C';
        num2nuc[2] = 'G';
        num2nuc[3] = 'T';
    }
}