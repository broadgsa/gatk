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
import java.io.PrintStream;
import java.io.FileNotFoundException;

@WalkerName("CountCovariates")
public class CovariateCounterWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="MAX_READ_LENGTH", shortName="mrl", doc="max read length", required=false,defaultValue="101")
    public int MAX_READ_LENGTH;

    @Argument(fullName="MAX_QUAL_SCORE", shortName="mqs", doc="max quality score", required=false,defaultValue="63")
    public int MAX_QUAL_SCORE;

    @Argument(fullName="OUTPUT_FILEROOT", shortName="outroot", required=false, defaultValue="output", doc="Filename root for the outputted logistic regression training files")
    public String OUTPUT_FILEROOT;

    @Argument(fullName="CREATE_TRAINING_DATA", shortName="trainingdata", required=false, defaultValue="false", doc="Create training data files for logistic regression")
    public boolean CREATE_TRAINING_DATA;

    int NDINUCS = 16;
    RecalData[][][] data = new RecalData[MAX_READ_LENGTH+1][MAX_QUAL_SCORE+1][NDINUCS];
    //RecalData[][][] data = new RecalData;
    ArrayList<RecalData> flattenData = new ArrayList();

    static int nuc2num[];
    static char num2nuc[];

    String dinuc_root = "dinuc";
    ArrayList<PrintStream> dinuc_outs = new ArrayList<PrintStream>();
    PrintStream covars_out = new PrintStream("covars.out");
    PrintStream ByQualFile; // = new PrintStream("quality_empirical_vs_observed.csv");
    PrintStream ByCycleFile;
    PrintStream ByDinucFile;

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

        public void inc(long incN, long incB) {
            N += incN;
            B += incB;
        }


        public void inc(char curBase, char ref) {
            inc(1, nuc2num[curBase] == nuc2num[ref] ? 0 : 1);
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
        if ( dbsnp == null || !dbsnp.isSNP() ) {
            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();
            for (int i =0; i < reads.size(); i++ ) {
                
                SAMRecord read = reads.get(i);
                int offset = offsets.get(i);
                int numBases = read.getReadLength();
                if ( offset > 0 && offset < (numBases-1) ) { // skip first and last bases because they suck and they don't have a dinuc count
                    char base = (char)read.getReadBases()[offset];
                    int qual = (int)read.getBaseQualities()[offset];
                    //out.printf("%d %d %d\n", offset, qual, bases2dinucIndex(prevBase,base));
                    // previous base is the next base in terms of machine chemistry if this is a negative strand
                    char prevBase = (char)read.getReadBases()[offset + (read.getReadNegativeStrandFlag() ? 1 : -1)];
                    int dinuc_index = bases2dinucIndex(prevBase,base);
                    if (qual > 0 && qual <= MAX_QUAL_SCORE) {
                        // Convert offset into cycle position which means reversing the position of reads on the negative strand
                        int cycle = read.getReadNegativeStrandFlag() ? numBases - offset - 1 : offset;
                        data[cycle][qual][dinuc_index].inc(base,ref);
                        if (CREATE_TRAINING_DATA)
                            dinuc_outs.get(dinuc_index).format("%d,%d,%d\n", qual, cycle, nuc2num[base]==nuc2num[ref] ? 0 : 1);
                    }
                }
            }
        }
        return 1;
    }

    public void onTraversalDone(Integer result) {
        if (flattenData.size() > 0)
            covars_out.println(flattenData.get(0).headerString());
        for ( RecalData datum : flattenData ) {
            covars_out.println(datum);
        }

        qualityEmpiricalObserved();
        qualityDiffVsCycle();
        qualityDiffVsDinucleotide();

        // Close dinuc filehandles
        if (CREATE_TRAINING_DATA)
            for ( PrintStream dinuc_stream : this.dinuc_outs)
                dinuc_stream.close();

    }

    class MeanReportedQuality {
        double Qn = 0;
        long n = 0;

        void inc(double Q, long n) {
            this.n += n;
            Qn += Q * n;
        }

        double result() {
            return Qn / n;
        }
    }

    public void qualityDiffVsCycle() {
        ArrayList<RecalData> ByCycle = new ArrayList<RecalData>();
        ArrayList<MeanReportedQuality> ByCycleReportedQ = new ArrayList<MeanReportedQuality>();
        ByCycleFile.printf("cycle,Qemp-obs,Qemp,Qobs,B,N%n");
        for (int c=0; c < MAX_READ_LENGTH+1; c++)  {
            ByCycle.add(new RecalData(c, -1, "-"));
            ByCycleReportedQ.add(new MeanReportedQuality());
        }

        for ( RecalData datum: flattenData ) {
            ByCycle.get(datum.pos).inc(datum.N, datum.B);
            ByCycleReportedQ.get(datum.pos).inc(datum.qual, datum.N);
        }

        for (int c=0; c < MAX_READ_LENGTH+1; c++) {
            double empiricalQual = -10 * Math.log10((double)ByCycle.get(c).B / ByCycle.get(c).N);
            double reportedQual = ByCycleReportedQ.get(c).result();
            ByCycleFile.printf("%d, %f, %f, %f, %d, %d%n", c, empiricalQual-reportedQual, empiricalQual, reportedQual, ByCycle.get(c).B, ByCycle.get(c).N);
        }
    }

    public void qualityDiffVsDinucleotide() {
        ArrayList<RecalData> ByCycle = new ArrayList<RecalData>();
        ArrayList<MeanReportedQuality> ByCycleReportedQ = new ArrayList<MeanReportedQuality>();
        ByDinucFile.printf("dinuc,Qemp-obs,Qemp,Qobs,B,N%n");
        for (int c=0; c < NDINUCS; c++) {
            ByCycle.add(new RecalData(-1, -1, dinucIndex2bases(c)));
            ByCycleReportedQ.add(new MeanReportedQuality());
        }

        for ( RecalData datum: flattenData ) {
            int dinucIndex = bases2dinucIndex(datum.dinuc.charAt(0), datum.dinuc.charAt(1));
            ByCycle.get(dinucIndex).inc(datum.N, datum.B);
            ByCycleReportedQ.get(dinucIndex).inc(datum.qual, datum.N);
        }

        for (int c=0; c < NDINUCS; c++) {
            double empiricalQual = -10 * Math.log10((double)ByCycle.get(c).B / ByCycle.get(c).N);
            double reportedQual = ByCycleReportedQ.get(c).result();
            ByDinucFile.printf("%s, %f, %f, %f, %d, %d%n", ByCycle.get(c).dinuc, empiricalQual-reportedQual, empiricalQual, reportedQual, ByCycle.get(c).B, ByCycle.get(c).N);
        }
    }

    public void qualityEmpiricalObserved() {

        ArrayList<RecalData> ByQ = new ArrayList<RecalData>();
        ByQualFile.printf("Qrep,Qemp,B,N%n");
        for (int q=0; q<MAX_QUAL_SCORE+1; q++)
            ByQ.add(new RecalData(-1,q,"-"));

        for ( RecalData datum: flattenData ){
            ByQ.get(datum.qual).inc(datum.N, datum.B);
        }

        for (int q=0; q<MAX_QUAL_SCORE; q++) {
            double empiricalQual = -10 * Math.log10((double)ByQ.get(q).B / ByQ.get(q).N);
            ByQualFile.printf("%d, %f, %d, %d%n", q, empiricalQual, ByQ.get(q).B, ByQ.get(q).N);
            //out.printf("%3d,%s,%3d,%5.1f,%5.1f,%6d,%6d", pos, dinuc, qual, empiricalQual, qual-empiricalQual, N, B);
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer a, Integer b) {
        return 0;
    }

    public int bases2dinucIndex(char prevBase, char base) {
        return nuc2num[prevBase] * 4 + nuc2num[base];
    }

    public String dinucIndex2bases(int index) {
        char data[] = {num2nuc[index / 4], num2nuc[index % 4]};
        return new String( data );
    }

    static {
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

    // Print out data for regression
    public CovariateCounterWalker()  throws FileNotFoundException {
        if (CREATE_TRAINING_DATA)
            for ( int i=0; i<NDINUCS; i++) {
                PrintStream next_dinuc = new PrintStream( dinuc_root+"."+dinucIndex2bases(i)+".csv");
                next_dinuc.println("logitQ,pos,indicator");
                dinuc_outs.add(next_dinuc);
            }
        ByQualFile = new PrintStream(OUTPUT_FILEROOT+".empirical_v_reported_quality.csv");
        ByCycleFile = new PrintStream(OUTPUT_FILEROOT+".quality_difference_v_cycle.csv");
        ByDinucFile = new PrintStream(OUTPUT_FILEROOT+".quality_difference_v_dinucleotide.csv");
    }
}