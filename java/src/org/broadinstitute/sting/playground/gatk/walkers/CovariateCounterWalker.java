package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.HashMap;
import java.io.PrintStream;
import java.io.FileNotFoundException;

@WalkerName("CountCovariates")
public class CovariateCounterWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="MAX_READ_LENGTH", shortName="mrl", doc="max read length", required=false)
    public int MAX_READ_LENGTH = 101;

    @Argument(fullName="MAX_QUAL_SCORE", shortName="mqs", doc="max quality score", required=false)
    public int MAX_QUAL_SCORE = 63;

    @Argument(fullName="OUTPUT_FILEROOT", shortName="outroot", required=false, doc="Filename root for the outputted logistic regression training files")
    public String OUTPUT_FILEROOT = "output";

    @Argument(fullName="CREATE_TRAINING_DATA", shortName="trainingdata", required=false, doc="Create training data files for logistic regression")
    public boolean CREATE_TRAINING_DATA = true;

    @Argument(fullName="DOWNSAMPLE_FRACTION", shortName="sample", required=false, doc="Fraction of bases to randomly sample")
    public float DOWNSAMPLE_FRACTION=1.0f;

    @Argument(fullName="MIN_MAPPING_QUALITY", shortName="minmap", required=false, doc="Only use reads with at least this quality score")
    public int MIN_MAPPING_QUALITY = 0;

    @Argument(fullName="READ_GROUP", shortName="rg", required=false, doc="Only use reads with this read group (@RG)")
    public String READ_GROUP = "none";

    int NDINUCS = 16;
    ArrayList<RecalData> flattenData = new ArrayList<RecalData>();
    HashMap<String, RecalData[][][]> data = new HashMap<String, RecalData[][][]>();
    //RecalData[][][] data;

    static int nuc2num[];
    static char num2nuc[];

    long counted_sites = 0; // number of sites used to count covariates
    long skipped_sites = 0; // number of sites skipped because of a dbSNP entry

    private class RecalData {
        long N;
        long B;
        int pos;
        int qual;
        String readGroup;
        String dinuc;

        public RecalData(int pos, int qual, String readGroup, String dinuc ) {
            this.pos = pos;
            this.qual = qual;
            this.readGroup = readGroup;
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
        for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
            data.put(readGroup.getReadGroupId(), new RecalData[MAX_READ_LENGTH+1][MAX_QUAL_SCORE+1][NDINUCS]);
            for ( int i = 0; i < MAX_READ_LENGTH+1; i++) {
                for ( int j = 0; j < MAX_QUAL_SCORE+1; j++) {
                    for ( int k = 0; k < NDINUCS; k++) {
                        String dinuc = dinucIndex2bases(k);
                        RecalData datum = new RecalData(i, j, readGroup.getReadGroupId(), dinuc);
                        data.get(readGroup.getReadGroupId())[i][j][k] = datum;
                        flattenData.add(datum);
                    }
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
                SAMReadGroupRecord readGroup = read.getHeader().getReadGroup((String)read.getAttribute("RG"));
                if ( "ILLUMINA".equalsIgnoreCase(readGroup.getAttribute("PL").toString()) &&
                    !read.getReadNegativeStrandFlag() &&
                    (READ_GROUP.equals("none") || read.getAttribute("RG") != null && read.getAttribute("RG").equals(READ_GROUP)) &&
                    (read.getMappingQuality() >= MIN_MAPPING_QUALITY)) {
                    //(random_genrator.nextFloat() <= DOWNSAMPLE_FRACTION)
                    int offset = offsets.get(i);
                    int numBases = read.getReadLength();
                    if ( offset > 0 && offset < (numBases-1) ) { // skip first and last bases because they suck and they don't have a dinuc count
                        int qual = (int)read.getBaseQualities()[offset];
                        if (qual > 0 && qual <= MAX_QUAL_SCORE) {
                            // previous base is the next base in terms of machine chemistry if this is a negative strand
                            char base = (char)read.getReadBases()[offset];
                            char prevBase = (char)read.getReadBases()[offset -1];
                            int dinuc_index = bases2dinucIndex(prevBase, base, false);
                            //char prevBase = (char)read.getReadBases()[offset + (read.getReadNegativeStrandFlag() ? 1 : -1)];
                            //int dinuc_index = bases2dinucIndex(prevBase, base, read.getReadNegativeStrandFlag());

                            // Convert offset into cycle position which means reversing the position of reads on the negative strand
                            //int cycle = read.getReadNegativeStrandFlag() ? numBases - offset - 1 : offset;
                            //data[cycle][qual][dinuc_index].inc(base,ref);
                            data.get(readGroup.getReadGroupId())[offset][qual][dinuc_index].inc(base,ref);
                        }
                    }
                }
            }
            counted_sites += 1;
        }else{
            skipped_sites += 1;
            //System.out.println(dbsnp.toSimpleString()+" "+new ReadBackedPileup(ref, context).getPileupString());
        }
        return 1;
    }

    public void onTraversalDone(Integer result) {
        PrintStream covars_out;
        try {
            covars_out = new PrintStream(OUTPUT_FILEROOT+".covars.out");
            if (flattenData.size() > 0)
                covars_out.println(flattenData.get(0).headerString());
            for ( RecalData datum : flattenData ) {
                covars_out.println(datum);
            }
        } catch (FileNotFoundException e) {
            System.err.println("FileNotFoundException: " + e.getMessage());
        }

        qualityEmpiricalObserved();
        qualityDiffVsCycle();
        qualityDiffVsDinucleotide();

        out.printf("Counted sites: %d%n", counted_sites);
        out.printf("Skipped sites: %d%n", skipped_sites);
        out.printf("Fraction skipped: 1/%.0f%n", (double)counted_sites / skipped_sites);

        if (CREATE_TRAINING_DATA) writeTrainingData();
    }

    void writeTrainingData() {
        PrintStream dinuc_out = null;
        try {
            dinuc_out = new PrintStream( OUTPUT_FILEROOT+".covariate_counts.csv");
            dinuc_out.println("rg,dn,logitQ,pos,indicator,count");            
            for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
                for ( int dinuc_index=0; dinuc_index<NDINUCS; dinuc_index++) {
                    for ( RecalData datum: flattenData ) {
                        if (datum.readGroup.equals(readGroup.getReadGroupId()) && string2dinucIndex(datum.dinuc) == dinuc_index) {
                            if ((datum.N - datum.B) > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup.getReadGroupId(), dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 0, datum.N - datum.B);
                            if (datum.B > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup.getReadGroupId(), dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 1, datum.B);
                        }
                    }
                }
            }
        }
        catch (FileNotFoundException e) {
            System.err.println("FileNotFoundException: " + e.getMessage());
            return;
        }
        finally {
            if (dinuc_out != null)
                dinuc_out.close();
        }

    }

    class MeanReportedQuality {
        double Qn = 0;
        long n = 0;
        double sumErrors = 0; // Count of estimated number of errors

        void inc(double Q, long n) {
            this.n += n;
            //Qn += Q * n; // wrong calculation that worked did not account for Qs being in log space
            sumErrors += QualityUtils.qualToErrorProb((byte)Q) * n;
        }

        double result() {
            //return Qn / n;
            return -10 * Math.log10(sumErrors / n);
            //return QualityUtils.probToQual(1.0 - (sumErrors / n));
        }
    }

    public void qualityDiffVsCycle() {
        for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
            PrintStream ByCycleFile = null;
            try {
                ByCycleFile = new PrintStream(OUTPUT_FILEROOT+".RG_"+readGroup.getReadGroupId()+".quality_difference_v_cycle.csv");
            } catch (FileNotFoundException e){
                System.out.println("Could not open output files based on OUTPUT_FILEROOT option: " + OUTPUT_FILEROOT);
                System.exit(1);
            }
            ArrayList<RecalData> ByCycle = new ArrayList<RecalData>();
            ArrayList<MeanReportedQuality> ByCycleReportedQ = new ArrayList<MeanReportedQuality>();
            ByCycleFile.printf("cycle,Qemp-obs,Qemp,Qobs,B,N%n");
            RecalData All = new RecalData(0,0,readGroup.getReadGroupId(),"");
            MeanReportedQuality AllReported = new MeanReportedQuality();
            for (int c=0; c < MAX_READ_LENGTH+1; c++)  {
                ByCycle.add(new RecalData(c, -1,readGroup.getReadGroupId(),"-"));
                ByCycleReportedQ.add(new MeanReportedQuality());
            }

            for ( RecalData datum: flattenData ) {
                ByCycle.get(datum.pos).inc(datum.N, datum.B);
                ByCycleReportedQ.get(datum.pos).inc(datum.qual, datum.N);
                All.inc(datum.N, datum.B);
                AllReported.inc(datum.qual, datum.N);
            }

            for (int c=0; c < MAX_READ_LENGTH+1; c++) {
                double empiricalQual = -10 * Math.log10((double)ByCycle.get(c).B / ByCycle.get(c).N);
                double reportedQual = ByCycleReportedQ.get(c).result();
                ByCycleFile.printf("%d, %f, %f, %f, %d, %d%n", c, empiricalQual-reportedQual, empiricalQual, reportedQual, ByCycle.get(c).B, ByCycle.get(c).N);
            }
        }
        //System.out.printf("Cycle: N=%d, B=%d, Qemp=%.1f, ", All.N, All.B, -10 * Math.log10((double)All.B/All.N));
        //System.out.printf("Qrep=%.1f%n", AllReported.result());
    }

    public void qualityDiffVsDinucleotide() {
        for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
            PrintStream ByDinucFile = null;
            try {
                ByDinucFile = new PrintStream(OUTPUT_FILEROOT+".RG_"+readGroup.getReadGroupId()+".quality_difference_v_dinucleotide.csv");
            } catch (FileNotFoundException e){
                System.out.println("Could not open output files based on OUTPUT_FILEROOT option: " + OUTPUT_FILEROOT);
                System.exit(1);
            }
            ArrayList<RecalData> ByCycle = new ArrayList<RecalData>();
            ArrayList<MeanReportedQuality> ByCycleReportedQ = new ArrayList<MeanReportedQuality>();
            ByDinucFile.printf("dinuc,Qemp-obs,Qemp,Qobs,B,N%n");
            RecalData All = new RecalData(0,0,readGroup.getReadGroupId(),"");
            MeanReportedQuality AllReported = new MeanReportedQuality();
            for (int c=0; c < NDINUCS; c++) {
                ByCycle.add(new RecalData(-1, -1,readGroup.getReadGroupId(),dinucIndex2bases(c)));
                ByCycleReportedQ.add(new MeanReportedQuality());
            }

            for ( RecalData datum: flattenData ) {
                int dinucIndex = string2dinucIndex(datum.dinuc); //bases2dinucIndex(datum.dinuc.charAt(0), datum.dinuc.charAt(1), false);
                ByCycle.get(dinucIndex).inc(datum.N, datum.B);
                ByCycleReportedQ.get(dinucIndex).inc(datum.qual, datum.N);
                All.inc(datum.N, datum.B);
                AllReported.inc(datum.qual, datum.N);
            }

            for (int c=0; c < NDINUCS; c++) {
                double empiricalQual = -10 * Math.log10((double)ByCycle.get(c).B / ByCycle.get(c).N);
                double reportedQual = ByCycleReportedQ.get(c).result();
                ByDinucFile.printf("%s, %f, %f, %f, %d, %d%n", ByCycle.get(c).dinuc, empiricalQual-reportedQual, empiricalQual, reportedQual, ByCycle.get(c).B, ByCycle.get(c).N);
            }
        }
        //System.out.printf("Dinuc: N=%d, B=%d, Qemp=%.1f, ", All.N, All.B, -10 * Math.log10((double)All.B/All.N));
        //System.out.printf("Qrep=%.1f%n", AllReported.result());
    }

    public void qualityEmpiricalObserved() {
        for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
            PrintStream ByQualFile  = null;
            try {
                ByQualFile = new PrintStream(OUTPUT_FILEROOT+".RG_"+readGroup.getReadGroupId()+".empirical_v_reported_quality.csv");
            } catch (FileNotFoundException e){
                System.out.println("Could not open output files based on OUTPUT_FILEROOT option: " + OUTPUT_FILEROOT);
                System.exit(1);
            }
            ArrayList<RecalData> ByQ = new ArrayList<RecalData>();
            ArrayList<MeanReportedQuality> ByQReportedQ = new ArrayList<MeanReportedQuality>();
            ByQualFile.printf("Qrep,Qemp,Qrep_avg,B,N%n");
            RecalData All = new RecalData(0,0,readGroup.getReadGroupId(),"");
            MeanReportedQuality AllReported = new MeanReportedQuality();
            for (int q=0; q<MAX_QUAL_SCORE+1; q++) {
                ByQ.add(new RecalData(-1,q,readGroup.getReadGroupId(),"-"));
                ByQReportedQ.add(new MeanReportedQuality());
            }

            for ( RecalData datum: flattenData ){
                ByQ.get(datum.qual).inc(datum.N, datum.B);
                ByQReportedQ.get(datum.qual).inc(datum.qual, datum.N);
                All.inc(datum.N, datum.B);
                AllReported.inc(datum.qual, datum.N);
                //out.printf("%2d%6d%3d %2d %s%n", datum.qual, datum.N, datum.pos, datum.qual, datum.dinuc);
            }

            for (int q=0; q<MAX_QUAL_SCORE; q++) {
                double empiricalQual = -10 * Math.log10((double)ByQ.get(q).B / ByQ.get(q).N);
                ByQualFile.printf("%d, %f, %.0f, %d, %d%n", q, empiricalQual, ByQReportedQ.get(q).result(), ByQ.get(q).B, ByQ.get(q).N);
                //out.printf("%3d,%s,%3d,%5.1f,%5.1f,%6d,%6d", pos, dinuc, qual, empiricalQual, qual-empiricalQual, N, B);                                                                                      n
            }
        }
        //System.out.printf("Emp-Obs: N=%d, B=%d, Qemp=%.1f, ", All.N, All.B, -10 * Math.log10((double)All.B/All.N));
        //System.out.printf("Qrep=%.1f%n", AllReported.result());
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer a, Integer b) {
        return 0;
    }

    public int bases2dinucIndex(char prevBase, char base, boolean Complement) {
        if (!Complement) {
            return nuc2num[prevBase] * 4 + nuc2num[base];
        }else{
            return (3 - nuc2num[prevBase]) * 4 + (3 - nuc2num[base]);
        }
    }

    public String dinucIndex2bases(int index) {
        char data[] = {num2nuc[index / 4], num2nuc[index % 4]};
        return new String( data );
    }

    public int string2dinucIndex(String s) {
        return bases2dinucIndex(s.charAt(0), s.charAt(1), false);
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
    Random random_genrator;
    // Print out data for regression
    public CovariateCounterWalker()  throws FileNotFoundException {
        random_genrator = new Random(123454321); // keep same random seed while debugging
    }
}
