package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.playground.gatk.walkers.RecalData;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.HashMap;
import java.util.Collections;
import java.io.PrintStream;
import java.io.FileNotFoundException;

@WalkerName("CountCovariates")
public class CovariateCounterWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="maxReadLen", shortName="mrl", doc="max read length", required=false)
    public int maxReadLen = 101;

    //@Argument(fullName="MAX_QUAL_SCORE", shortName="mqs", doc="max quality score", required=false)
    //public int MAX_QUAL_SCORE = 63;

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

    @Argument(fullName="MAX_READ_GROUPS", shortName="mrg", required=false, doc="Abort if number of read groups in input file exceeeds this count.")
    public int MAX_READ_GROUPS = 100;

    @Argument(fullName="PLATFORM", shortName="pl", required=false, doc="Only calibrate read groups generated from the given platform (default = Illumina)")
    public List<String> platforms = Collections.singletonList("ILLUMINA");

    @Argument(fullName="rawData", shortName="raw", required=false, doc="If true, raw mismatch observations will be output to a file")
    public boolean outputRawData = true;

    int NDINUCS = 16;
    //ArrayList<RecalData> flattenData = new ArrayList<RecalData>();
    //HashMap<String, RecalData[][][]> data = new HashMap<String, RecalData[][][]>();
    HashMap<String, RecalDataManager> data = new HashMap<String, RecalDataManager>();
    //RecalData[][][] data;
    boolean trackPos = true;
    boolean trackDinuc = true;

    long counted_sites = 0; // number of sites used to count covariates
    long counted_bases = 0; // number of bases used to count covariates
    long skipped_sites = 0; // number of sites skipped because of a dbSNP entry

    public void initialize() {
        if( getToolkit().getEngine().getSAMHeader().getReadGroups().size() > MAX_READ_GROUPS )
            Utils.scareUser("Number of read groups in the specified file exceeds the number that can be processed in a reasonable amount of memory." +
                            "To override this limit, use the --MAX_READ_GROUPS (-mrg) parameter");

        for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
            if( readGroup.getAttribute("PL") == null )
                Utils.warnUser(String.format("PL attribute for read group %s is unset; assuming all reads are supported",readGroup.getReadGroupId()));
            if( !isSupportedReadGroup(readGroup) )
                continue;
            String rg = readGroup.getReadGroupId();
            RecalDataManager manager = new RecalDataManager(rg, maxReadLen, QualityUtils.MAX_QUAL_SCORE+1, NDINUCS, trackPos, trackDinuc );
            //data.put(rg, new RecalData[maxReadLen+1][QualityUtils.MAX_QUAL_SCORE+1][NDINUCS]);
            data.put(rg, manager);
        }
    }

    private RecalData getRecalData(String readGroup, int pos, int qual, int dinuc_index) {
        return data.get(readGroup).expandingGetRecalData(pos, qual, dinuc_index, true);
    }

    private List<RecalData> getRecalData(String readGroup) {
        return data.get(readGroup).getAll();
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        //System.out.printf("%s %c%n", context.getLocation(), ref);
        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp == null || !dbsnp.isSNP() ) {
            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();
            for (int i =0; i < reads.size(); i++ ) {
                SAMRecord read = reads.get(i);

                if ( read.getReadLength() > maxReadLen ) {
                    throw new RuntimeException("Expectedly long read, please increase maxium read len with maxReadLen parameter: " + read.format());
                }

                SAMReadGroupRecord readGroup = read.getHeader().getReadGroup((String)read.getAttribute("RG"));
                if ( isSupportedReadGroup(readGroup) &&
                    //!read.getReadNegativeStrandFlag() &&
                    (READ_GROUP.equals("none") || read.getAttribute("RG") != null && read.getAttribute("RG").equals(READ_GROUP)) &&
                    (read.getMappingQuality() >= MIN_MAPPING_QUALITY)) {
                    int offset = offsets.get(i);
                    int numBases = read.getReadLength();
                    if ( offset > 0 && offset < (numBases-1) ) { // skip first and last bases because they suck and they don't have a dinuc count
                        counted_bases += updateDataFromRead(readGroup.getReadGroupId(), read, offset, ref);
                    }
                }
            }
            counted_sites += 1;
        } else {
            skipped_sites += 1;
            //System.out.println(dbsnp.toSimpleString()+" "+new ReadBackedPileup(ref, context).getPileupString());
        }
        return 1;
    }

    private int updateDataFromRead( String rg, SAMRecord read, int offset, char ref ) {
        int cycle = offset;
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();

        char base = (char)bases[offset];
        char prevBase = (char)bases[offset - 1];        

        if (read.getReadNegativeStrandFlag()) {
            ref = (char)BaseUtils.simpleComplement(ref);
            base = (char)BaseUtils.simpleComplement(base);
            prevBase = (char)BaseUtils.simpleComplement((char)bases[offset+1]);
            cycle = read.getReadLength() - (offset + 1);
        }

        int qual = quals[offset];
        if ( qual > 0 && qual <= QualityUtils.MAX_QUAL_SCORE ) {
            // previous base is the next base in terms of machine chemistry if this is a negative strand
            int dinuc_index = RecalData.bases2dinucIndex(prevBase, base, false);
            //System.out.printf("Adding b_offset=%c offset=%d cycle=%d qual=%d dinuc=%c%c ref_match=%c comp=%c%n", (char)read.getReadBases()[offset], offset, cycle, qual, prevBase, base, ref, (char)BaseUtils.simpleComplement(ref));
            getRecalData(rg, cycle, qual, dinuc_index).inc(base,ref);
            return 1;
        } else {
            return 0;
        }
    }

    public void onTraversalDone(Integer result) {
        PrintStream covars_out;
        try {
            covars_out = new PrintStream(OUTPUT_FILEROOT+".covars.out");
            covars_out.println(RecalData.headerString());
            for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
                for ( RecalData datum : getRecalData(readGroup.getReadGroupId()) ) {
                    covars_out.println(datum);
                }
            }
        } catch (FileNotFoundException e) {
            System.err.println("FileNotFoundException: " + e.getMessage());
        }

        qualityEmpiricalObserved();
        qualityDiffVsCycle();
        qualityDiffVsDinucleotide();

        out.printf("Counted sites: %d%n", counted_sites);
        out.printf("Counted bases: %d%n", counted_bases);
        out.printf("Skipped sites: %d%n", skipped_sites);
        out.printf("Fraction skipped: 1/%.0f%n", (double)counted_sites / skipped_sites);

        if (CREATE_TRAINING_DATA) writeTrainingData();
    }

    void writeTrainingData() {
        PrintStream dinuc_out = null;
        PrintStream table_out = null;
        try {
            dinuc_out = new PrintStream( OUTPUT_FILEROOT+".covariate_counts.csv");
            dinuc_out.println("rg,dn,logitQ,pos,indicator,count");
            for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
                for ( int dinuc_index=0; dinuc_index<NDINUCS; dinuc_index++) {
                    for ( RecalData datum: getRecalData(readGroup.getReadGroupId()) ) {
                        if ( RecalData.string2dinucIndex(datum.dinuc) == dinuc_index ) {
                            if ((datum.N - datum.B) > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup.getReadGroupId(), RecalData.dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 0, datum.N - datum.B);
                            if (datum.B > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup.getReadGroupId(), RecalData.dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 1, datum.B);
                        }
                    }
                }
            }

            if ( outputRawData ) {
                table_out = new PrintStream( OUTPUT_FILEROOT+".raw_data.csv");
                for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
                    for ( RecalData datum: getRecalData(readGroup.getReadGroupId()) ) {
                        if ( datum.N > 0 )
                            table_out.format("%s%n", datum.toCSVString());
                    }
                }
            }
        }
        catch (FileNotFoundException e) {
            System.err.println("FileNotFoundException: " + e.getMessage());
            return;
        }
        finally {
            if (dinuc_out != null) dinuc_out.close();
            if (table_out != null) table_out.close();
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
            for (int c=0; c < maxReadLen+1; c++)  {
                ByCycle.add(new RecalData(c, -1,readGroup.getReadGroupId(),"-"));
                ByCycleReportedQ.add(new MeanReportedQuality());
            }

            for ( RecalData datum: getRecalData(readGroup.getReadGroupId()) ) {
                ByCycle.get(datum.pos).inc(datum.N, datum.B);
                ByCycleReportedQ.get(datum.pos).inc(datum.qual, datum.N);
                All.inc(datum.N, datum.B);
                AllReported.inc(datum.qual, datum.N);
            }

            for (int c=0; c < maxReadLen+1; c++) {
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
                ByCycle.add(new RecalData(-1, -1,readGroup.getReadGroupId(),RecalData.dinucIndex2bases(c)));
                ByCycleReportedQ.add(new MeanReportedQuality());
            }

            for ( RecalData datum: getRecalData(readGroup.getReadGroupId()) ) {
                int dinucIndex = RecalData.string2dinucIndex(datum.dinuc); //bases2dinucIndex(datum.dinuc.charAt(0), datum.dinuc.charAt(1), false);
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
            for (int q=0; q<QualityUtils.MAX_QUAL_SCORE+1; q++) {
                ByQ.add(new RecalData(-1,q,readGroup.getReadGroupId(),"-"));
                ByQReportedQ.add(new MeanReportedQuality());
            }

            for ( RecalData datum: getRecalData(readGroup.getReadGroupId()) ){
                ByQ.get(datum.qual).inc(datum.N, datum.B);
                ByQReportedQ.get(datum.qual).inc(datum.qual, datum.N);
                All.inc(datum.N, datum.B);
                AllReported.inc(datum.qual, datum.N);
                //out.printf("%2d%6d%3d %2d %s%n", datum.qual, datum.N, datum.pos, datum.qual, datum.dinuc);
            }

            for (int q=0; q<QualityUtils.MAX_QUAL_SCORE; q++) {
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

    Random random_genrator;
    // Print out data for regression
    public CovariateCounterWalker()  throws FileNotFoundException {
        random_genrator = new Random(123454321); // keep same random seed while debugging
    }

    /**
     * Check to see whether this read group should be processed.
     * @param readGroup
     * @return
     */
    private boolean isSupportedReadGroup( SAMReadGroupRecord readGroup ) {
        for( String platform: platforms ) {
            platform = platform.trim();
            if( readGroup.getAttribute("PL") == null || readGroup.getAttribute("PL").toString().equalsIgnoreCase(platform) )
                return true;
        }

        return false;
    }
}
