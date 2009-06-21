package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.*;
import java.io.PrintStream;
import java.io.FileNotFoundException;

@WalkerName("CountCovariates")
public class CovariateCounterWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="buggyMaxReadLen", doc="If we see a read longer than this, we assume there's a bug and abort", required=false)
    public int buggyMaxReadLen = 100000;

    @Argument(fullName="OUTPUT_FILEROOT", shortName="outroot", required=false, doc="Filename root for the outputted logistic regression training files")
    public String OUTPUT_FILEROOT = "output";

    @Argument(fullName="MIN_MAPPING_QUALITY", shortName="minmap", required=false, doc="Only use reads with at least this quality score")
    public int MIN_MAPPING_QUALITY = 1;

    @Argument(fullName="READ_GROUP", shortName="rg", required=false, doc="Only use reads with this read group (@RG)")
    public String READ_GROUP = "none";

    //@Argument(fullName="MAX_READ_GROUPS", shortName="mrg", required=false, doc="Abort if number of read groups in input file exceeeds this count.")
    //public int MAX_READ_GROUPS = 100;

    @Argument(fullName="PLATFORM", shortName="pl", required=false, doc="Only calibrate read groups generated from the given platform (default = * for all platforms)")
    public List<String> platforms = Collections.singletonList("*");
    //public List<String> platforms = Collections.singletonList("ILLUMINA");

    @Argument(fullName="collapsePos", shortName="collapsePos", required=false, doc="")
    public boolean collapsePos = false;

    @Argument(fullName="collapseDinuc", shortName="collapseDinuc", required=false, doc="")
    public boolean collapseDinuc = false;

    HashMap<String, RecalDataManager> data = new HashMap<String, RecalDataManager>();

    long counted_sites = 0; // number of sites used to count covariates
    long counted_bases = 0; // number of bases used to count covariates
    long skipped_sites = 0; // number of sites skipped because of a dbSNP entry

    public void initialize() {
        //if( getToolkit().getEngine().getSAMHeader().getReadGroups().size() > MAX_READ_GROUPS )
        //    Utils.scareUser("Number of read groups in the specified file exceeds the number that can be processed in a reasonable amount of memory." +
        //                    "To override this limit, use the --MAX_READ_GROUPS (-mrg) parameter");

        for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
            if( readGroup.getAttribute("PL") == null )
                Utils.warnUser(String.format("PL attribute for read group %s is unset; assuming all reads are supported",readGroup.getReadGroupId()));
            if( !isSupportedReadGroup(readGroup) )
                continue;
            String rg = readGroup.getReadGroupId();
            //RecalDataManager manager = new RecalDataManager(rg, maxReadLen, QualityUtils.MAX_QUAL_SCORE+1, RecalData.NDINUCS, ! collapsePos, ! collapseDinuc );
            RecalDataManager manager = new RecalDataManager(rg, ! collapsePos, ! collapseDinuc );
            data.put(rg, manager);
        }
        out.printf("Created recalibration data collectors for %d read group(s)%n", data.size());
    }

    private RecalData getRecalData(String readGroup, int pos, int qual, char prevBase, char base) {
        byte[] cs = {(byte)prevBase, (byte)base};
        String s = new String(cs);
        return data.get(readGroup).expandingGetRecalData(pos, qual, s, true);
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

                if ( read.getReadLength() > buggyMaxReadLen ) {
                    throw new RuntimeException("Expectedly long read, please increase maxium read len with maxReadLen parameter: " + read.format());
                }

                SAMReadGroupRecord readGroup = read.getHeader().getReadGroup((String)read.getAttribute("RG"));
                if ( isSupportedReadGroup(readGroup) &&
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
        if ( qual > 0 ) { // && qual <= QualityUtils.MAX_QUAL_SCORE ) {
            // previous base is the next base in terms of machine chemistry if this is a negative strand
            // if ( qual == 2 )
            //if ( read.getReadName().equals("30PTAAAXX090126:5:14:132:764#0") )
            //    System.out.printf("Adding neg?=%b b_offset=%c offset=%d cycle=%d qual=%d dinuc=%c%c ref_match=%c comp=%c name=%s%n",
            //            read.getReadNegativeStrandFlag(), (char)read.getReadBases()[offset], offset, cycle, qual, prevBase, base,
            //            ref, (char)BaseUtils.simpleComplement(ref), read.getReadName());
            RecalData datum = getRecalData(rg, cycle, qual, prevBase, base);
            if (datum != null) datum.inc(base,ref);
            return 1;
        } else {
            return 0;
        }
    }

    public void onTraversalDone(Integer result) {
        printInfo(out);

        out.printf("Writing raw recalibration data%n");
        writeRecalTable();
        out.printf("...done%n");

        //out.printf("Writing logistic recalibration data%n");
        //writeLogisticRecalibrationTable();
        //out.printf("...done%n");
    }

    private void printInfo(PrintStream out) {
        out.printf("# date          %s%n", new Date());
        out.printf("# collapsed_pos %b%n", collapsePos);
        out.printf("# collapsed_dinuc %b%n", collapseDinuc);
        out.printf("# counted_sites %d%n", counted_sites);
        out.printf("# counted_bases %d%n", counted_bases);
        out.printf("# skipped_sites %d%n", skipped_sites);
        out.printf("# fraction_skipped 1/%.0f%n", (double)counted_sites / skipped_sites);
    }

    private void writeLogisticRecalibrationTable() {
        PrintStream dinuc_out = null;
        try {
            dinuc_out = new PrintStream( OUTPUT_FILEROOT+".covariate_counts.csv");
            dinuc_out.println("rg,dn,logitQ,pos,indicator,count");
            for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
                for ( int dinuc_index=0; dinuc_index<RecalData.NDINUCS; dinuc_index++) {
                    for ( RecalData datum: getRecalData(readGroup.getReadGroupId()) ) {
                        if ( RecalData.dinucIndex(datum.dinuc) == dinuc_index ) {
                            if ((datum.N - datum.B) > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup.getReadGroupId(), RecalData.dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 0, datum.N - datum.B);
                            if (datum.B > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup.getReadGroupId(), RecalData.dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 1, datum.B);
                        }
                    }
                }
            }
        }
        catch (FileNotFoundException e) {
            System.err.println("FileNotFoundException: " + e.getMessage());
        }
        finally {
            if (dinuc_out != null) dinuc_out.close();
        }
    }

     private void writeRecalTable() {
         PrintStream table_out = null;
         try {
             table_out = new PrintStream( OUTPUT_FILEROOT+".recal_data.csv");
             printInfo(table_out);
             table_out.println("rg,dn,Qrep,pos,NBases,MMismatches,Qemp");
             for (SAMReadGroupRecord readGroup : this.getToolkit().getEngine().getSAMHeader().getReadGroups()) {
                 for ( RecalData datum: getRecalData(readGroup.getReadGroupId()) ) {
                     if ( datum.N > 0 )
                         table_out.format("%s%n", datum.toCSVString(collapsePos));
                 }
             }
         } catch (FileNotFoundException e ) {
             System.err.println("FileNotFoundException: " + e.getMessage());
         }
         finally {
             if (table_out != null) table_out.close();
         }
     }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer a, Integer b) {
        return 0;
    }

    /**
     * Check to see whether this read group should be processed.
     * @param readGroup
     * @return
     */
    private boolean isSupportedReadGroup( SAMReadGroupRecord readGroup ) {
        for( String platform: platforms ) {
            platform = platform.trim();
            if( readGroup.getAttribute("PL") == null ||
                    platform.equals("*") || 
                    readGroup.getAttribute("PL").toString().equalsIgnoreCase(platform) )
                return true;
        }

        return false;
    }
}
