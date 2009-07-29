package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;

import java.util.*;
import java.io.PrintStream;
import java.io.FileNotFoundException;

@WalkerName("CountCovariates")
public class CovariateCounterWalker extends LocusWalker<Integer, PrintStream> {
    @Argument(fullName="buggyMaxReadLen", doc="If we see a read longer than this, we assume there's a bug and abort", required=false)
    public int buggyMaxReadLen = 100000;

    @Argument(fullName="OUTPUT_FILEROOT", shortName="outroot", required=false, doc="Filename root for the outputted logistic regression training files")
    public String OUTPUT_FILEROOT = "output";

    @Argument(fullName="MIN_MAPPING_QUALITY", shortName="minmap", required=false, doc="Only use reads with at least this quality score")
    public int MIN_MAPPING_QUALITY = 1;

    @Argument(fullName="PLATFORM", shortName="pl", required=false, doc="Only calibrate read groups generated from the given platform (default = * for all platforms)")
    public List<String> platforms = Collections.singletonList("*");
    //public List<String> platforms = Collections.singletonList("ILLUMINA");

    @Argument(fullName="assumeFaultyHeader", required=false, doc="")
    public boolean assumeFaultyHeader = false;

    //@Argument(fullName="collapsePos", shortName="collapsePos", required=false, doc="")
    public boolean collapsePos = false;

    //@Argument(fullName="collapseDinuc", shortName="collapseDinuc", required=false, doc="")
    public boolean collapseDinuc = false;

    private CovariateCounter covariateCounter = null;

    private long counted_sites = 0; // number of sites used to count covariates
    private long counted_bases = 0; // number of bases used to count covariates
    private long skipped_sites = 0; // number of sites skipped because of a dbSNP entry


    /**
     * Initialize the system.  Setup the data CovariateCountry for the read groups in our header
     */
    public void initialize() {
        Set<String> readGroups = new HashSet<String>();
        for (SAMReadGroupRecord readGroup : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            if( readGroup.getAttribute("PL") == null )
                Utils.warnUser(String.format("PL attribute for read group %s is unset; assuming all reads are supported",readGroup.getReadGroupId()));
            if( !isSupportedReadGroup(readGroup) )
                continue;
            readGroups.add(readGroup.getReadGroupId());
        }

        covariateCounter = new CovariateCounter(readGroups, collapsePos, collapseDinuc, assumeFaultyHeader);
        logger.info(String.format("Created recalibration data collectors for %d read group(s)", covariateCounter.getNReadGroups()));
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Walk over each read in the locus pileup and update the covariate counts based on these bases and their
     * matching (or not) with ref.  dbSNP aware, so avoids sites that are known as SNPs in DBSNP.
     *
     * @param tracker
     * @param ref
     * @param context
     * @return
     */
    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp == null || !dbsnp.isSNP() ) {
            // We aren't at a dbSNP position that's a SNP, so update the read

            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();
            for (int i =0; i < reads.size(); i++ ) {
                SAMRecord read = reads.get(i);

                if ( read.getReadLength() > buggyMaxReadLen ) {
                    throw new RuntimeException("Expectedly long read, please increase maxium read len with maxReadLen parameter: " + read.format());
                }

                final String readGroupString = ((String)read.getAttribute("RG"));
                SAMReadGroupRecord readGroup = read.getHeader().getReadGroup(readGroupString);

                if ( readGroupString == null ) {
                    throw new RuntimeException("No read group annotation found for read " + read.format());
                }

                if ((read.getMappingQuality() >= MIN_MAPPING_QUALITY && isSupportedReadGroup(readGroup) )) {
                    int offset = offsets.get(i);
                    if ( offset > 0 && offset < (read.getReadLength() - 1) ) { // skip first and last bases because they suck and they don't have a dinuc count
                        counted_bases += covariateCounter.updateDataFromRead(readGroupString, read, offset, ref);
                    }
                }
            }
            counted_sites += 1;
        } else {
            skipped_sites += 1;
        }
        return 1;
    }


    /**
     * Check to see whether this read group should be processed.  Returns true if the
     * read group is in the list of platforms to process or the platform == *, indicating
     * that all platforms should be processed.
     *
     * @param readGroup
     * @return
     */
    private boolean isSupportedReadGroup( SAMReadGroupRecord readGroup ) {
        for( String platform: platforms ) {
            platform = platform.trim();
            if( platform.equals("*") || readGroup == null ||
                    readGroup.getAttribute("PL") == null ||
                    readGroup.getAttribute("PL").toString().equalsIgnoreCase(platform) )
                return true;
        }

        return false;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Reduce
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Provide an initial value for reduce computations.
     * @return Initial value of reduce.
     */
    public PrintStream reduceInit() {
        try {
            return new PrintStream( OUTPUT_FILEROOT+".recal_data.csv" );
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException("Couldn't open output file", e);
        }
    }

    public void onTraversalDone(PrintStream recalTableStream) {
        printInfo(out);

        out.printf("Writing raw recalibration data..."); out.flush();
        writeRecalTable(recalTableStream);
        out.printf("...done%n");

        //out.printf("Writing logistic recalibration data%n");
        //writeLogisticRecalibrationTable();
        //out.printf("...done%n");

        recalTableStream.close();
    }

    /**
     * Prints some basic information about the CountCovariates run to the output stream out
     * @param out
     */
    private void printInfo(PrintStream out) {
        out.printf("# date             \"%s\"%n", new Date());
        out.printf("# collapsed_pos    %b%n", collapsePos);
        out.printf("# collapsed_dinuc  %b%n", collapseDinuc);
        out.printf("# counted_sites    %d%n", counted_sites);
        out.printf("# counted_bases    %d%n", counted_bases);
        out.printf("# skipped_sites    %d%n", skipped_sites);
        out.printf("# fraction_skipped 1 / %.0f bp%n", (double)counted_sites / skipped_sites);
    }

    @Deprecated
    private void writeLogisticRecalibrationTable() {
        PrintStream dinuc_out = null;
        try {
            dinuc_out = new PrintStream( OUTPUT_FILEROOT+".covariate_counts.csv");
            dinuc_out.println("rg,dn,logitQ,pos,indicator,count");
            for (String readGroup : covariateCounter.getReadGroups()) {
                for ( int dinuc_index=0; dinuc_index<RecalData.NDINUCS; dinuc_index++) {
                    for ( RecalData datum: covariateCounter.getRecalData(readGroup) ) {
                        if ( RecalData.dinucIndex(datum.dinuc) == dinuc_index ) {
                            if ((datum.N - datum.B) > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup, RecalData.dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 0, datum.N - datum.B);
                            if (datum.B > 0)
                                dinuc_out.format("%s,%s,%d,%d,%d,%d%n", readGroup, RecalData.dinucIndex2bases(dinuc_index), datum.qual, datum.pos, 1, datum.B);
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

    /**
     * Writes out the key recalibration data collected from the reads.  Dumps this recalibration data
     * as a CVS string to the recalTableOut PrintStream.  Emits the data for all read groups into this file.
     */
    private void writeRecalTable(PrintStream recalTableStream) {
        printInfo(recalTableStream);

        recalTableStream.println("rg,pos,Qrep,dn,nBases,nMismatches,Qemp");
        for (String readGroup : new TreeSet<String>(covariateCounter.getReadGroups()) ) {
            for ( RecalData datum: RecalData.sort(covariateCounter.getRecalData(readGroup)) ) {
                if ( datum.N > 0 )
                    recalTableStream.println(datum.toCSVString(collapsePos));
            }
        }
    }

    /**
     * Doesn't do anything
     *
     * @param empty
     * @param recalTableStream
     * @return
     */
    public PrintStream reduce(Integer empty, PrintStream recalTableStream) {
        return recalTableStream;
    }
}
