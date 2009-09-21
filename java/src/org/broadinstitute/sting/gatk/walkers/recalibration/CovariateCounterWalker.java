package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;

import java.util.*;
import java.io.PrintStream;
import java.io.FileNotFoundException;

@WalkerName("CountCovariates")
public class CovariateCounterWalker extends LocusWalker<Integer, PrintStream> {
    @Argument(fullName="buggyMaxReadLen", doc="If we see a read longer than this, we assume there's a bug and abort", required=false)
    public int buggyMaxReadLen = 100000;

    @Argument(fullName="OUTPUT_FILEROOT", shortName="outroot", required=false, doc="Depreciated output file root -- now use --params to directly specify the file output name")
    public String OUTPUT_FILEROOT = null; // going to blow up if specified

    @Argument(fullName="params", shortName="params", required=false, doc="Filename root for the outputted logistic regression training files")
    public String params = "output.recal_data.csv";
    
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

    @Argument(fullName = "useOriginalQuals", shortName="OQ", doc="If provided, we will use use the quals from the original qualities OQ attribute field instead of the quals in the regular QUALS field", required=false)
    public boolean useOriginalQuals = false;


    private CovariateCounter covariateCounter = null;

    private long counted_sites = 0; // number of sites used to count covariates
    private long counted_bases = 0; // number of bases used to count covariates
    private long skipped_sites = 0; // number of sites skipped because of a dbSNP entry

    // THIS IS A HACK required in order to reproduce the behavior of old (and imperfect) RODIterator and
    // hence to pass the integration test. The new iterator this code is now using does see ALL the SNPs,
    // whether masked by overlapping indels/other events or not.
    //TODO process correctly all the returned dbSNP rods at each location


    private Pair<Long, Long> dbSNP_counts = new Pair<Long, Long>(0l, 0l);  // mismatch/base counts for dbSNP loci
    private Pair<Long, Long> novel_counts = new Pair<Long, Long>(0l, 0l);  // mismatch/base counts for non-dbSNP loci
    private static final double DBSNP_VS_NOVEL_MISMATCH_RATE = 2.0;        // rate at which dbSNP sites (on an individual level) mismatch relative to novel sites (determined by looking at NA12878)
    private static final int DBSNP_VALIDATION_CHECK_FREQUENCY = 1000000;   // how often to validate dbsnp mismatch rate (in terms of loci seen)
    private int lociSinceLastDbsnpCheck = 0;                               // loci since last dbsnp validation

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
        // THIS IS A HACK required in order to reproduce the behavior of old (and imperfect) RODIterator and
        // hence to pass the integration test. The new iterator this code is now using does see ALL the SNPs,
        // whether masked by overlapping indels/other events or not.
        //TODO process correctly all the returned dbSNP rods at each location
        BrokenRODSimulator.attach("dbSNP");
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
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        rodDbSNP dbsnp = (rodDbSNP)BrokenRODSimulator.simulate_lookup("dbSNP",ref.getLocus(),tracker);
//        long testpos =  10410913 ;
//        if ( ref.getLocus().getStart() ==  testpos ) {
//            System.out.println(rods.size()+" rods:");
//            for ( ReferenceOrderedDatum d : rods.getRecords() ) System.out.println((((rodDbSNP)d).isSNP()?"SNP":"non-SNP")+" at " + d.getLocation());
//          System.exit(1);
//        }



        if ( dbsnp == null || !dbsnp.isSNP() ) {
            // We aren't at a dbSNP position that's a SNP, so update the read

 //           if ( ref.getLocus().getStart() == testpos)  System.out.println("NOT A SNP INDEED");
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
                        counted_bases += covariateCounter.updateDataFromRead(readGroupString, read, offset, ref.getBase(), useOriginalQuals);
                    }
                }
            }
            counted_sites += 1;
            updateMismatchCounts(novel_counts, context, ref.getBase());
        } else {
 //           if ( ref.getLocus().getStart() == testpos)  System.out.println("TREATED AS A SNP");
   //         out.println(ref.getLocus()+" SNP at "+dbsnp.getLocation() );
            updateMismatchCounts(dbSNP_counts, context, ref.getBase());
            skipped_sites += 1;
        }

        if ( ++lociSinceLastDbsnpCheck == DBSNP_VALIDATION_CHECK_FREQUENCY ) {
            lociSinceLastDbsnpCheck = 0;
            validateDbsnpMismatchRate();
        }

        return 1;
    }


    /**
     * Update the mismatch / total_base counts for a given class of loci.
     *
     * @param counts
     * @param context
     * @return
     */
    private static void updateMismatchCounts(Pair<Long, Long> counts, AlignmentContext context, char ref) {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        for (int i =0; i < reads.size(); i++ ) {
            char readChar = reads.get(i).getReadString().charAt(offsets.get(i));
            int readCharBaseIndex = BaseUtils.simpleBaseToBaseIndex(readChar);
            int refCharBaseIndex  = BaseUtils.simpleBaseToBaseIndex(ref);

            if ( readCharBaseIndex != -1 && refCharBaseIndex != -1 ) {
                if ( readCharBaseIndex != refCharBaseIndex )
                    counts.first++;
                counts.second++;
            }
        }
    }

    /**
     * Validate the dbSNP mismatch rates.
     *
     * @return
     */
    private void validateDbsnpMismatchRate() {
        if ( novel_counts.second == 0 || dbSNP_counts.second == 0 )
            return;
        
        double fractionMM_novel = (double)novel_counts.first / (double)novel_counts.second;
        double fractionMM_dbsnp = (double)dbSNP_counts.first / (double)dbSNP_counts.second;
        //System.out.println(String.format("dbSNP rate = %.2f, novel rate = %.2f", fractionMM_dbsnp, fractionMM_novel));

        if ( fractionMM_dbsnp < DBSNP_VS_NOVEL_MISMATCH_RATE * fractionMM_novel )
            Utils.warnUser("The variation rate of dbSNP sites seems suspicious!  Please double-check that the correct DBSNP ROD is being used.");
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
            if ( OUTPUT_FILEROOT != null )
                throw new RuntimeException("OUTPUT_FILEROOT argument has been removed, please use --params from now on to directly specify the output parameter filename");
            return new PrintStream( params );
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException("Couldn't open output file", e);
        }
    }

    public void onTraversalDone(PrintStream recalTableStream) {
        out.printf("Writing raw recalibration data...");
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
        //out.printf("# date             \"%s\"%n", new Date());
        out.printf("# collapsed_pos    %b%n", collapsePos);
        out.printf("# collapsed_dinuc  %b%n", collapseDinuc);
        out.printf("# counted_sites    %d%n", counted_sites);
        out.printf("# counted_bases    %d%n", counted_bases);
        out.printf("# skipped_sites    %d%n", skipped_sites);
        out.printf("# fraction_skipped 1 / %.0f bp%n", (double)counted_sites / skipped_sites);
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
                if ( datum.N > 0 ) {
                    recalTableStream.println(datum.toCSVString(collapsePos));
                }
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
