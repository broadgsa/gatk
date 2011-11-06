/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;
import net.sf.samtools.util.RuntimeIOException;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.NWaySAMFileWriter;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Performs local realignment of reads based on misalignments due to the presence of indels.
 *
 * <p>
 * The local realignment tool is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases
 * is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion
 * or deletion (indels) in the individual's genome with respect to the reference genome.  Such alignment artifacts result in many bases mismatching
 * the reference near the misalignment, which are easily mistaken as SNPs.  Moreover, since read mapping algorithms operate on each read independently,
 * it is impossible to place reads on the reference genome such at mismatches are minimized across all reads.  Consequently, even when some reads are
 * correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel,
 * also requiring realignment.  Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus
 * indel suitable for standard variant discovery approaches.  Unlike most mappers, this walker uses the full alignment context to determine whether an
 * appropriate alternate reference (i.e. indel) exists.  Following local realignment, the GATK tool Unified Genotyper can be used to sensitively and
 * specifically identify indels.
 * <p>
 *     <ol>There are 2 steps to the realignment process:
 *     <li>Determining (small) suspicious intervals which are likely in need of realignment (see the RealignerTargetCreator tool)</li>
 *     <li>Running the realigner over those intervals (IndelRealigner)</li>
 *     </ol>
 *     <p>
 * An important note: the input bam(s), reference, and known indel file(s) should be the same ones used for the RealignerTargetCreator step.
 * <p>
 * Another important note: because reads produced from the 454 technology inherently contain false indels, the realigner will not currently work with them
 * (or with reads from similar technologies).
 *
 * <h2>Input</h2>
 * <p>
 * One or more aligned BAM files and optionally one or more lists of known indels.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A realigned version of your input BAM file(s).
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -I input.bam \
 *   -R ref.fasta \
 *   -T IndelRealigner \
 *   -targetIntervals intervalListFromRTC.intervals \
 *   -o realignedBam.bam \
 *   [--known /path/to/indels.vcf] \
 *   [-compress 0]    (this argument recommended to speed up the process *if* this is only a temporary file; otherwise, use the default value)
 * </pre>
 *
 * @author ebanks
 */
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.ON_OUTPUT)
public class IndelRealigner extends ReadWalker<Integer, Integer> {

    public static final String ORIGINAL_CIGAR_TAG = "OC";
    public static final String ORIGINAL_POSITION_TAG = "OP";
    public static final String PROGRAM_RECORD_NAME = "GATK IndelRealigner";

    public enum ConsensusDeterminationModel {
        /**
         * Uses only indels from a provided ROD of known indels.
         */
        KNOWNS_ONLY,
        /**
         * Additionally uses indels already present in the original alignments of the reads.
         */
        USE_READS,
        /**
         * Additionally uses 'Smith-Waterman' to generate alternate consenses.
         */
        USE_SW
    }

    /**
     * Any number of VCF files representing known indels to be used for constructing alternate consenses.
     * Could be e.g. dbSNP and/or official 1000 Genomes indel calls.  Non-indel variants in these files will be ignored.
     */
    @Input(fullName="knownAlleles", shortName = "known", doc="Input VCF file(s) with known indels", required=false)
    public List<RodBinding<VariantContext>> known = Collections.emptyList();

    /**
     * The interval list output from the RealignerTargetCreator tool using the same bam(s), reference, and known indel file(s).
     */
    @Input(fullName="targetIntervals", shortName="targetIntervals", doc="intervals file output from RealignerTargetCreator", required=true)
    protected IntervalBinding<Feature> intervalsFile = null;

    /**
     * This term is equivalent to "significance" - i.e. is the improvement significant enough to merit realignment? Note that this number
     * should be adjusted based on your particular data set. For low coverage and/or when looking for indels with low allele frequency,
     * this number should be smaller.
     */
    @Argument(fullName="LODThresholdForCleaning", shortName="LOD", doc="LOD threshold above which the cleaner will clean", required=false)
    protected double LOD_THRESHOLD = 5.0;

    /**
     * The realigned bam file.
     */
    @Output(required=false, doc="Output bam")
    protected StingSAMFileWriter writer = null;
    protected ConstrainedMateFixingManager manager = null;
    protected SAMFileWriter writerToUse = null;

    /**
     * We recommend that users run with USE_READS when trying to realign high quality longer read data mapped with a gapped aligner;
     * Smith-Waterman is really only necessary when using an ungapped aligner (e.g. MAQ in the case of single-end read data).
     */
    @Argument(fullName = "consensusDeterminationModel", shortName = "model", doc = "Determines how to compute the possible alternate consenses", required = false)
    public ConsensusDeterminationModel consensusModel = ConsensusDeterminationModel.USE_READS;


    // ADVANCED OPTIONS FOLLOW

    /**
     * For expert users only!  This is similar to the argument in the RealignerTargetCreator walker. The point here is that the realigner
     * will only proceed with the realignment (even above the given threshold) if it minimizes entropy among the reads (and doesn't simply
     * push the mismatch column to another position). This parameter is just a heuristic and should be adjusted based on your particular data set.
     */
    @Advanced
    @Argument(fullName="entropyThreshold", shortName="entropy", doc="percentage of mismatches at a locus to be considered having high entropy", required=false)
    protected double MISMATCH_THRESHOLD = 0.15;

    /**
     * For expert users only!  To minimize memory consumption you can lower this number (but then the tool may skip realignment on regions with too much coverage;
     * and if the number is too low, it may generate errors during realignment). Just make sure to give Java enough memory! 4Gb should be enough with the default value.
     */
    @Advanced
    @Argument(fullName="maxReadsInMemory", shortName="maxInMemory", doc="max reads allowed to be kept in memory at a time by the SAMFileWriter", required=false)
    protected int MAX_RECORDS_IN_MEMORY = 150000;

    /**
     * For expert users only!
     */
    @Advanced
    @Argument(fullName="maxIsizeForMovement", shortName="maxIsize", doc="maximum insert size of read pairs that we attempt to realign", required=false)
    protected int MAX_ISIZE_FOR_MOVEMENT = 3000;

    /**
     * For expert users only!
     */
    @Advanced
    @Argument(fullName="maxPositionalMoveAllowed", shortName="maxPosMove", doc="maximum positional move in basepairs that a read can be adjusted during realignment", required=false)
    protected int MAX_POS_MOVE_ALLOWED = 200;

    /**
     * For expert users only!  If you need to find the optimal solution regardless of running time, use a higher number.
     */
    @Advanced
    @Argument(fullName="maxConsensuses", shortName="maxConsensuses", doc="max alternate consensuses to try (necessary to improve performance in deep coverage)", required=false)
    protected int MAX_CONSENSUSES = 30;

    /**
     * For expert users only!  If you need to find the optimal solution regardless of running time, use a higher number.
     */
    @Advanced
    @Argument(fullName="maxReadsForConsensuses", shortName="greedy", doc="max reads used for finding the alternate consensuses (necessary to improve performance in deep coverage)", required=false)
    protected int MAX_READS_FOR_CONSENSUSES = 120;

    /**
     * For expert users only!  If this value is exceeded at a given interval, realignment is not attempted and the reads are passed to the output file(s) as-is.
     * If you need to allow more reads (e.g. with very deep coverage) regardless of memory, use a higher number.
     */
    @Advanced
    @Argument(fullName="maxReadsForRealignment", shortName="maxReads", doc="max reads allowed at an interval for realignment", required=false)
    protected int MAX_READS = 20000;

    @Advanced
    @Argument(fullName="noOriginalAlignmentTags", shortName="noTags", required=false, doc="Don't output the original cigar or alignment start tags for each realigned read in the output bam")
    protected boolean NO_ORIGINAL_ALIGNMENT_TAGS = false;

    /**
     * Reads from all input files will be realigned together, but then each read will be saved in the output file corresponding to the input file that
     * the read came from. There are two ways to generate output bam file names: 1) if the value of this argument is a general string (e.g. '.cleaned.bam'),
     * then extensions (".bam" or ".sam") will be stripped from the input file names and the provided string value will be pasted on instead; 2) if the
     * value ends with a '.map' (e.g. input_output.map), then the two-column tab-separated file with the specified name must exist and list unique output
     * file name (2nd column) for each input file name (1st column).
     */
    @Argument(fullName="nWayOut", shortName="nWayOut", required=false, doc="Generate one output file for each input (-I) bam file")
    protected String N_WAY_OUT = null;

    @Hidden
    @Argument(fullName="generate_nWayOut_md5s",doc="Generate md5sums for BAMs")
    protected boolean generateMD5s = false;

    // DEBUGGING OPTIONS FOLLOW

    @Hidden
    @Argument(fullName="check_early",shortName="check_early",required=false,doc="Do early check of reads against existing consensuses")
    protected boolean CHECKEARLY = false;

    @Hidden
    @Argument(fullName="noPGTag", shortName="noPG", required=false,
            doc="Don't output the usual PG tag in the realigned bam file header. FOR DEBUGGING PURPOSES ONLY.  This option is required in order to pass integration tests.")
    protected boolean NO_PG_TAG = false;

    @Hidden
    @Argument(fullName="keepPGTags", shortName="keepPG", required=false,
            doc="Keep older PG tags left in the bam header by previous runs of this tool (by default, all these "+
                    "historical tags will be replaced by the latest tag generated in the current run).")
    protected boolean KEEP_ALL_PG_RECORDS = false;

    @Hidden
    @Output(fullName="indelsFileForDebugging", shortName="indels", required=false, doc="Output file (text) for the indels found; FOR DEBUGGING PURPOSES ONLY")
    protected String OUT_INDELS = null;

    @Hidden
    @Output(fullName="statisticsFileForDebugging", shortName="stats", doc="print out statistics (what does or doesn't get cleaned); FOR DEBUGGING PURPOSES ONLY", required=false)
    protected String OUT_STATS = null;

    @Hidden
    @Output(fullName="SNPsFileForDebugging", shortName="snps", doc="print out whether mismatching columns do or don't get cleaned out; FOR DEBUGGING PURPOSES ONLY", required=false)
    protected String OUT_SNPS = null;

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

    // the intervals input by the user
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;
    private boolean sawReadInCurrentInterval = false;

    // the reads and known indels that fall into the current interval
    private final ReadBin readsToClean = new ReadBin();
    private final ArrayList<GATKSAMRecord> readsNotToClean = new ArrayList<GATKSAMRecord>();
    private final ArrayList<VariantContext> knownIndelsToTry = new ArrayList<VariantContext>();
    private final HashSet<Object> indelRodsSeen = new HashSet<Object>();
    private final HashSet<GATKSAMRecord> readsActuallyCleaned = new HashSet<GATKSAMRecord>();

    private static final int MAX_QUAL = 99;

    // fraction of mismatches that need to no longer mismatch for a column to be considered cleaned
    private static final double MISMATCH_COLUMN_CLEANED_FRACTION = 0.75;

    private static final double SW_MATCH = 30.0;      // 1.0;
    private static final double SW_MISMATCH = -10.0;  //-1.0/3.0;
    private static final double SW_GAP = -10.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -2.0; //-1.0/.0;

    // reference base padding size
    // TODO -- make this a command-line argument if the need arises
    private static final int REFERENCE_PADDING = 30;

    // other output files
    private FileWriter indelOutput = null;
    private FileWriter statsOutput = null;
    private FileWriter snpsOutput = null;

    //###protected Map<SAMReaderID, ConstrainedMateFixingManager> nwayWriters = null;


    // debug info for lazy SW evaluation:
    private long exactMatchesFound = 0; // how many reads exactly matched a consensus we already had
    private long SWalignmentRuns = 0; // how many times (=for how many reads) we ran SW alignment
    private long SWalignmentSuccess = 0; // how many SW alignments were "successful" (i.e. found a workable indel and resulted in non-null consensus)

    private Map<String,String> loadFileNameMap(String mapFile) {
        Map<String,String> fname_map = new HashMap<String,String>();

        try {

             XReadLines reader = new XReadLines(new File(mapFile),true);
             for ( String line : reader ) {
                 if ( line.length() == 0 ) continue;

                 String fields[] = line.split("\t");

                 if ( fields.length != 2 )
                     throw new UserException.BadInput("Input-output map file must have exactly two columns. Offending line:\n"+line);
                 if ( fields[0].length() == 0 || fields[1].length() == 0 )
                     throw new UserException.BadInput("Input-output map file can not have empty strings in either column. Offending line:\n"+line);

                 if ( fname_map.containsKey(fields[0]) )
                     throw new UserException.BadInput("Input-output map file contains duplicate entries for input name "+fields[0]);
                 if ( fname_map.containsValue(fields[1]) )
                     throw new UserException.BadInput("Input-output map file maps multiple entries onto single output name "+fields[1]);

                 fname_map.put(fields[0],fields[1]);
             }
        } catch (IOException e) {
            throw new StingException("I/O Error while reading input-output map file "+N_WAY_OUT+": "+e.getMessage());
        }
       return fname_map;
    }

    public void initialize() {

        if ( N_WAY_OUT == null && writer == null ) {
            throw new UserException.CommandLineException("Either -o or -nWayOut must be specified");
        }
        if ( N_WAY_OUT != null && writer != null ) {
            throw new UserException.CommandLineException("-o and -nWayOut can not be used simultaneously");
        }
        if ( LOD_THRESHOLD < 0.0 )
            throw new RuntimeException("LOD threshold cannot be a negative number");
        if ( MISMATCH_THRESHOLD <= 0.0 || MISMATCH_THRESHOLD > 1.0 )
            throw new RuntimeException("Entropy threshold must be a fraction between 0 and 1");

        try {
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }

        intervals = intervalsFile.getIntervals(getToolkit()).iterator();

        currentInterval = intervals.hasNext() ? intervals.next() : null;

        writerToUse = writer;

        if ( N_WAY_OUT != null ) {
            boolean createIndex =  true;

            if ( N_WAY_OUT.toUpperCase().endsWith(".MAP") ) {
                writerToUse = new NWaySAMFileWriter(getToolkit(),loadFileNameMap(N_WAY_OUT),
                            SAMFileHeader.SortOrder.coordinate,true, createIndex, generateMD5s,createProgramRecord(),KEEP_ALL_PG_RECORDS);
            } else {
                writerToUse = new NWaySAMFileWriter(getToolkit(),N_WAY_OUT,SAMFileHeader.SortOrder.coordinate,true,
                        createIndex, generateMD5s,createProgramRecord(),KEEP_ALL_PG_RECORDS);
            }
        }   else {

            // set up the output writer
            setupWriter(getToolkit().getSAMFileHeader());
        }
        manager = new ConstrainedMateFixingManager(writerToUse, getToolkit().getGenomeLocParser(), MAX_ISIZE_FOR_MOVEMENT, MAX_POS_MOVE_ALLOWED, MAX_RECORDS_IN_MEMORY);

        if ( OUT_INDELS != null ) {
            try {
                indelOutput = new FileWriter(new File(OUT_INDELS));
            } catch (Exception e) {
                logger.error("Failed to create output file "+ OUT_INDELS+". Indel output will be suppressed");
                logger.error(e.getMessage());
                indelOutput = null;
            }
        }
        if ( OUT_STATS != null ) {
            try {
                statsOutput = new FileWriter(new File(OUT_STATS));
            } catch (Exception e) {
                logger.error("Failed to create output file "+ OUT_STATS+". Cleaning stats output will be suppressed");
                logger.error(e.getMessage());
                statsOutput = null;
            }
        }
        if ( OUT_SNPS != null ) {
            try {
                snpsOutput = new FileWriter(new File(OUT_SNPS));
            } catch (Exception e) {
                logger.error("Failed to create output file "+ OUT_SNPS+". Cleaning snps output will be suppressed");
                logger.error(e.getMessage());
                snpsOutput = null;
            }
        }
    }

    private void setupWriter(SAMFileHeader header) {
        
        if ( !NO_PG_TAG ) {
            final SAMProgramRecord programRecord = createProgramRecord();

            List<SAMProgramRecord> oldRecords = header.getProgramRecords();
            List<SAMProgramRecord> newRecords = new ArrayList<SAMProgramRecord>(oldRecords.size()+1);
            for ( SAMProgramRecord record : oldRecords ) {
                if ( !record.getId().startsWith(PROGRAM_RECORD_NAME) || KEEP_ALL_PG_RECORDS )
                    newRecords.add(record);
            }
            newRecords.add(programRecord);
            header.setProgramRecords(newRecords);
        }

        writer.writeHeader(header);
        writer.setPresorted(true);
    }


    private SAMProgramRecord createProgramRecord() {
        if ( NO_PG_TAG ) return null;

        final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
        final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        try {
            final String version = headerInfo.getString("org.broadinstitute.sting.gatk.version");
            programRecord.setProgramVersion(version);
        } catch (MissingResourceException e) {}
        programRecord.setCommandLine(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));
        return programRecord;
    }

    private void emit(final SAMRecord read) {

        // check to see whether the read was modified by looking at the temporary tag
        boolean wasModified = readsActuallyCleaned.contains(read);

        try {
            manager.addRead(read, wasModified);
        } catch (RuntimeIOException e) {
            throw new UserException.ErrorWritingBamFile(e.getMessage());
        }
    }

    private void emitReadLists() {
        // pre-merge lists to sort them in preparation for constrained SAMFileWriter
        readsNotToClean.addAll(readsToClean.getReads());
        ReadUtils.coordinateSortReads(readsNotToClean);
        manager.addReads(readsNotToClean, readsActuallyCleaned);
        readsToClean.clear();
        readsNotToClean.clear();
        readsActuallyCleaned.clear();
    }

    public Integer map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        if ( currentInterval == null ) {
            emit(read);
            return 0;
        }

        // edge case: when the last target interval abuts the end of the genome, we'll get one of the
        //   unmapped reads while the currentInterval still isn't null.  We need to trigger the cleaning
        //   at this point without trying to create a GenomeLoc.
        if ( read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX ) {
            cleanAndCallMap(ref, read, metaDataTracker, null);
            return 0;
        }

        GenomeLoc readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(read);
        // hack to get around unmapped reads having screwy locations
        if ( readLoc.getStop() == 0 )
            readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(readLoc.getContig(), readLoc.getStart(), readLoc.getStart());

        if ( readLoc.isBefore(currentInterval) ) {
            if ( !sawReadInCurrentInterval )
                emit(read);
            else
                readsNotToClean.add(read);
        }
        else if ( readLoc.overlapsP(currentInterval) ) {
            sawReadInCurrentInterval = true;

            if ( doNotTryToClean(read) ) {
                readsNotToClean.add(read);
            } else {
                readsToClean.add(read);

                // add the rods to the list of known variants
                populateKnownIndels(metaDataTracker, ref);
            }

            if ( readsToClean.size() + readsNotToClean.size() >= MAX_READS ) {
                logger.info("Not attempting realignment in interval " + currentInterval + " because there are too many reads.");
                abortCleanForCurrentInterval();
            }
        }
        else {  // the read is past the current interval
            cleanAndCallMap(ref, read, metaDataTracker, readLoc);
        }

        return 0;
    }

    private void abortCleanForCurrentInterval() {
        emitReadLists();
        currentInterval = intervals.hasNext() ? intervals.next() : null;
        sawReadInCurrentInterval = false;
    }

    private boolean doNotTryToClean(SAMRecord read) {
        return read.getReadUnmappedFlag() ||
                read.getNotPrimaryAlignmentFlag() ||
                read.getReadFailsVendorQualityCheckFlag() ||
                read.getMappingQuality() == 0 ||
                read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ||
                ConstrainedMateFixingManager.iSizeTooBigToMove(read, MAX_ISIZE_FOR_MOVEMENT) ||
                ReadUtils.is454Read(read);
        // TODO -- it would be nice if we could use indels from 454 reads as alternate consenses
    }

    private void cleanAndCallMap(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker, GenomeLoc readLoc) {
        if ( readsToClean.size() > 0 ) {
            GenomeLoc earliestPossibleMove = getToolkit().getGenomeLocParser().createGenomeLoc(readsToClean.getReads().get(0));
            if ( manager.canMoveReads(earliestPossibleMove) )
                clean(readsToClean);
        }
        knownIndelsToTry.clear();
        indelRodsSeen.clear();

        emitReadLists();
        try {
            do {
                currentInterval = intervals.hasNext() ? intervals.next() : null;

            } while ( currentInterval != null && (readLoc == null || currentInterval.isBefore(readLoc)) );
        } catch (ReviewedStingException e) {
            throw new UserException.MissortedFile(new File(intervalsFile.getSource()), " *** Are you sure that your interval file is sorted? If not, you must use the --targetIntervalsAreNotSorted argument. ***", e);
        }
        sawReadInCurrentInterval = false;

        // call back into map now that the state has been updated
        map(ref, read, metaDataTracker);
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        if ( readsToClean.size() > 0 ) {
            GenomeLoc earliestPossibleMove = getToolkit().getGenomeLocParser().createGenomeLoc(readsToClean.getReads().get(0));
            if ( manager.canMoveReads(earliestPossibleMove) )
                clean(readsToClean);
            emitReadLists();
        } else if ( readsNotToClean.size() > 0 ) {
            emitReadLists();                            
        }

        knownIndelsToTry.clear();
        indelRodsSeen.clear();

        if ( OUT_INDELS != null ) {
            try {
                indelOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_INDELS+" gracefully. Data may be corrupt.");
            }
        }
        if ( OUT_STATS != null ) {
            try {
                statsOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_STATS+" gracefully. Data may be corrupt.");
            }
        }
        if ( OUT_SNPS != null ) {
            try {
                snpsOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_SNPS+" gracefully. Data may be corrupt.");
            }
        }

        manager.close();
        if ( N_WAY_OUT != null ) writerToUse.close();

        if ( CHECKEARLY ) {
            logger.info("SW alignments runs: "+SWalignmentRuns);
            logger.info("SW alignments successfull: "+SWalignmentSuccess + " ("+SWalignmentSuccess/SWalignmentRuns+"% of SW runs)");
            logger.info("SW alignments skipped (perfect match): "+exactMatchesFound);
            logger.info("Total reads SW worked for: "+(SWalignmentSuccess + exactMatchesFound)+
                    " ("+(SWalignmentSuccess+exactMatchesFound)/(SWalignmentRuns+exactMatchesFound)+"% of all reads requiring SW)");
        }
    }

    private void populateKnownIndels(ReadMetaDataTracker metaDataTracker, ReferenceContext ref) {
        for ( Collection<GATKFeature> rods : metaDataTracker.getContigOffsetMapping().values() ) {
            Iterator<GATKFeature> rodIter = rods.iterator();
            while ( rodIter.hasNext() ) {
                Object rod = rodIter.next().getUnderlyingObject();
                if ( indelRodsSeen.contains(rod) )
                    continue;
                indelRodsSeen.add(rod);
                if ( rod instanceof VariantContext )
                    knownIndelsToTry.add((VariantContext)rod);
            }
        }
    }

    private static int mismatchQualitySumIgnoreCigar(final AlignedRead aRead, final byte[] refSeq, int refIndex, int quitAboveThisValue) {
        final byte[] readSeq = aRead.getReadBases();
        final byte[] quals = aRead.getBaseQualities();
        int sum = 0;
        for (int readIndex = 0 ; readIndex < readSeq.length ; refIndex++, readIndex++ ) {
            if ( refIndex >= refSeq.length ) {
                sum += MAX_QUAL;
                // optimization: once we pass the threshold, stop calculating
                if ( sum > quitAboveThisValue )
                    return sum;
            } else {
                byte refChr = refSeq[refIndex];
                byte readChr = readSeq[readIndex];
                if ( !BaseUtils.isRegularBase(readChr) || !BaseUtils.isRegularBase(refChr) )
                    continue; // do not count Ns/Xs/etc ?
                if ( readChr != refChr ) {
                    sum += (int)quals[readIndex];
                    // optimization: once we pass the threshold, stop calculating
                    if ( sum > quitAboveThisValue )
                        return sum;
                }
            }
        }
        return sum;
    }

    private void clean(ReadBin readsToClean) {

        final List<GATKSAMRecord> reads = readsToClean.getReads();
        if ( reads.size() == 0 )
            return;

        byte[] reference = readsToClean.getReference(referenceReader);
        int leftmostIndex = readsToClean.getLocation().getStart();

        final ArrayList<GATKSAMRecord> refReads = new ArrayList<GATKSAMRecord>();                 // reads that perfectly match ref
        final ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();               // reads that don't perfectly match
        final LinkedList<AlignedRead> altAlignmentsToTest = new LinkedList<AlignedRead>();  // should we try to make an alt consensus from the read?
        final Set<Consensus> altConsenses = new LinkedHashSet<Consensus>();               // list of alt consenses

        // if there are any known indels for this region, get them and create alternate consenses
        generateAlternateConsensesFromKnownIndels(altConsenses, leftmostIndex, reference);

        // decide which reads potentially need to be cleaned;
        // if there are reads with a single indel in them, add that indel to the list of alternate consenses
        long totalRawMismatchSum = determineReadsThatNeedCleaning(reads, refReads, altReads, altAlignmentsToTest, altConsenses, leftmostIndex, reference);

        // use 'Smith-Waterman' to create alternate consenses from reads that mismatch the reference, using totalRawMismatchSum as the random seed
        if ( consensusModel == ConsensusDeterminationModel.USE_SW )
            generateAlternateConsensesFromReads(altAlignmentsToTest, altConsenses, reference, leftmostIndex);

        // if ( debugOn ) System.out.println("------\nChecking consenses...\n--------\n");

        Consensus bestConsensus = null;
        Iterator<Consensus> iter = altConsenses.iterator();

        while ( iter.hasNext() ) {
            Consensus consensus = iter.next();
            //logger.debug("Trying new consensus: " + consensus.cigar + " " + new String(consensus.str));

//            if ( DEBUG ) {
//                System.out.println("Checking consensus with alignment at "+consensus.positionOnReference+" cigar "+consensus.cigar);
//                System.out.println(new String(consensus.str));
//                int z = 0;
//                for ( ; z < consensus.positionOnReference; z++ )  System.out.print('.');
//                for ( z=0 ; z < consensus.cigar.getCigarElement(0).getLength() ; z++ ) System.out.print('.');
//                if ( consensus.cigar.getCigarElement(1).getOperator() == CigarOperator.I ) for ( z= 0; z < consensus.cigar.getCigarElement(1).getLength(); z++ )  System.out.print('I');
//                System.out.println();
//            }

            // if ( debugOn ) System.out.println("Consensus: "+consensus.str);

            for ( int j = 0; j < altReads.size(); j++ ) {
                AlignedRead toTest = altReads.get(j);
                Pair<Integer, Integer> altAlignment = findBestOffset(consensus.str, toTest, leftmostIndex);

                // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                int myScore = altAlignment.second;

                if ( myScore > toTest.getAlignerMismatchScore() || myScore >= toTest.getMismatchScoreToReference() )
                    myScore = toTest.getMismatchScoreToReference();
                // keep track of reads that align better to the alternate consensus.
                // By pushing alignments with equal scores to the alternate, it means we'll over-call (het -> hom non ref) but are less likely to under-call (het -> ref, het non ref -> het)
                else
                    consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.first));

                //logger.debug(consensus.cigar +  " vs. " + toTest.getRead().getReadName() + "-" + toTest.getRead().getReadString() + " => " + myScore + " vs. " + toTest.getMismatchScoreToReference());
                if ( !toTest.getRead().getDuplicateReadFlag() )
                    consensus.mismatchSum += myScore;

                // optimization: once the mismatch sum is higher than the best consensus, quit since this one can't win
                //  THIS MUST BE DISABLED IF WE DECIDE TO ALLOW MORE THAN ONE ALTERNATE CONSENSUS!
                if ( bestConsensus != null && consensus.mismatchSum > bestConsensus.mismatchSum )
                    break;
            }

            //logger.debug("Mismatch sum of new consensus: " + consensus.mismatchSum);
            if ( bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                // we do not need this alt consensus, release memory right away!!
                if ( bestConsensus != null )
                    bestConsensus.readIndexes.clear();
                bestConsensus = consensus;
                //logger.debug("New consensus " + bestConsensus.cigar +  " is now best consensus");
            } else {
                // we do not need this alt consensus, release memory right away!!
                consensus.readIndexes.clear();
            }
        }

        // if:
        // 1) the best alternate consensus has a smaller sum of quality score mismatches than the aligned version of the reads,
        // 2) beats the LOD threshold for the sum of quality score mismatches of the raw version of the reads,
        // 3) didn't just move around the mismatching columns (i.e. it actually reduces entropy), 
        // then clean!
        final double improvement = (bestConsensus == null ? -1 : ((double)(totalRawMismatchSum - bestConsensus.mismatchSum))/10.0);
        if ( improvement >= LOD_THRESHOLD ) {

            bestConsensus.cigar = AlignmentUtils.leftAlignIndel(bestConsensus.cigar, reference, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference);

           // start cleaning the appropriate reads
            for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes ) {
                AlignedRead aRead = altReads.get(indexPair.first);
                if ( !updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.second, aRead, leftmostIndex) )
                    return;
            }
            if ( consensusModel != ConsensusDeterminationModel.KNOWNS_ONLY && !alternateReducesEntropy(altReads, reference, leftmostIndex) ) {
                if ( statsOutput != null ) {
                    try {
                        statsOutput.write(currentInterval.toString());
                        statsOutput.write("\tFAIL (bad indel)\t"); // if improvement > LOD_THRESHOLD *BUT* entropy is not reduced (SNPs still exist)
                        statsOutput.write(Double.toString(improvement));
                        statsOutput.write("\n");
                        statsOutput.flush();
                    } catch (Exception e) {
                        throw new UserException.CouldNotCreateOutputFile("statsOutput", "Failed to write stats output file", e);
                    }
                }
            } else {
                //logger.debug("CLEAN: " + bestConsensus.cigar + " " + bestConsensus.str.toString() + " " + bestConsensus.cigar.numCigarElements() );
                if ( indelOutput != null && bestConsensus.cigar.numCigarElements() > 1 ) {
                    // NOTE: indels are printed out in the format specified for the low-coverage pilot1
                    //  indel calls (tab-delimited): chr position size type sequence
                    StringBuilder str = new StringBuilder();
                    str.append(reads.get(0).getReferenceName());
                    int position = bestConsensus.positionOnReference + bestConsensus.cigar.getCigarElement(0).getLength();
                    str.append("\t" + (leftmostIndex + position - 1));
                    CigarElement ce = bestConsensus.cigar.getCigarElement(1);
                    str.append("\t" + ce.getLength() + "\t" + ce.getOperator() + "\t");
                    int length = ce.getLength();
                    if ( ce.getOperator() == CigarOperator.D ) {
                        for ( int i = 0; i < length; i++)
                            str.append((char)reference[position+i]);
                    } else {
                        for ( int i = 0; i < length; i++)
                            str.append((char)bestConsensus.str[position+i]);
                    }
                    str.append("\t" + (((double)(totalRawMismatchSum - bestConsensus.mismatchSum))/10.0) + "\n");
                    try {
                        indelOutput.write(str.toString());
                        indelOutput.flush();
                    } catch (Exception e) {
                        throw new UserException.CouldNotCreateOutputFile("indelOutput", "Failed to write indel output file", e);
                    }
                }
                if ( statsOutput != null ) {
                    try {
                        statsOutput.write(currentInterval.toString());
                        statsOutput.write("\tCLEAN"); // if improvement > LOD_THRESHOLD *AND* entropy is reduced
                        if ( bestConsensus.cigar.numCigarElements() > 1 )
                            statsOutput.write(" (found indel)");
                        statsOutput.write("\t");
                        statsOutput.write(Double.toString(improvement));
                        statsOutput.write("\n");
                        statsOutput.flush();
                    } catch (Exception e) {
                        throw new UserException.CouldNotCreateOutputFile("statsOutput", "Failed to write stats output file", e);
                    }
                }

                // finish cleaning the appropriate reads
                for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes ) {
                    final AlignedRead aRead = altReads.get(indexPair.first);
                    if ( aRead.finalizeUpdate() ) {
                        // We need to update the mapping quality score of the cleaned reads;
                        // however we don't have enough info to use the proper MAQ scoring system.
                        // For now, we will just arbitrarily add 10 to the mapping quality. [EB, 6/7/2010].
                        // TODO -- we need a better solution here
                        GATKSAMRecord read = aRead.getRead();
                        if ( read.getMappingQuality() != 255 ) // 255 == Unknown, so don't modify it
                            read.setMappingQuality(Math.min(aRead.getRead().getMappingQuality() + 10, 254));

                        // before we fix the attribute tags we first need to make sure we have enough of the reference sequence
                        int neededBasesToLeft = leftmostIndex - read.getAlignmentStart();
                        int neededBasesToRight = read.getAlignmentEnd() - leftmostIndex - reference.length + 1;
                        int neededBases = Math.max(neededBasesToLeft, neededBasesToRight);
                        if ( neededBases > 0 ) {
                            int padLeft = Math.max(leftmostIndex-neededBases, 1);
                            int padRight = Math.min(leftmostIndex+reference.length+neededBases, referenceReader.getSequenceDictionary().getSequence(currentInterval.getContig()).getSequenceLength());
                            reference = referenceReader.getSubsequenceAt(currentInterval.getContig(), padLeft, padRight).getBases();
                            leftmostIndex = padLeft;
                        }

                        // now, fix the attribute tags
                        // TODO -- get rid of this try block when Picard does the right thing for reads aligned off the end of the reference
                        try {
                            if ( read.getAttribute(SAMTag.NM.name()) != null )
                                read.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(read, reference, leftmostIndex-1));
                            if ( read.getAttribute(SAMTag.UQ.name()) != null )
                                read.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(read, reference, leftmostIndex-1));
                        } catch (Exception e) {
                            // ignore it
                        }
                        // TODO -- this is only temporary until Tim adds code to recalculate this value
                        if ( read.getAttribute(SAMTag.MD.name()) != null )
                            read.setAttribute(SAMTag.MD.name(), null);

                        // mark that it was actually cleaned
                        readsActuallyCleaned.add(read);
                    }
                }
            }

            // END IF ( improvement >= LOD_THRESHOLD )

        } else if ( statsOutput != null ) {
            try {
                statsOutput.write(String.format("%s\tFAIL\t%.1f%n",
                        currentInterval.toString(), improvement));
                statsOutput.flush();
            } catch (Exception e) {
                throw new UserException.CouldNotCreateOutputFile("statsOutput", "Failed to write stats output file", e);
            }
        }
    }

    private void generateAlternateConsensesFromKnownIndels(final Set<Consensus> altConsensesToPopulate, final int leftmostIndex, final byte[] reference) {
        for ( VariantContext knownIndel : knownIndelsToTry ) {
            if ( knownIndel == null || !knownIndel.isIndel() || knownIndel.isComplexIndel() )
                continue;
            byte[] indelStr = knownIndel.isSimpleInsertion() ? knownIndel.getAlternateAllele(0).getBases() : Utils.dupBytes((byte)'-', knownIndel.getReference().length());
            int start = knownIndel.getStart() - leftmostIndex + 1;
            Consensus c = createAlternateConsensus(start, reference, indelStr, knownIndel);
            if ( c != null )
                altConsensesToPopulate.add(c);
        }
    }

    private long determineReadsThatNeedCleaning(final List<GATKSAMRecord> reads,
                                                final ArrayList<GATKSAMRecord> refReadsToPopulate,
                                                final ArrayList<AlignedRead> altReadsToPopulate,
                                                final LinkedList<AlignedRead> altAlignmentsToTest,
                                                final Set<Consensus> altConsenses,
                                                final int leftmostIndex,
                                                final byte[] reference) {

        long totalRawMismatchSum = 0L;

        for ( final GATKSAMRecord read : reads ) {

            // we can not deal with screwy records
            if ( read.getCigar().numCigarElements() == 0 ) {
                refReadsToPopulate.add(read);
                continue;
            }

            final AlignedRead aRead = new AlignedRead(read);

            // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
            int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
            if ( numBlocks == 2 ) {
                Cigar newCigar = AlignmentUtils.leftAlignIndel(unclipCigar(read.getCigar()), reference, read.getReadBases(), read.getAlignmentStart()-leftmostIndex, 0);
                aRead.setCigar(newCigar, false);
            }

            final int startOnRef = read.getAlignmentStart()-leftmostIndex;
            final int rawMismatchScore = mismatchQualitySumIgnoreCigar(aRead, reference, startOnRef, Integer.MAX_VALUE);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( rawMismatchScore > 0 ) {
                altReadsToPopulate.add(aRead);
                //logger.debug("Adding " + read.getReadName() + " with raw mismatch score " + rawMismatchScore + " to non-ref reads");

                if ( !read.getDuplicateReadFlag() )
                    totalRawMismatchSum += rawMismatchScore;
                aRead.setMismatchScoreToReference(rawMismatchScore);
                aRead.setAlignerMismatchScore(AlignmentUtils.mismatchingQualities(aRead.getRead(), reference, startOnRef));

                // if it has an indel, let's see if that's the best consensus
                if ( consensusModel != ConsensusDeterminationModel.KNOWNS_ONLY && numBlocks == 2 )  {
                    Consensus c = createAlternateConsensus(startOnRef, aRead.getCigar(), reference, aRead.getReadBases());
                    if ( c != null )
                        altConsenses.add(c);
                } else {
                    altAlignmentsToTest.add(aRead);
                }
            }
            // otherwise, we can emit it as is
            else {
                //logger.debug("Adding " + read.getReadName() + " with raw mismatch score " + rawMismatchScore + " to ref reads");
                refReadsToPopulate.add(read);
            }
        }

        return totalRawMismatchSum;
    }

    private void generateAlternateConsensesFromReads(final LinkedList<AlignedRead> altAlignmentsToTest,
                                                     final Set<Consensus> altConsensesToPopulate,
                                                     final byte[] reference,
                                                     final int leftmostIndex) {

        // if we are under the limit, use all reads to generate alternate consenses
        if ( altAlignmentsToTest.size() <= MAX_READS_FOR_CONSENSUSES ) {
            for ( AlignedRead aRead : altAlignmentsToTest ) {
                if ( CHECKEARLY ) createAndAddAlternateConsensus1(aRead, altConsensesToPopulate, reference,leftmostIndex);
                else createAndAddAlternateConsensus(aRead.getReadBases(), altConsensesToPopulate, reference);
            }
        }
        // otherwise, choose reads for alternate consenses randomly
        else {
            int readsSeen = 0;
            while ( readsSeen++ < MAX_READS_FOR_CONSENSUSES && altConsensesToPopulate.size() <= MAX_CONSENSUSES) {
                int index = GenomeAnalysisEngine.getRandomGenerator().nextInt(altAlignmentsToTest.size());
                AlignedRead aRead = altAlignmentsToTest.remove(index);
                if ( CHECKEARLY ) createAndAddAlternateConsensus1(aRead, altConsensesToPopulate, reference,leftmostIndex);
                else createAndAddAlternateConsensus(aRead.getReadBases(), altConsensesToPopulate, reference);
            }
        }
    }

    private void createAndAddAlternateConsensus(final byte[] read, final Set<Consensus> altConsensesToPopulate, final byte[] reference) {

        // do a pairwise alignment against the reference
         SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, read, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
         Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, read);
         if ( c != null )
             altConsensesToPopulate.add(c);
    }

    private void createAndAddAlternateConsensus1(AlignedRead read, final Set<Consensus> altConsensesToPopulate,
                                                 final byte[] reference, final int leftmostIndex) {

         for ( Consensus known : altConsensesToPopulate ) {
              Pair<Integer, Integer> altAlignment = findBestOffset(known.str, read, leftmostIndex);
              // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
              int myScore = altAlignment.second;
              if ( myScore == 0 ) {exactMatchesFound++; return; }// read matches perfectly to a known alt consensus - no need to run SW, we already know the answer
         }
         // do a pairwise alignment against the reference
         SWalignmentRuns++;
         SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, read.getReadBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
         Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, read.getReadBases());
         if ( c != null ) {
             altConsensesToPopulate.add(c);
             SWalignmentSuccess++;
         }
    }

    // create a Consensus from cigar/read strings which originate somewhere on the reference
    private Consensus createAlternateConsensus(final int indexOnRef, final Cigar c, final byte[] reference, final byte[] readStr) {
        if ( indexOnRef < 0 )
            return null;

        // if there are no indels, we do not need this consensus, can abort early:
        if ( c.numCigarElements() == 1 && c.getCigarElement(0).getOperator() == CigarOperator.M ) return null;

        // create the new consensus
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(c.numCigarElements()-1);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < indexOnRef; i++)
            sb.append((char)reference[i]);

        int indelCount = 0;
        int altIdx = 0;
        int refIdx = indexOnRef;
        boolean ok_flag = true;
        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
            CigarElement ce = c.getCigarElement(i);
            int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
            case D:
                refIdx += elementLength;
                indelCount++;
                elements.add(ce);
                break;
            case M:
                altIdx += elementLength;
            case N:
                if ( reference.length < refIdx + elementLength )
                    ok_flag = false;
                else  {
                    for (int j = 0; j < elementLength; j++)
                        sb.append((char)reference[refIdx+j]);
                }
                refIdx += elementLength;
                elements.add(new CigarElement(elementLength, CigarOperator.M));
                break;
            case I:
                for (int j = 0; j < elementLength; j++) {
                    if ( ! BaseUtils.isRegularBase(readStr[altIdx+j]) ) {
                        // Insertions with N's in them cause real problems sometimes; it's better to drop them altogether
                        ok_flag=false;
                        break;
                    }
                    sb.append((char)readStr[altIdx + j]);
                }
                altIdx += elementLength;
                indelCount++;
                elements.add(ce);
                break;
            case S:
            default:
                break;
            }
        }
        // make sure that there is at most only a single indel and it aligns appropriately!
        if ( !ok_flag || indelCount != 1 || reference.length < refIdx )
            return null;

        for (int i = refIdx; i < reference.length; i++)
            sb.append((char)reference[i]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, new Cigar(elements), indexOnRef);
    }

    // create a Consensus from just the indel string that falls on the reference
    private Consensus createAlternateConsensus(final int indexOnRef, final byte[] reference, final byte[] indelStr, final VariantContext indel) {
        if ( indexOnRef < 0 || indexOnRef >= reference.length )
            return null;

        // create the new consensus
        StringBuilder sb = new StringBuilder();
        Cigar cigar = new Cigar();
        int refIdx;

        for (refIdx = 0; refIdx < indexOnRef; refIdx++)
            sb.append((char)reference[refIdx]);
        if ( indexOnRef > 0 )
            cigar.add(new CigarElement(indexOnRef, CigarOperator.M));

        if ( indel.isSimpleDeletion() ) {
            refIdx += indelStr.length;
            cigar.add(new CigarElement(indelStr.length, CigarOperator.D));
        }
        else if ( indel.isSimpleInsertion() ) {
            for ( byte b : indelStr )
                sb.append((char)b);
            cigar.add(new CigarElement(indelStr.length, CigarOperator.I));
        } else {
            throw new IllegalStateException("Creating an alternate consensus from a complex indel is not allows");
        }

        if ( reference.length - refIdx > 0 )
            cigar.add(new CigarElement(reference.length - refIdx, CigarOperator.M));        
        for (; refIdx < reference.length; refIdx++)
            sb.append((char)reference[refIdx]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, cigar, 0);
    }

    private Pair<Integer, Integer> findBestOffset(final byte[] ref, final AlignedRead read, final int leftmostIndex) {

        // optimization: try the most likely alignment first (to get a low score to beat)
        int originalAlignment = read.getOriginalAlignmentStart() - leftmostIndex;
        int bestScore = mismatchQualitySumIgnoreCigar(read, ref, originalAlignment, Integer.MAX_VALUE);
        int bestIndex = originalAlignment;

        // optimization: we can't get better than 0, so we can quit now
        if ( bestScore == 0 )
            return new Pair<Integer, Integer>(bestIndex, 0);

        // optimization: the correct alignment shouldn't be too far from the original one (or else the read wouldn't have aligned in the first place)
        for ( int i = 0; i < originalAlignment; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        final int maxPossibleStart = ref.length - read.getReadLength();
        for ( int i = originalAlignment + 1; i <= maxPossibleStart; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        return new Pair<Integer, Integer>(bestIndex, bestScore);
    }


    private boolean updateRead(final Cigar altCigar, final int altPosOnRef, final int myPosOnAlt, final AlignedRead aRead, final int leftmostIndex) {
        Cigar readCigar = new Cigar();

        // special case: there is no indel
        if ( altCigar.getCigarElements().size() == 1 ) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.setCigar(readCigar);
            return true;
        }

        CigarElement altCE1 = altCigar.getCigarElement(0);
        CigarElement altCE2 = altCigar.getCigarElement(1);

        int leadingMatchingBlockLength = 0; // length of the leading M element or 0 if the leading element is I

        CigarElement indelCE;
        if ( altCE1.getOperator() == CigarOperator.I  ) {
            indelCE=altCE1;
            if ( altCE2.getOperator() != CigarOperator.M  ) {
                logger.warn("When the first element of the alt consensus is I, the second one must be M. Actual: " + altCigar.toString() + ".  Skipping this site...");
                return false;
            }
        }
        else {
            if ( altCE1.getOperator() != CigarOperator.M  ) {
                logger.warn("First element of the alt consensus cigar must be M or I. Actual: " + altCigar.toString() + ".  Skipping this site...");
                return false;
            }
            if ( altCE2.getOperator() == CigarOperator.I  || altCE2.getOperator() == CigarOperator.D ) {
                indelCE=altCE2;
            } else {
                logger.warn("When first element of the alt consensus is M, the second one must be I or D. Actual: " + altCigar.toString() + ".  Skipping this site...");
                return false;
            }
            leadingMatchingBlockLength = altCE1.getLength();
        }

        // the easiest thing to do is to take each case separately
        int endOfFirstBlock = altPosOnRef + leadingMatchingBlockLength;
        boolean sawAlignmentStart = false;

        // for reads starting before the indel
        if ( myPosOnAlt < endOfFirstBlock) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            sawAlignmentStart = true;

            // for reads ending before the indel
            if ( myPosOnAlt + aRead.getReadLength() <= endOfFirstBlock) {
                //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
                //aRead.setCigar(readCigar);
                aRead.setCigar(null); // reset to original alignment
                return true;
            }
            readCigar.add(new CigarElement(endOfFirstBlock - myPosOnAlt, CigarOperator.M));
        }

        // forward along the indel
        //int indelOffsetOnRef = 0, indelOffsetOnRead = 0;
        if ( indelCE.getOperator() == CigarOperator.I ) {
            // for reads that end in an insertion
            if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.getLength() ) {
                int partialInsertionLength = myPosOnAlt + aRead.getReadLength() - endOfFirstBlock;
                // if we also started inside the insertion, then we need to modify the length
                if ( !sawAlignmentStart )
                    partialInsertionLength = aRead.getReadLength();
                readCigar.add(new CigarElement(partialInsertionLength, CigarOperator.I));
                aRead.setCigar(readCigar);
                return true;
            }

            // for reads that start in an insertion
            if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.getLength() ) {
                aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(indelCE.getLength() - (myPosOnAlt - endOfFirstBlock), CigarOperator.I));
                //indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
                sawAlignmentStart = true;
            } else if ( sawAlignmentStart ) {
                readCigar.add(indelCE);
                //indelOffsetOnRead = indelCE.getLength();
            }
        } else if ( indelCE.getOperator() == CigarOperator.D ) {
            if ( sawAlignmentStart )
                readCigar.add(indelCE);
            //indelOffsetOnRef = indelCE.getLength();
        }

        // for reads that start after the indel
        if ( !sawAlignmentStart ) {
            //aRead.setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
            //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            //aRead.setCigar(readCigar);
            aRead.setCigar(null); // reset to original alignment
            return true;
        }

        int readRemaining = aRead.getReadBases().length;
        for ( CigarElement ce : readCigar.getCigarElements() ) {
            if ( ce.getOperator() != CigarOperator.D )
                readRemaining -= ce.getLength();
        }
        if ( readRemaining > 0 )
            readCigar.add(new CigarElement(readRemaining, CigarOperator.M));
        aRead.setCigar(readCigar);

        return true;
    }

    private boolean alternateReducesEntropy(final List<AlignedRead> reads, final byte[] reference, final int leftmostIndex) {
        final int[] originalMismatchBases = new int[reference.length];
        final int[] cleanedMismatchBases = new int[reference.length];
        final int[] totalOriginalBases = new int[reference.length];
        final int[] totalCleanedBases = new int[reference.length];

        // set to 1 to prevent dividing by zero
        for ( int i=0; i < reference.length; i++ )
            originalMismatchBases[i] = totalOriginalBases[i] = cleanedMismatchBases[i] = totalCleanedBases[i] = 0;

        for (int i=0; i < reads.size(); i++) {
            final AlignedRead read = reads.get(i);
            if ( read.getRead().getAlignmentBlocks().size() > 1 )
                 continue;

            int refIdx = read.getOriginalAlignmentStart() - leftmostIndex;
            final byte[] readStr = read.getReadBases();
            final byte[] quals = read.getBaseQualities();

            for (int j=0; j < readStr.length; j++, refIdx++ ) {
                if ( refIdx < 0 || refIdx >= reference.length ) {
                    //System.out.println( "Read: "+read.getRead().getReadName() + "; length = " + readStr.length() );
                    //System.out.println( "Ref left: "+ leftmostIndex +"; ref length=" + reference.length() + "; read alignment start: "+read.getOriginalAlignmentStart() );
                    break;
                }
                totalOriginalBases[refIdx] += quals[j];
                if ( readStr[j] != reference[refIdx] )
                    originalMismatchBases[refIdx] += quals[j];
            }

            // reset and now do the calculation based on the cleaning
            refIdx = read.getAlignmentStart() - leftmostIndex;
            int altIdx = 0;
            Cigar c = read.getCigar();
            for (int j = 0 ; j < c.numCigarElements() ; j++) {
                CigarElement ce = c.getCigarElement(j);
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case M:
                        for (int k = 0 ; k < elementLength ; k++, refIdx++, altIdx++ ) {
                            if ( refIdx >= reference.length )
                                break;
                            totalCleanedBases[refIdx] += quals[altIdx];
                            if ( readStr[altIdx] != reference[refIdx] )
                                cleanedMismatchBases[refIdx] += quals[altIdx];
                        }
                        break;
                    case I:
                        altIdx += elementLength;
                        break;
                    case D:
                        refIdx += elementLength;
                        break;
                    case S:
                    default:
                        break;
                }
            }
        }

        int originalMismatchColumns = 0, cleanedMismatchColumns = 0;
        StringBuilder sb = new StringBuilder();
        for ( int i=0; i < reference.length; i++ ) {
            if ( cleanedMismatchBases[i] == originalMismatchBases[i] )
                continue;
            boolean didMismatch = false, stillMismatches = false;
            if ( originalMismatchBases[i] > totalOriginalBases[i] * MISMATCH_THRESHOLD )  {
                didMismatch = true;
                originalMismatchColumns++;
                if ( totalCleanedBases[i] > 0 && ((double)cleanedMismatchBases[i] / (double)totalCleanedBases[i]) > ((double)originalMismatchBases[i] / (double)totalOriginalBases[i]) * (1.0 - MISMATCH_COLUMN_CLEANED_FRACTION) ) {
                    stillMismatches = true;
                    cleanedMismatchColumns++;
                }
            } else if ( cleanedMismatchBases[i] > totalCleanedBases[i] * MISMATCH_THRESHOLD ) {
                cleanedMismatchColumns++;
            }
            if ( snpsOutput != null ) {
                    if ( didMismatch ) {
                        sb.append(reads.get(0).getRead().getReferenceName() + ":");
                        sb.append((leftmostIndex + i));
                        if ( stillMismatches )
                            sb.append(" SAME_SNP\n");
                        else
                            sb.append(" NOT_SNP\n");
                    }
            }
        }

        //logger.debug("Original mismatch columns = " + originalMismatchColumns + "; cleaned mismatch columns = " + cleanedMismatchColumns);

        final boolean reduces = (originalMismatchColumns == 0 || cleanedMismatchColumns < originalMismatchColumns);
        if ( reduces && snpsOutput != null ) {
            try {
                snpsOutput.write(sb.toString());
                snpsOutput.flush();
            } catch (Exception e) {
                throw new UserException.CouldNotCreateOutputFile("snpsOutput", "Failed to write SNPs output file", e);
            }
        }
        return reduces;
    }

    protected static Cigar unclipCigar(Cigar cigar) {
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(cigar.numCigarElements());
        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) )
                elements.add(ce);
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    protected static Cigar reclipCigar(Cigar cigar, SAMRecord read) {
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>();

        int i = 0;
        int n = read.getCigar().numCigarElements();
        while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            elements.add(read.getCigar().getCigarElement(i++));

        elements.addAll(cigar.getCigarElements());

        i++;
        while ( i < n && !isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            i++;

        while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
            elements.add(read.getCigar().getCigarElement(i++));

        return new Cigar(elements);
    }

    private class AlignedRead {
        private final GATKSAMRecord read;
        private byte[] readBases = null;
        private byte[] baseQuals = null;
        private Cigar newCigar = null;
        private int newStart = -1;
        private int mismatchScoreToReference = 0;
        private long alignerMismatchScore = 0;

        public AlignedRead(GATKSAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public GATKSAMRecord getRead() {
               return read;
        }

        public int getReadLength() {
            return readBases != null ? readBases.length : read.getReadLength();
        }

        public byte[] getReadBases() {
            if ( readBases == null )
                getUnclippedBases();
            return readBases;
        }

        public byte[] getBaseQualities() {
            if ( baseQuals == null )
                getUnclippedBases();
            return baseQuals;
        }

        // pull out the bases that aren't clipped out
        private void getUnclippedBases() {
            readBases = new byte[getReadLength()];
            baseQuals = new byte[getReadLength()];
            byte[] actualReadBases = read.getReadBases();
            byte[] actualBaseQuals = read.getBaseQualities();
            int fromIndex = 0, toIndex = 0;

            for ( CigarElement ce : read.getCigar().getCigarElements() ) {
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case S:
                        fromIndex += elementLength;
                        break;
                    case M:
                    case I:
                        System.arraycopy(actualReadBases, fromIndex, readBases, toIndex, elementLength);
                        System.arraycopy(actualBaseQuals, fromIndex, baseQuals, toIndex, elementLength);
                        fromIndex += elementLength;
                        toIndex += elementLength;
                    default:
                        break;
                }
            }

            // if we got clipped, trim the array
            if ( fromIndex != toIndex ) {
                byte[] trimmedRB = new byte[toIndex];
                byte[] trimmedBQ = new byte[toIndex];
                System.arraycopy(readBases, 0, trimmedRB, 0, toIndex);
                System.arraycopy(baseQuals, 0, trimmedBQ, 0, toIndex);
                readBases = trimmedRB;
                baseQuals = trimmedBQ;
            }
        }

        public Cigar getCigar() {
            return (newCigar != null ? newCigar : read.getCigar());
        }

        public void setCigar(Cigar cigar) {
            setCigar(cigar, true);
        }

        // tentatively sets the new Cigar, but it needs to be confirmed later
        public void setCigar(Cigar cigar, boolean fixClippedCigar) {
            if ( cigar == null ) {
                newCigar = null;
                return;
            }

            if ( fixClippedCigar && getReadBases().length < read.getReadLength() )
                cigar = reclipCigar(cigar);

            // no change?
            if ( read.getCigar().equals(cigar) ) {
                newCigar = null;
                return;
            }

            // no indel?
            String str = cigar.toString();
            if ( !str.contains("D") && !str.contains("I") ) {
                logger.debug("Modifying a read with no associated indel; although this is possible, it is highly unlikely.  Perhaps this region should be double-checked: " + read.getReadName() + " near " + read.getReferenceName() + ":" + read.getAlignmentStart());
                //    newCigar = null;
                //    return;
            }

            newCigar = cigar;
        }

        // pull out the bases that aren't clipped out
        private Cigar reclipCigar(Cigar cigar) {
            return IndelRealigner.reclipCigar(cigar, read);
        }

        // tentatively sets the new start, but it needs to be confirmed later
        public void setAlignmentStart(int start) {
            newStart = start;
        }

        public int getAlignmentStart() {
            return (newStart != -1 ? newStart : read.getAlignmentStart());
        }

        public int getOriginalAlignmentStart() {
            return read.getAlignmentStart();
        }

        // finalizes the changes made.
        // returns true if this record actually changes, false otherwise
        public boolean finalizeUpdate() {
            // if we haven't made any changes, don't do anything
            if ( newCigar == null )
                return false;
            if ( newStart == -1 )
                newStart = read.getAlignmentStart();
            else if ( Math.abs(newStart - read.getAlignmentStart()) > MAX_POS_MOVE_ALLOWED ) {
                logger.debug(String.format("Attempting to realign read %s at %d more than %d bases to %d.", read.getReadName(), read.getAlignmentStart(), MAX_POS_MOVE_ALLOWED, newStart));
                return false;
            }

            // annotate the record with the original cigar (and optionally the alignment start)
            if ( !NO_ORIGINAL_ALIGNMENT_TAGS ) {
                read.setAttribute(ORIGINAL_CIGAR_TAG, read.getCigar().toString());
                if ( newStart != read.getAlignmentStart() )
                    read.setAttribute(ORIGINAL_POSITION_TAG, read.getAlignmentStart());
            }

            read.setCigar(newCigar);
            read.setAlignmentStart(newStart);

            return true;
        }

        public void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }

        public int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }

        public void setAlignerMismatchScore(long score) {
            alignerMismatchScore = score;
        }

        public long getAlignerMismatchScore() {
            return alignerMismatchScore;
        }
    }

    private static class Consensus {
        public final byte[] str;
        public final ArrayList<Pair<Integer, Integer>> readIndexes;
        public final int positionOnReference;
        public int mismatchSum;
        public Cigar cigar;

        public Consensus(byte[] str, Cigar cigar, int positionOnReference) {
            this.str = str;
            this.cigar = cigar;
            this.positionOnReference = positionOnReference;
            mismatchSum = 0;
            readIndexes = new ArrayList<Pair<Integer, Integer>>();
        }

        @Override
        public boolean equals(Object o) {
            return ( this == o || (o instanceof Consensus && Arrays.equals(this.str,(((Consensus)o).str)) ) );
        }

        public boolean equals(Consensus c) {
            return ( this == c || Arrays.equals(this.str,c.str) ) ;
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(this.str);
        }
    }

    private class ReadBin implements HasGenomeLocation {

        private final ArrayList<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        private byte[] reference = null;
        private GenomeLoc loc = null;

        public ReadBin() { }

        // Return false if we can't process this read bin because the reads are not correctly overlapping.
        // This can happen if e.g. there's a large known indel with no overlapping reads.
        public void add(GATKSAMRecord read) {

            GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(read);
            if ( loc == null )
                loc = locForRead;
            else if ( locForRead.getStop() > loc.getStop() )
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());

            reads.add(read);
        }

        public List<GATKSAMRecord> getReads() { return reads; }

        public byte[] getReference(IndexedFastaSequenceFile referenceReader) {
            // set up the reference if we haven't done so yet
            if ( reference == null ) {
                // first, pad the reference to handle deletions in narrow windows (e.g. those with only 1 read)
                int padLeft = Math.max(loc.getStart()-REFERENCE_PADDING, 1);
                int padRight = Math.min(loc.getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(loc.getContig()).getSequenceLength());
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), padLeft, padRight);
                reference = referenceReader.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
                StringUtil.toUpperCase(reference);
            }

            return reference;
        }

        public GenomeLoc getLocation() { return loc; }

        public int size() { return reads.size(); }

        public void clear() {
            reads.clear();
            reference = null;
            loc = null;
        }

    }
}
