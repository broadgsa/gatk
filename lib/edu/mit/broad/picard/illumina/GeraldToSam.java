/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Iterator;

import edu.mit.broad.picard.util.*;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.CommandLineParser;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMFileWriter;
import edu.mit.broad.sam.SAMFileWriterFactory;
import edu.mit.broad.sam.SAMProgramRecord;
import edu.mit.broad.sam.SAMReadGroupRecord;
import edu.mit.broad.sam.SAMRecord;

/**
 * Read alignments for a lane (paired or unpaired) from Gerald directory and write to SAM file.
 */
public class GeraldToSam extends CommandLineProgram {

    // These are all written to the SAM header
    private static final String DEFAULT_CN = "broad";
    private static final String DEFAULT_PL = "illumina";
    private static final String PROGRAM_VERSION = "1.0";
    private static final String READ_GROUP_ID = "0";
    private static final String PROGRAM_RECORD_ID = "0";
    private static final String UNKNOWN_SAMPLE = "N/A";

    private static final Log log = Log.getInstance(GeraldToSam.class);

    // The following attributes define the command-line arguments
    @Usage(programVersion=PROGRAM_VERSION)
    public String USAGE =
            getStandardUsagePreamble() +
                    "Read Gerald alignments for the given lane, and write in SAM format, coordinate sorted.\n";

    @Option(shortName = "G", doc = "Location of Gerald files.")
    public File GERALD_DIR;

    @Option(shortName = "L")
    public Integer LANE;

    @Option(shortName = "M", doc = "Translates from Gerald alignment coordinates to genomic coordinates.")
    public File SQUASHED_MAP;

    @Option(shortName = "D", doc = "Input SAM or BAM file defining the names, sizes and order of the reference contig, " +
                                    "and other reference metadata.")
    public File SEQUENCE_DICT;

    @Option(shortName = "O", doc = "SAM or BAM file to be written (file extension determines format).")
    public File OUTPUT;

    @Option(doc = "Populates SM field of read group.  Use pool name when a pool is being sequenced.  " +
            "If any other read group fields are specified, then this is required.")
    public String SAMPLE = UNKNOWN_SAMPLE;

    @Option(doc = "Populates LB field of read group.")
    public String LIBRARY;

    @Option(doc = "Populates DS field of read group.", optional = true)
    public String DESCRIPTION;

    @Option(doc = "Flowcell.lane.  Populates PU field of read group.")
    public String RUN;

    @Option(doc = "Predicted median insert size (may be different from the actual median insert size.  " +
                    "Populates the PI field of read group.", optional = true)
    public Integer PI;

    @Option(doc = "Sequencing center that produced the reads.  Populates CN field of read group.")
    public String CN = DEFAULT_CN;

    @Option(doc = "Date the run was produced.  Populates the DT field of read group.")
    public Date RUN_DATE;

    @Option(doc = "Platform/technology used to produce the reads.  Populates the PL field of read group")
    public String PL = DEFAULT_PL;

    @Option(shortName = "JUMPING", doc = "True if this is a jumping library")
    public Boolean JUMPING_LIBRARY = Boolean.FALSE;

    @Option(doc = "String to put in the PG:CL header field.  If not present, the GeraldToSam command line is put there",
        optional = true)
    public String ALIGNMENT_COMMAND;

    @Option(doc = "Write no more than this number of alignment records.  Default: Write all the alignment records",
    optional = true)
    public Integer MAX_ALIGNMENTS;

    private SAMFileWriter writer;
    SAMFileHeader header;
    private boolean paired;


    public static void main(final String[] argv) {
        System.exit(new GeraldToSam().instanceMain(argv));
    }

    @Override
	public int doWork() {
        makeHeader(clp.getArgv());
        writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);
        writeAlignments();
        writer.close();
        return 0;
    }

    /**
     * If any of the read group options are specified on the command line, then SAMPLE must be specified.
     * This is currently not doing anything because SAMPLE has a non-null default value.
     * @return false if there is a problem with the command line
     */
    @Override
	protected boolean customCommandLineValidation() {
        if (SAMPLE == null &&
                (LIBRARY != null || DESCRIPTION != null || RUN != null || PI != null || !CN.equals(DEFAULT_CN)
                        || RUN_DATE != null || !PL.equals(DEFAULT_PL)
                )) {
            System.err.println("SAMPLE must be specified if any read group options are used.");
            clp.usage(System.err);
            return false;
        }
        return true;
    }

    /**
     * Create the SAMFileHeader given the cmd-line args
     * @param argv
     */
    private void makeHeader(final String[] argv) {
        header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_ID);
        programRecord.setProgramVersion(PROGRAM_VERSION);
        String commandLine = ALIGNMENT_COMMAND;
        if (commandLine == null) {
            commandLine = StringUtil.join(" ", argv);
        }
        programRecord.setCommandLine(commandLine);
        header.addProgramRecord(programRecord);

        final SAMFileReader sequenceDictionary = new SAMFileReader(SEQUENCE_DICT);
        final SAMFileHeader sequenceDictionaryHeader = sequenceDictionary.getFileHeader();
        header.setSequences(sequenceDictionaryHeader.getSequences());

        if (SAMPLE != null) {
            final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(READ_GROUP_ID);
            final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
            readGroups.add(readGroup);
            readGroup.setSample(SAMPLE);
            if (LIBRARY != null) {
                readGroup.setLibrary(LIBRARY);
            }
            setRGAttributeIfNotNull(readGroup, DESCRIPTION, "DS");
            setRGAttributeIfNotNull(readGroup, RUN, "PU");
            setRGAttributeIfNotNull(readGroup, PI, SAMReadGroupRecord.PREDICTED_MEDIAN_INSERT_SIZE_TAG);
            setRGAttributeIfNotNull(readGroup, CN, "CN");
            setRGAttributeIfNotNull(readGroup, RUN_DATE, SAMReadGroupRecord.DATE_RUN_PRODUCED_TAG);
            setRGAttributeIfNotNull(readGroup, PL, "PL");
            header.setReadGroups(readGroups);
        }
    }

    private void setRGAttributeIfNotNull(final SAMReadGroupRecord readGroup, final Object value, final String key) {
        if (value == null) {
            return;
        }
        readGroup.setAttribute(key, value);
    }

    /**
     * Iterate through the Gerald output and write alignments.  eland_extended.txt and export.txt are
     * iterated together using PasteParser.  If paired end lane, then two PasteParsers are iterated in tandem,
     * so that mate info is available when a SAMRecord is created.
     */
    private void writeAlignments() {
        final GeraldParserFactory geraldParserFactory = new GeraldParserFactory(GERALD_DIR, LANE, SQUASHED_MAP);
        paired = geraldParserFactory.isPairedRun();
        final GeraldParser firstEndIterator = geraldParserFactory.makeParser(paired ? 1: null);
        GeraldParser secondEndIterator = null;
        if (paired) {
            secondEndIterator = geraldParserFactory.makeParser(2);
        }
        int numAlignmentsOrPairsWritten = 0;
        while (firstEndIterator.hasNext()) {
            final GeraldParser.GeraldAlignment firstEnd = firstEndIterator.next();
            GeraldParser.GeraldAlignment secondEnd = null;
            if (paired) {
                hasNextAssert(secondEndIterator);
                secondEnd = secondEndIterator.next();
            }
            final SAMRecord firstEndAlignment = createSAMRecordFromGerald(firstEnd);
            SAMRecord secondEndAlignment = null;
            if (paired) {
                secondEndAlignment = createSAMRecordFromGerald(secondEnd);
                setMateInfo(secondEndAlignment, firstEnd);
                setMateInfo(firstEndAlignment, secondEnd);
                secondEndAlignment.setSecondOfPairFlag(true);
                firstEndAlignment.setFirstOfPairFlag(true);
                final boolean properPair = SamPairUtil.isProperPair(firstEndAlignment, secondEndAlignment, JUMPING_LIBRARY);
                firstEndAlignment.setProperPairFlag(properPair);
                secondEndAlignment.setProperPairFlag(properPair);
                int insertSize = SamPairUtil.computeInsertSize(firstEndAlignment, secondEndAlignment);
                firstEndAlignment.setInferredInsertSize(insertSize);
                secondEndAlignment.setInferredInsertSize(-insertSize);
            }

            writer.addAlignment(firstEndAlignment);
            if (secondEndAlignment != null) {
                writer.addAlignment(secondEndAlignment);
            }
            ++numAlignmentsOrPairsWritten;
            if (MAX_ALIGNMENTS != null && numAlignmentsOrPairsWritten >= MAX_ALIGNMENTS) {
                break;
            }
            if (numAlignmentsOrPairsWritten % 500000 == 0) {
                log.info("Loaded " + numAlignmentsOrPairsWritten + " reads");
            }
        }
        if (MAX_ALIGNMENTS == null) {
            noMoreAssert(firstEndIterator);
            if (paired) {
                noMoreAssert(secondEndIterator);
            }
        }
        log.info("Done loading " + numAlignmentsOrPairsWritten + " reads");
    }

    /**
     * Write into the samRecord the mate info from the mate gerald alignment
     */
    private void setMateInfo(final SAMRecord samRecord, final GeraldParser.GeraldAlignment mateGeraldAlignment) {
        final boolean isMapped = mateGeraldAlignment.getPrimaryChrom() != null;
        if (isMapped) {
            samRecord.setMateReferenceName(mateGeraldAlignment.getPrimaryChrom());
            samRecord.setMateAlignmentStart((int)mateGeraldAlignment.getPrimaryStart());
            samRecord.setMateNegativeStrandFlag(isNegativeStrand(mateGeraldAlignment));
        } else {
            samRecord.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            samRecord.setMateUnmappedFlag(true);
        }
    }

    private boolean isNegativeStrand(final GeraldParser.GeraldAlignment alignment) {
        final String orientation = alignment.getOrientation();
        if (orientation.equals("F")) {
            return false;
        } else if (orientation.equals("R")) {
            return true;
        } else {
            throw new RuntimeException("Strange orientation in eland_extended file");
        }
    }

    private SAMRecord createSAMRecordFromGerald(final GeraldParser.GeraldAlignment alignment) {
        final SAMRecord samRecord = new SAMRecord();
        // Consider an alignment with a negative start (i.e. that hangs off the beginning of the contig)
        // to be unmapped.
        final boolean isMapped = alignment.getPrimaryChrom() != null && alignment.getPrimaryStart() >= 0;

        String readName = alignment.getReadName();
        if (readName.endsWith("/1") || readName.endsWith("/2")) {
            readName = readName.substring(0, readName.length() - 2);
        }
        samRecord.setReadName(readName);

        // Set all the flags
        samRecord.setReadPairedFlag(paired);
        samRecord.setReadUmappedFlag(!isMapped);
        if (isMapped) {
            samRecord.setReadNegativeStrandFlag(isNegativeStrand(alignment));
        }
        // For now we are only taking the primary alignment
        samRecord.setNotPrimaryAlignmentFlag(false);
        String readBases = alignment.getReadBases();
        if (samRecord.getReadNegativeStrandFlag()) {
            readBases = SequenceUtil.reverseComplement(readBases);
        }
        samRecord.setReadString(readBases);
        final byte[] phredQualities = alignment.getPhredQualities();
        if (isMapped && samRecord.getReadNegativeStrandFlag()) {
            ArrayUtil.reverseArray(phredQualities);
        }
        samRecord.setBaseQualities(phredQualities);
        if (isMapped) {
            /*
        if ("23".equals(geraldReferenceName)) {
                        geraldReferenceName = "X";
                    } else if ("24".equals(geraldReferenceName)) {
                        geraldReferenceName = "Y";
                    }
                    return REFERENCE_PREFIX + geraldReferenceName;
            */
            samRecord.setReferenceName(alignment.getPrimaryChrom());
            samRecord.setAlignmentStart((int)alignment.getPrimaryStart());
            samRecord.setMappingQuality(SAMRecord.UNKNOWN_MAPPING_QUALITY);
            // CIGAR is trivial because there are no indels or clipping in Gerald
            final String cigar = Integer.toString(alignment.getReadLength()) + "M";
            samRecord.setCigarString(cigar);
            // We've decided not to bother with this, and just load the reference
            // if we want to determine mismatches.
            // samRecord.setAttribute("MD", alignment.getMismatchString());
        } else {
            samRecord.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            samRecord.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            samRecord.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
            samRecord.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
        }

        if (SAMPLE != null) {
            // There is a read group (id = READ_GROUP_ID)
            samRecord.setAttribute("RG", READ_GROUP_ID);
        }

        samRecord.setAttribute("PG", PROGRAM_RECORD_ID);
        return samRecord;
    }

    private void hasNextAssert(final Iterator iterator) {
        if (!iterator.hasNext()) {
            throw new RuntimeException("gerald output file ends unexpectedly.");

        }
    }

    private void noMoreAssert(final Iterator iterator) {
        if (iterator.hasNext()) {
            throw new RuntimeException("gerald output file has more lines than expected.");
        }
    }

}
