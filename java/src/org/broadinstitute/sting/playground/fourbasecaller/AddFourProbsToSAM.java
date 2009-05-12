package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;

import java.io.File;
import java.util.HashMap;

import net.sf.samtools.*;

public class AddFourProbsToSAM extends CommandLineProgram {
    public static AddFourProbsToSAM Instance = null;

    public File UNALIGNED_SAM;
    public File ALIGNED_SAM;
    public File FINAL_SAM;
    public int END;
    public Boolean DEBUG = false;

    public static void main(String[] argv) {
        Instance = new AddFourProbsToSAM();
        start(Instance, argv);
    }

    protected void setupArgs() {
        //m_parser.addRequiredArg("unaligned_sam", "U", "Unaligned SAM file", "UNALIGNED_SAM");
        //m_parser.addRequiredArg("aligned_sam",   "A", "Aligned SAM file",   "ALIGNED_SAM");
        //m_parser.addRequiredArg("final_sam",     "F", "Final SAM file",     "FINAL_SAM");
        //m_parser.addRequiredArg("end",           "E", "Pair end (0 - all, 1 - first, 2 - second)", "END");
        //m_parser.addOptionalFlag("debug",        "D", "Turn on debugging output", "DEBUG");
    }

    protected int execute() {
        int processed;

        SAMFileReader alignedSf = new SAMFileReader(ALIGNED_SAM);
        alignedSf.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        // First, hash the aligned records (because there are less of them than unaligned reads)
        System.err.println("Hashing aligned records...");

        HashMap<String, SAMRecord> records = new HashMap<String, SAMRecord>(10000000);
        processed = 0;
        for (SAMRecord alignedSr : alignedSf) {
            if (END == 0 || (END == 1 && alignedSr.getSecondOfPairFlag() == false) || (END == 2 && alignedSr.getSecondOfPairFlag() == true)) {
                if (!alignedSr.getReadUnmappedFlag()) {
                    records.put(alignedSr.getReadName(), alignedSr);

                    if (processed % 100000 == 0) { System.err.print("\tProcessed " + processed + " records.\r"); }
                    processed++;
                }
            }
        }

        // Now, iterate over the unaligned SAM file and stick the four-base probs in.
        System.err.println("\nInterating over unaligned records...");

        SAMFileReader unalignedSf = new SAMFileReader(UNALIGNED_SAM);
        unalignedSf.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        SAMFileHeader swhead = alignedSf.getFileHeader();
        swhead.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter sw = new SAMFileWriterFactory().makeSAMOrBAMWriter(swhead, true, FINAL_SAM);

        processed = 0;
        for (SAMRecord unalignedSr : unalignedSf) {
            if (records.containsKey(unalignedSr.getReadName())) {
                SAMRecord alignedSr = records.get(unalignedSr.getReadName());

                byte[] sq = (byte[]) unalignedSr.getAttribute("SQ");
                if (alignedSr.getReadNegativeStrandFlag()) {
                    sq = QualityUtils.reverseComplementCompressedQualityArray(sq);
                }
                
                alignedSr.setAttribute("SQ", sq);
                alignedSr.setAttribute("KB", unalignedSr.getReadBases());
                alignedSr.setAttribute("KQ", unalignedSr.getBaseQualities());
                sw.addAlignment(alignedSr);

                if (DEBUG) {
                    System.out.println(alignedSr.format());
                }

                if (processed % 100000 == 0) { System.err.print("\tProcessed " + processed + " records.\r"); }
                processed++;
            }
        }

        sw.close();
        alignedSf.close();
        unalignedSf.close();

        return 0;
    }
}
