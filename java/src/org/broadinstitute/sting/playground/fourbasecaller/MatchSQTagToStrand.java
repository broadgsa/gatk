package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.QualityUtils;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

public class MatchSQTagToStrand extends CommandLineProgram {
    public static MatchSQTagToStrand Instance = null;

    public File SAM_IN;
    public File SAM_OUT;

    public static void main(String[] argv) {
        Instance = new MatchSQTagToStrand();
        start(Instance, argv);
    }

    protected void setupArgs() {
        m_parser.addRequiredArg("sam_in",  "I", "Input SAM file",  "SAM_IN");
        m_parser.addRequiredArg("sam_out", "O", "Output SAM file", "SAM_OUT");
    }

    protected int execute() {
        SAMFileReader sf = new SAMFileReader(SAM_IN);
        sf.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        SAMFileWriter sw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sf.getFileHeader(), true, SAM_OUT);

        for (SAMRecord sr : sf) {
            if (sr.getReadNegativeStrandFlag()) {
                byte[] sq = (byte[]) sr.getAttribute("SQ");
                sq = QualityUtils.reverseComplementCompressedQualityArray(sq);

                sr.setAttribute("SQ", sq);
            }

            sw.addAlignment(sr);
        }

        sw.close();

        return 0;
    }
}
