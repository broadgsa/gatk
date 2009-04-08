package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.QualityUtils;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

public class RepairBadlyCombinedSamFile extends CommandLineProgram {
    public static RepairBadlyCombinedSamFile Instance = null;

    public File ASAM;
    public File OSAM;

    public static void main(String[] argv) {
        Instance = new RepairBadlyCombinedSamFile();
        start(Instance, argv);
    }

    protected void setupArgs() {
        m_parser.addRequiredArg("aligned_sam",   "A", "Second input SAM file", "ASAM");
        m_parser.addRequiredArg("sam_out",       "O", "Output SAM file",       "OSAM");
    }

    protected int execute() {
        SAMFileReader asamr = new SAMFileReader(ASAM);
        asamr.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        SAMFileWriter osamw = new SAMFileWriterFactory().makeSAMOrBAMWriter(asamr.getFileHeader(), false, OSAM);


        for (SAMRecord samRecord : asamr) {
            if (samRecord.getReadNegativeStrandFlag()) {
                byte[] bases = samRecord.getReadBases();
                byte[] quals = samRecord.getBaseQualities();
                byte[] sq    = (byte[]) samRecord.getAttribute("SQ");

                byte[] newbases = new byte[bases.length];
                byte[] newquals = new byte[bases.length];
                byte[] newsq    = new byte[bases.length];

                for (int i = 0; i < bases.length; i++) {
                    byte newbase;
                    switch (bases[bases.length - i - 1]) {
                        case 'A': newbase = 'T'; break;
                        case 'C': newbase = 'G'; break;
                        case 'G': newbase = 'C'; break;
                        case 'T': newbase = 'A'; break;
                        default: newbase = '.';
                    }
                    newbases[i] = newbase;
                    newquals[i] = quals[bases.length - i - 1];
                    newsq[i]    = QualityUtils.complementCompressedQuality(newquals[bases.length - i - 1]);
                }

                samRecord.setReadBases(newbases);
                samRecord.setBaseQualities(newquals);
                samRecord.setAttribute("SQ", newsq);
            }
            
            osamw.addAlignment(samRecord);
        }

        return 0;
    }
}
