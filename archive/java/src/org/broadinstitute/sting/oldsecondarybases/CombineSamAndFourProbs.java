package org.broadinstitute.sting.secondarybase;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class CombineSamAndFourProbs extends CommandLineProgram {
    public static CombineSamAndFourProbs Instance = null;

    @Argument(fullName="sam",      shortName="S", doc="Input SAM file")
    public File SAM;
    @Argument(fullName="fourprob", shortName="F", doc="Input text file := read_name sq_field")
    public File FOURPROBS;
    @Argument(fullName="sam_out",  shortName="O", doc="Output SAM file")
    public File SAM_OUT;

    public static void main(String[] argv) {
        Instance = new CombineSamAndFourProbs();
        start(Instance, argv);
    }

    protected int execute() {
        BufferedReader fpreader = null;

        try {
            fpreader = new BufferedReader(new FileReader(FOURPROBS));

            HashMap fourprobMap = new HashMap(27000000);

            String fourprobLine;
            int processed = 0;
            while ((fourprobLine = fpreader.readLine()) != null) {
                String[] fourprobPieces = fourprobLine.split("\\s+");
                String[] sqfield = fourprobPieces[1].split(":");
                byte[] sq = hexStringToBytes(sqfield[2]);

                fourprobMap.put(fourprobPieces[0], sq);

                if (processed % 1000000 == 0) {
                    System.out.println("Processed " + processed);
                }
                processed++;
            }

            fpreader.close();

            SAMFileReader sf = new SAMFileReader(SAM);
            sf.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

            SAMFileWriter sw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sf.getFileHeader(), true, SAM_OUT);

            for (SAMRecord sr : sf) {
                String readname = sr.getReadName();
                byte[] sq = (byte[]) fourprobMap.get(readname);

                if (sq != null) {
                    if (sr.getReadNegativeStrandFlag()) {
                        sq = QualityUtils.reverseComplementCompressedQualityArray(sq);
                    }

                    sr.setAttribute("SQ", sq);
                    sw.addAlignment(sr);
                }
            }

            sf.close();
            sw.close();
        } catch (IOException e) {
            System.err.println("There was an error.");
            System.exit(1);
        }

        return 0;
    }

    static String bytesToHexString(final byte[] data) {
        final char[] chars = new char[2 * data.length];
        for (int i = 0; i < data.length; i++) {
            final byte b = data[i];
            chars[2*i] = toHexDigit((b >> 4) & 0xF);
            chars[2*i+1] = toHexDigit(b & 0xF);
        }
        return new String(chars);
    }
    
    static byte[] hexStringToBytes(final String s)  throws NumberFormatException {
        if (s.length() % 2 != 0) {
            throw new NumberFormatException("Hex representation of byte string does not have even number of hex chars: " + s);
        }
        final byte[] ret = new byte[s.length() / 2];
        for (int i = 0; i < ret.length; ++i) {
            ret[i] = (byte) ((fromHexDigit(s.charAt(i * 2)) << 4) | fromHexDigit(s.charAt(i * 2 + 1)));
        }
        return ret;
    }

    private static char toHexDigit(final int value) {
        return (char) ((value < 10) ? ('0' + value) : ('A' + value - 10));
    }
    
    private static int fromHexDigit(final char c) throws NumberFormatException {
        final int ret = Character.digit(c, 16);
        if (ret == -1) {
            throw new NumberFormatException("Not a valid hex digit: " + c);
        }
        return ret;
    }
}

