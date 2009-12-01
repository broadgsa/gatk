package org.broadinstitute.sting.playground.piecemealannotator;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.File;
import java.util.HashMap;

public class MergeAlignedAndSecondarySAMs extends CommandLineProgram {
    public static MergeAlignedAndSecondarySAMs instance = null;

    @Argument(fullName="unaligned_sam", shortName="US", doc="unaligned SAM with secondary bases") public File USAM;
    @Argument(fullName="aligned_sam", shortName="AS", doc="aligned SAM without secondary bases") public File ASAM;
    @Argument(fullName="sam_out", shortName="SO", doc="output SAM file") public File OSAM;

    public static void main(String[] argv) {
        instance = new MergeAlignedAndSecondarySAMs();
        start(instance, argv);
    }

    protected int execute() {
        // Hash unaligned file
        HashMap<String, SAMRecord> usamhash = new HashMap<String, SAMRecord>(500000);

        SAMFileReader usamr = new SAMFileReader(USAM);
        for (SAMRecord usr : usamr) {
            String key = String.format("%s:%b", usr.getReadName(), usr.getFirstOfPairFlag());

            usamhash.put(key, usr);
        }
        usamr.close();

        // Annotate aligned file
        SAMFileReader asamr = new SAMFileReader(ASAM);

        SAMFileHeader asamh = asamr.getFileHeader();
        asamh.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter osamw = new SAMFileWriterFactory().makeSAMOrBAMWriter(asamh, true, OSAM);

        for (SAMRecord asr : asamr) {
            String key = String.format("%s:%b", asr.getReadName(), asr.getFirstOfPairFlag());

            if (usamhash.containsKey(key)) {
                String abases = asr.getReadString();
                String ubases1 = usamhash.get(key).getReadString();
                String ubases2 = (String) usamhash.get(key).getAttribute("SB");

                if (asr.getReadNegativeStrandFlag()) {
                    ubases1 = BaseUtils.simpleReverseComplement(ubases1);
                    ubases2 = BaseUtils.simpleReverseComplement(ubases2);
                }

                byte[] sqbases = new byte[abases.length()];
                for (int cycle = 0; cycle < abases.length(); cycle++) {
                    char sqbase = (abases.charAt(cycle) == ubases1.charAt(cycle)) ? ubases2.charAt(cycle) : ubases1.charAt(cycle);
                    int sqBaseIndex = BaseUtils.simpleBaseToBaseIndex(sqbase);
                    sqbases[cycle] = QualityUtils.baseAndProbDiffToCompressedQuality(sqBaseIndex, 0.0);
                }

                asr.setAttribute("SQ", sqbases);

                osamw.addAlignment(asr);
            }
        }

        osamw.close();
        asamr.close();

        return 0;
    }
}
