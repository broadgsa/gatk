package org.broadinstitute.sting.playground.piecemealannotator;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.File;
import java.util.Date;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SplitAlignedSAMByTile extends CommandLineProgram {
    public static SplitAlignedSAMByTile instance = null;

    @Argument(fullName="sam_in", shortName="SI", doc="SAM file to split") public File SAM_IN;
    @Argument(fullName="sam_out_prefix", shortName="SOP", doc="Output prefix for split SAMs") public String SAM_OUT_PREFIX;

    public static void main(String[] argv) {
        instance = new SplitAlignedSAMByTile();
        start(instance, argv);
    }

    protected int execute() {
        SAMFileReader sreader = new SAMFileReader(SAM_IN);

        SAMFileHeader sheader = sreader.getFileHeader();
        sheader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        Vector<SAMFileWriter> swriters = new Vector<SAMFileWriter>();
        swriters.setSize(101);
        for (int tile = 1; tile <= 100; tile++) {
            File tileFile = new File(String.format("%s.%d.sam", SAM_OUT_PREFIX, tile));
            swriters.add(tile, new SAMFileWriterFactory().makeSAMOrBAMWriter(sheader, true, tileFile));
        }

        int reads = 0;
        for (SAMRecord sr : sreader) {
            Pattern p = Pattern.compile(":\\d:(\\d+):\\d+:\\d+#");
            Matcher m = p.matcher(sr.getReadName());

            if (m.find()) {
                int tile = Integer.valueOf(m.group(1));

                swriters.get(tile).addAlignment(sr);

                if (reads % 1000000 == 0) {
                    System.out.printf("%s: processed %d reads ... \n", (new Date()).toString(), reads);
                }
                reads++;
            }
        }

        sreader.close();
        for (SAMFileWriter sfw : swriters) {
            if (sfw != null) {
                sfw.close();
            }
        }

        return 0;
    }
}
