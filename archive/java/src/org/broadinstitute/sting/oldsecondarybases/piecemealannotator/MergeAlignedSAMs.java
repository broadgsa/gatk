package org.broadinstitute.sting.playground.piecemealannotator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.File;

public class MergeAlignedSAMs extends CommandLineProgram {
    public static MergeAlignedSAMs instance = null;

    @Argument(fullName="sam_tile_prefix", shortName="STP", doc="SAM tile prefix") public String SAM_TILE_PREFIX;
    @Argument(fullName="sam_tile_suffix", shortName="STS", doc="SAM tile suffix") public String SAM_TILE_SUFFIX;
    @Argument(fullName="sam_out", shortName="SO", doc="SAM output file") public File SAM_OUT;

    public static void main(String[] argv) {
        instance = new MergeAlignedSAMs();
        start(instance, argv);
    }

    protected int execute() {
        SAMFileWriter swriter = null;

        for (int tile = 1; tile <= 100; tile++) {
            File tileFile = new File(String.format("%s.%d.%s", SAM_TILE_PREFIX, tile, SAM_TILE_SUFFIX));

            SAMFileReader sreader = new SAMFileReader(tileFile);

            if (swriter == null) {
                swriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(sreader.getFileHeader(), true, SAM_OUT);
            }

            System.out.println("Processing " + tileFile.getName() + " ...");

            for (SAMRecord sr : sreader) {
                swriter.addAlignment(sr);
            }

            sreader.close();
        }

        if (swriter != null) {
            swriter.close();
        }

        return 0;
    }
}
