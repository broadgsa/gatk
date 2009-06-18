package org.broadinstitute.sting.playground.piecemealannotator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.File;

public class VerifySAM extends CommandLineProgram {
    public static VerifySAM instance = null;

    @Argument(fullName="sam_in", shortName="SI", doc="SAM file to check") public File SAM_IN;

    public static void main(String[] argv) {
        instance = new VerifySAM();
        start(instance, argv);
    }

    protected int execute() {
        SAMFileReader sfr = new SAMFileReader(SAM_IN);

        for ( SAMRecord sr : sfr ) {
            if (sr.getAttribute("SQ") == null) {
                throw new StingException(String.format("SAMRecord is missing an SQ tag\n%s", sr.format()));
            }
        }

        sfr.close();

        return 0;
    }
}
