package org.broadinstitute.sting.playground.piecemealannotator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.File;
import java.util.Date;
import java.util.HashMap;

public class MergePieces extends CommandLineProgram {
    public static MergePieces instance = null;

    @Argument(fullName="unaligned_sam_prefix", shortName="USP", doc="Prefix for unaligned SAM files") public String USAM_PREFIX;
    @Argument(fullName="unaligned_sam_suffix", shortName="USS", doc="Suffix for unaligned SAM files") public String USAM_SUFFIX;
    @Argument(fullName="tile_start", shortName="TS", doc="Starting tile number") public int TILE_START;
    @Argument(fullName="tile_end", shortName="TE", doc="Ending tile number (inclusive)") public int TILE_END;
    @Argument(fullName="aligned_sam_in", shortName="ASI", doc="Aligned SAM file") public File ALIGNED_SAM;
    @Argument(fullName="annotated_sam_out", shortName="ASO", doc="Annotated SAM output file") public File ANNOTATED_SAM;

    public static void main(String[] argv) {
        instance = new MergePieces();
        start(instance, argv);
    }

    protected int execute() {
        HashMap<String, byte[]> samHash = new HashMap<String, byte[]>(30000000);

        System.out.println("Hashing records...");
        for (int tile = TILE_START; tile <= TILE_END; tile++) {
            System.out.printf("  %s: Hashing SQ tags from tile %d...\n", (new Date()).toString(), tile);
            
            File uTileFile = new File(String.format("%s.%d.%s", USAM_PREFIX, tile, USAM_SUFFIX));

            SAMFileReader uTileReader = new SAMFileReader(uTileFile);

            for (SAMRecord sr : uTileReader) {
                String key = String.format("%s:%b", sr.getReadName(), sr.getReadPairedFlag() && sr.getFirstOfPairFlag());

                //System.out.printf("name=%s ispaired=%b negstrand=%b firstend=%b secondend=%b\n", sr.getReadName(), sr.getReadPairedFlag(), sr.getReadNegativeStrandFlag(), sr.getFirstOfPairFlag(), sr.getSecondOfPairFlag());

                samHash.put(key, (byte[]) sr.getAttribute("SQ"));
            }

            uTileReader.close();
        }

        SAMFileReader aReader = new SAMFileReader(ALIGNED_SAM);
        SAMFileWriter aWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(aReader.getFileHeader(), true, ANNOTATED_SAM);

        System.out.println("Annotating reads...");
        int annotatedReads = 0;
        for (SAMRecord sr : aReader) {
            String key = String.format("%s:%b", sr.getReadName(), sr.getReadPairedFlag() && sr.getFirstOfPairFlag());

            if (samHash.containsKey(key)) {
                byte[] sqtag = samHash.get(key);

                if (sr.getReadNegativeStrandFlag()) {
                    sqtag = QualityUtils.reverseComplementCompressedQualityArray(sqtag);
                }

                sr.setAttribute("SQ", sqtag);

                annotatedReads++;
                if (annotatedReads % 100000 == 0) {
                    System.out.printf("  %s: Annotated %d reads...\n", (new Date()).toString(), annotatedReads);
                }

                //aWriter.addAlignment(sr);
            }
            
            aWriter.addAlignment(sr);
        }

        aReader.close();
        aWriter.close();

        return 0;
    }
}
