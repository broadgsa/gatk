package org.broadinstitute.sting.gatk.refdata.indexer;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.utils.Utils;

import java.io.File;
import java.io.FileOutputStream;
import java.util.Map;

/**
 * a utility class that can create an index, written to a target location.  This is useful when you're unable to write to the directory
 * in which an index is located, or if you'd like to pre-index files to save time.
 */
public class RMDIndexer extends CommandLineProgram {
    @Argument(shortName="in", fullName="inputFile", doc="The reference meta data file to index", required = true)
    File inputFileSource = null;

    @Argument(shortName="t", fullName="type", doc="The reference meta data file format (e.g. vcf, bed)", required = true)
    String inputFileType = null;

    @Input(fullName = "referenceSequence", shortName = "R", doc = "The reference to use when indexing; this sequence will be set in the index", required = true)
    public File referenceFile = null;

    @Input(shortName = "i", fullName = "indexFile", doc = "Where to write the index to (as a file), if not supplied we write to <inputFile>.idx", required = false)
    public File indexFile = null;

    @Argument(shortName = "ba", fullName = "balanceApproach", doc="the index balancing approach to take", required=false)
    IndexFactory.IndexBalanceApproach approach = IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME;

    private static Logger logger = Logger.getLogger(RMDIndexer.class);

    @Override
    protected int execute() throws Exception {
        logger.info(String.format("attempting to index file:   %s",inputFileSource));
        logger.info(String.format("using reference:            %s",referenceFile.getAbsolutePath()));
        logger.info(String.format("using type:                 %s",inputFileType));
        logger.info(String.format("writing to location:        %s",indexFile.getAbsolutePath()));

        // check parameters
        // ---------------------------------------------------------------------------------

        // check the input parameters
        if (referenceFile == null || !referenceFile.canRead())
            throw new IllegalArgumentException("We can't read the reference file: "
                + referenceFile + ", check that it exists, and that you have permissions to read it");

        // create a reference file reader
        IndexedFastaSequenceFile seq = new IndexedFastaSequenceFile(referenceFile);

        // check that we can create the output file
        if (indexFile == null || indexFile.exists())
            throw new IllegalArgumentException("We can't write to the index file location: "
                + indexFile + ", the index exists");


        // try to index the file
        // ---------------------------------------------------------------------------------

        // get a track builder
        RMDTrackBuilder builder = new RMDTrackBuilder();

        // find the types available to the track builders
        Map<String,Class> typeMapping = builder.getAvailableTrackNamesAndTypes();

        // check that the type is valid
        if (!typeMapping.containsKey(inputFileType))
            throw new IllegalArgumentException("The type specified " + inputFileType + " is not a valid type.  Valid type list: " + Utils.join(",",typeMapping.keySet()));

        // create the codec
        FeatureCodec codec = builder.createByType(typeMapping.get(inputFileType));

        // get some timing info
        long currentTime = System.currentTimeMillis();

        Index index = IndexFactory.createIndex(inputFileSource, codec, approach);

        // create the output stream, and write the index
        LittleEndianOutputStream stream = new LittleEndianOutputStream(new FileOutputStream(indexFile));
        index.write(stream);
        stream.close();

        // report and exit
        logger.info("Successfully wrote the index to location: " + indexFile + " in " + ((System.currentTimeMillis() - currentTime)/1000) + " seconds");
        return 0;  // return successfully
    }


    /**
     * the generic call execute main
     * @param argv the arguments from the command line
     */
    public static void main(String[] argv) {
        try {
            RMDIndexer instance = new RMDIndexer();
            start(instance, argv);
            System.exit(CommandLineProgram.result);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }
}
