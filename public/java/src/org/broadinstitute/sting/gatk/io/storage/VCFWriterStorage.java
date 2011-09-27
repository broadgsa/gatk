package org.broadinstitute.sting.gatk.io.storage;

import net.sf.samtools.util.BlockCompressedOutputStream;
import org.apache.log4j.Logger;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub;
import org.broadinstitute.sting.utils.codecs.vcf.StandardVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * Provides temporary and permanent storage for genotypes in VCF format.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFWriterStorage implements Storage<VCFWriterStorage>, VCFWriter {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(VCFWriterStorage.class);

    protected final File file;
    protected OutputStream stream;
    protected final VCFWriter writer;

    /**
     * Constructs an object which will write directly into the output file provided by the stub.
     * Intentionally delaying the writing of the header -- this should be filled in by the walker.
     * @param stub Stub to use when constructing the output file.
     */
    public VCFWriterStorage( VCFWriterStub stub )  {
        if ( stub.getFile() != null ) {
            this.file = stub.getFile();
            writer = vcfWriterToFile(stub,stub.getFile(),true);
        }
        else if ( stub.getOutputStream() != null ) {
            this.file = null;
            this.stream = stub.getOutputStream();
            writer = new StandardVCFWriter(stream, stub.getMasterSequenceDictionary(), stub.doNotWriteGenotypes());
        }
        else
            throw new ReviewedStingException("Unable to create target to which to write; storage was provided with neither a file nor a stream.");
    }

    /**
     * common initialization routine for multiple constructors
     * @param stub Stub to use when constructing the output file.
     * @param file Target file into which to write VCF records.
     * @param indexOnTheFly true to index the file on the fly.  NOTE: will be forced to false for compressed files.
     * @return A VCF writer for use with this class
     */
    private StandardVCFWriter vcfWriterToFile(VCFWriterStub stub, File file, boolean indexOnTheFly) {
        try {
            if ( stub.isCompressed() )
                stream = new BlockCompressedOutputStream(file);
            else
                stream = new PrintStream(file);
        }
        catch(IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, "Unable to open target output stream", ex);
        }

        // The GATK/Tribble can't currently index block-compressed files on the fly.  Disable OTF indexing even if the user explicitly asked for it.
        return new StandardVCFWriter(file, this.stream, stub.getMasterSequenceDictionary(), indexOnTheFly && !stub.isCompressed(), stub.doNotWriteGenotypes());
    }


    /**
     * Constructs an object which will redirect into a different file.
     * @param stub Stub to use when synthesizing file / header info.
     * @param tempFile File into which to direct the output data.
     */
    public VCFWriterStorage(VCFWriterStub stub, File tempFile) {
        logger.debug("Creating temporary VCF file " + tempFile.getAbsolutePath() + " for VCF output.");
        this.file = tempFile;
        this.writer = vcfWriterToFile(stub, file, false);
        writer.writeHeader(stub.getVCFHeader());
    }

    public void add(VariantContext vc) {
        writer.add(vc);
    }

    /**
     * initialize this VCF header
     *
     * @param header  the header
     */
    public void writeHeader(VCFHeader header) {
        writer.writeHeader(header);
    }

    /**
     * Close the VCF storage object.
     */
    public void close() {
        if(file != null)
            logger.debug("Closing temporary file " + file.getAbsolutePath());
        writer.close();
    }

    public void mergeInto(VCFWriterStorage target) {
        try {
            String sourceFilePath = file.getAbsolutePath();
            String targetFilePath = target.file != null ? target.file.getAbsolutePath() : "/dev/stdin";
            logger.debug(String.format("Merging %s into %s",sourceFilePath,targetFilePath));
            BasicFeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(file.getAbsolutePath(), new VCFCodec(), false);
            
            for ( VariantContext vc : source.iterator() ) {
                target.writer.add(vc);
            }

            source.close();
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(file, "Error reading file in VCFWriterStorage: ", e);
        }
    }
}
