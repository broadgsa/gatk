package org.broadinstitute.sting.gatk.io.storage;

import org.broad.tribble.readers.LineReader;
import org.broad.tribble.source.BasicFeatureSource;
import org.broad.tribble.vcf.*;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub;

import java.io.*;

import net.sf.samtools.util.BlockCompressedOutputStream;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

/**
 * Provides temporary and permanent storage for genotypes in VCF format.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFWriterStorage implements Storage<VCFWriterStorage>, VCFWriter {
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
            writer = new StandardVCFWriter(stream);
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
        return new StandardVCFWriter(file, this.stream, indexOnTheFly && !stub.isCompressed());
    }


    /**
     * Constructs an object which will redirect into a different file.
     * @param stub Stub to use when synthesizing file / header info.
     * @param tempFile File into which to direct the output data.
     */
    public VCFWriterStorage(VCFWriterStub stub, File tempFile) {
        this.file = tempFile;
        this.writer = vcfWriterToFile(stub, file, false);
        writer.writeHeader(stub.getVCFHeader());
    }

    public void add(VariantContext vc, byte ref) {
        writer.add(vc, ref);
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
        writer.close();
    }

    public void mergeInto(VCFWriterStorage target) {
        try {
            //System.out.printf("merging %s%n", file);
            BasicFeatureSource<VariantContext> source = BasicFeatureSource.getFeatureSource(file.getAbsolutePath(), new VCFCodec(), false);
            
            for ( VariantContext vc : source.iterator() ) {
                target.writer.add(vc, vc.getReferenceBaseForIndel());
            }

            source.close();
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(file, "Error reading file in VCFWriterStorage: ", e);
        }
    }
}
