package org.broadinstitute.sting.gatk.io.storage;

import org.broad.tribble.vcf.StandardVCFWriter;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub;

import java.io.*;
import java.util.Set;

import net.sf.samtools.util.BlockCompressedOutputStream;

/**
 * Provides temporary and permanent storage for genotypes in VCF format.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFWriterStorage implements Storage<VCFWriterStorage>, VCFWriter {
    protected final File file;
    protected final OutputStream stream;
    protected final VCFWriter writer;

    /**
     * Constructs an object which will write directly into the output file provided by the stub.
     * Intentionally delaying the writing of the header -- this should be filled in by the walker.
     * @param stub Stub to use when constructing the output file.
     */
    public VCFWriterStorage( VCFWriterStub stub )  {

        if ( stub.getFile() != null ) {
            file = stub.getFile();
            try {
                if ( stub.isCompressed() )
                    stream = new BlockCompressedOutputStream(file);
                else
                    stream = new PrintStream(file);
            }
            catch(IOException ex) {
                throw new StingException("Unable to open target output stream", ex);
            }
        }
        else if ( stub.getOutputStream() != null ) {
            this.file = null;
            this.stream = stub.getOutputStream();
        }
        else
            throw new StingException("Unable to create target to which to write; storage was provided with neither a file nor a stream.");

        writer = new StandardVCFWriter(stream);
    }

    /**
     * Constructs an object which will redirect into a different file.
     * @param stub Stub to use when synthesizing file / header info.
     * @param file File into which to direct the output data.
     */
    public VCFWriterStorage(VCFWriterStub stub, File file) {
        this.file = file;
        try {
            this.stream = new PrintStream(file);
        }
        catch(IOException ex) {
            throw new StingException("Unable to open target output stream",ex);
        }
        writer = new StandardVCFWriter(this.stream);
        Set<String> samples = SampleUtils.getSAMFileSamples(stub.getSAMFileHeader());
        writer.writeHeader(new VCFHeader(null, samples));
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

    /**
     * Merges the stream backing up this temporary storage into the target.
     * @param target Target stream for the temporary storage.  May not be null.
     */
    public void mergeInto(VCFWriterStorage target) {
        PrintStream formattingTarget = new PrintStream(target.stream);
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line = reader.readLine();
            while ( line != null ) {
                if (!VCFHeaderLine.isHeaderLine(line))
                    formattingTarget.printf("%s%n",line);
                line = reader.readLine();
            }

            reader.close();
        } catch (IOException e) {
            throw new StingException("Error reading file " + file + " in GATKVCFWriter: ", e);
        }
    }    
}
