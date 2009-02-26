/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.sam;

import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.SAMFileWriter;
import edu.mit.broad.sam.SAMFileWriterFactory;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.PicardException;

import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.math.BigInteger;

/**
 * Create a SAM/BAM file from a fasta containing reference sequence.  The output SAM file contains a header but no
 * SAMRecords, and the header contains only sequence records.
 */
public class CreateSequenceDictionary extends CommandLineProgram {

    private static final String PROGRAM_VERSION = "1.0";

    // The following attributes define the command-line arguments
    @Usage(programVersion=PROGRAM_VERSION)
    public String USAGE =
            "Usage: " + getClass().getName() + " [options]\n\n" +
                    "Read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary.\n";

    @Option(doc = "Input reference fasta or fasta.gz")
    public File REFERENCE;

    @Option(doc = "Output SAM or BAM file containing only the sequence dictionary")
    public File OUTPUT;

    @Option(doc = "Put into AS field of sequence dictionary entry if supplied", optional = true)
    public String GENOME_ASSEMBLY;

    @Option(doc = "Put into UIR field of sequence dictionary entry.  If not supplied, input reference file is used",
            optional = true)
    public String URI;

    @Option(doc = "Put into SP field of sequence dictionary entry", optional = true)
    public String SPECIES;

    private final MessageDigest md5;

    public CreateSequenceDictionary() {
        try {
            md5 = MessageDigest.getInstance("MD5");
        } catch (NoSuchAlgorithmException e) {
            throw new PicardException("MD5 algorithm not found", e);
        }
    }

    public static void main(final String[] argv) {
        System.exit(new CreateSequenceDictionary().instanceMain(argv));
    }

    /**
     * Use reference filename to create URI to go into header if URI was not passed on cmd line.
     */
    protected boolean customCommandLineValidation() {
        if (URI == null) {
            URI = "file:" + REFERENCE.getAbsolutePath();
        }
        return true;
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    protected int doWork() {
        final List<SAMSequenceRecord> sequences = makeSequenceDictionary(REFERENCE);
        final SAMFileHeader samHeader = new SAMFileHeader();
        samHeader.setSequences(sequences);
        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samHeader, false, OUTPUT);
        samWriter.close();
        return 0;
    }


    /**
     * Read all the sequences from the given reference file, and convert into SAMSequenceRecords
     * @param referenceFile fasta or fasta.gz
     * @return SAMSequenceRecords containing info from the fasta, plus from cmd-line arguments.
     */
    List<SAMSequenceRecord> makeSequenceDictionary(final File referenceFile) {
        final ReferenceSequenceFile refSeqFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile);
        ReferenceSequence refSeq;
        final List<SAMSequenceRecord> ret = new ArrayList<SAMSequenceRecord>();
        while ((refSeq = refSeqFile.nextSequence()) != null) {
            ret.add(makeSequenceRecord(refSeq));
        }
        return ret;
    }

    /**
     * Create one SAMSequenceRecord from a single fasta sequence
     */
    private SAMSequenceRecord makeSequenceRecord(final ReferenceSequence refSeq) {
        final SAMSequenceRecord ret = new SAMSequenceRecord(refSeq.getName());
        ret.setSequenceLength(refSeq.length());

        // Compute MD5 of upcased bases
        final byte[] bases = refSeq.getBases();
        for (int i = 0; i < bases.length; ++i) {
                bases[i] = (byte) (Character.toUpperCase(bases[i]) & 0xff);
            }

        ret.setAttribute(SAMSequenceRecord.MD5_TAG, md5Hash(bases));
        if (GENOME_ASSEMBLY != null) {
            ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, GENOME_ASSEMBLY);
        }
        ret.setAttribute(SAMSequenceRecord.URI_TAG, URI);
        if (SPECIES != null) {
                ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, SPECIES);
            }
        return ret;
    }

    private String md5Hash(final byte[] bytes) {
        md5.reset();
        md5.update(bytes);
        return new BigInteger(1, md5.digest()).toString(16);
    }
}
