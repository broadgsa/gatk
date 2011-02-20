package org.broadinstitute.sting.playground.tools;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.*;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * User: mdepristo
 *
 * Reorders a SAM/BAM input file according to the order of contigs in a second reference sequence
 */
public class ReorderSam extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Reorders reads in a SAM/BAM file to match the contig ordering in a provided reference file, as determined by exact name matching of contigs.  Reads mapped to contigs absent in the new reference are dropped";
    @Option(shortName="I", doc="Input file (bam or sam) to extract reads from.", optional=false)
    public File IN = null;
    @Option(shortName="O",doc="Output file (bam or sam) to write extracted reads to.", optional=false)
    public File OUT = null;
    @Option(shortName="R", doc="Reference sequence to reorder reads to match", optional=false)
    public File REFERENCE = null;
    @Option(shortName="S", doc="If true, then allows only a partial overlap of the BAM contigs with the new reference sequence contigs.  By default, this tool requires a corresponding contig in the new reference for each read contig", optional=true)
    public Boolean ALLOW_INCOMPLETE_DICT_CONCORDANCE = false;
    @Option(shortName="U", doc="If true, then permits mapping from a read contig to a new reference contig with the same name but a different length.  Highly dangerous, only use if you know what you are doing.", optional=true)
    public Boolean ALLOW_CONTIG_LENGTH_DISCORDANCE = false;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new ReorderSam().instanceMain(argv));
    }

    protected int doWork() {
        SAMFileReader inReader = new SAMFileReader(IN);

        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
        SAMSequenceDictionary refDict = reference.getSequenceDictionary();

        if ( refDict == null ) {
        	System.out.println("No reference sequence dictionary found. Aborting.");
        	inReader.close();
        	System.exit(1);
        }

        printDictionary("SAM/BAM file", inReader.getFileHeader().getSequenceDictionary());
        printDictionary("Reference", refDict);
        Map<Integer, Integer> newOrder = buildSequenceDictionaryMap(refDict, inReader.getFileHeader().getSequenceDictionary());

        // has to be after we create the newOrder
        SAMFileHeader outHeader = inReader.getFileHeader().clone();
        outHeader.setSequenceDictionary(refDict);
        SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader, true, OUT) ;

        //
        // write the reads in contig order
        //
        System.out.println("Writing reads...");
        for ( SAMSequenceRecord contig : refDict.getSequences() ) {
            SAMRecordIterator it = inReader.query(contig.getSequenceName(), 0, 0, false);
            writeReads(outWriter, it, newOrder, contig.getSequenceName());
        }
        // don't forget the unmapped reads
        writeReads( outWriter, inReader.queryUnmapped(), newOrder, "unmapped" );

        // cleanup
        inReader.close();
        outWriter.close();
        return 0;
    }

    /**
     * Low-level helper function that returns the new reference index for oldIndex according to the
     * ordering map newOrder.  Read is provided in case an error occurs, so that an informative message
     * can be made.
     *
     * @param read
     * @param oldIndex
     * @param newOrder
     * @return
     */
    private static int newOrderIndex(SAMRecord read, int oldIndex, Map<Integer, Integer> newOrder) {
        if ( oldIndex == -1 )
            return -1; // unmapped read
        else if ( ! newOrder.containsKey(oldIndex) )
            throw new PicardException("BUG: no mapping found for read " + read.format());
        else
            return newOrder.get(oldIndex);
    }

    /**
     * Helper function that writes reads from iterator it into writer out, updating each SAMRecord along the way
     * according to the newOrder mapping from dictionary index -> index.  Name is used for printing only.
     *
     * @param out
     * @param it
     * @param newOrder
     * @param name
     */
    private static void writeReads(SAMFileWriter out, SAMRecordIterator it, Map<Integer, Integer> newOrder, String name) {
        long counter = 0;
        System.out.print("  Processing " + name);
        System.out.flush();

        while ( it.hasNext() ) {
            counter++;
            SAMRecord read = it.next();
            read.setHeader(out.getFileHeader());
            //System.out.println("Writing read " + read.format());

            int oldRefIndex = read.getReferenceIndex();
            int newRefIndex = newOrderIndex(read, oldRefIndex, newOrder);
            read.setReferenceIndex(newRefIndex);

            int oldMateIndex = read.getMateReferenceIndex();
            int newMateIndex = newOrderIndex(read, oldMateIndex, newOrder);
            if ( oldMateIndex != -1 && newMateIndex == -1 ) { // becoming unmapped
                read.setMateAlignmentStart(0);
                read.setMateUnmappedFlag(true);
            }
            read.setMateReferenceIndex(newMateIndex);

            out.addAlignment(read);
        }

        it.close();
        System.out.println(" => Wrote " + counter + " reads");
    }

    /**
     * Constructs a mapping from read sequence records index -> new sequence dictionary index for use in
     * reordering the reference index and mate reference index in each read.  -1 means unmapped.
     * @param refDict
     * @param readsDict
     * @return
     */
    private Map<Integer, Integer> buildSequenceDictionaryMap(SAMSequenceDictionary refDict, SAMSequenceDictionary readsDict) {
        Map<Integer, Integer> newOrder = new HashMap<Integer, Integer>();

        System.out.println("Reordering SAM/BAM file:");
        for ( SAMSequenceRecord refRec : refDict.getSequences() ) {
            SAMSequenceRecord readsRec = readsDict.getSequence(refRec.getSequenceName());
            if ( readsRec != null ) {
                if ( refRec.getSequenceLength() != readsRec.getSequenceLength() ) {
                    String msg = String.format("Discordant contig lengths: read %s LN=%d, ref %s LN=%d",
                            readsRec.getSequenceName(), readsRec.getSequenceLength(),
                            refRec.getSequenceName(), refRec.getSequenceLength());
                    if ( ALLOW_CONTIG_LENGTH_DISCORDANCE ) {
                        System.err.println("WARN: " + msg);
                    } else {
                        throw new PicardException(msg);
                    }
                }

                System.out.printf("  Reordering read contig %s [index=%d] to => ref contig %s [index=%d]%n",
                        readsRec.getSequenceName(), readsRec.getSequenceIndex(),
                        refRec.getSequenceName(), refRec.getSequenceIndex());
                newOrder.put(readsRec.getSequenceIndex(), refRec.getSequenceIndex());
            }
        }

        for ( SAMSequenceRecord readsRec : readsDict.getSequences() ) {
            if ( ! newOrder.containsKey(readsRec.getSequenceIndex()) ) {
                if ( ALLOW_INCOMPLETE_DICT_CONCORDANCE )
                    newOrder.put(readsRec.getSequenceIndex(), -1);
                else
                    throw new PicardException("New reference sequence does not contain a matching contig for " + readsRec.getSequenceName());
            }
        }

        return newOrder;
    }

    /**
     * Helper function to print out a sequence dictionary
     * @param name
     * @param dict
     */
    private static void printDictionary(String name, SAMSequenceDictionary dict) {
        System.out.println(name);
        for ( SAMSequenceRecord contig : dict.getSequences() ) {
            System.out.printf( "  SN=%s LN=%d%n", contig.getSequenceName(), contig.getSequenceLength());
        }
    }
}


