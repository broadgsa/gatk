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
    @Option(shortName="S", doc="If true, then every contig in the reads must be present in the new reference sequence.  By default, reads without a corresponding contig in the new reference are made unmapped", optional=true)
    public Boolean REQUIRE_COMPLETE_DICT_CONCORDANCE = false;

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

        // write the reads in contig order
        for ( SAMSequenceRecord contig : refDict.getSequences() ) {
            System.out.println("Writing reads on contig " + contig.getSequenceName());
            writeReads(outWriter, inReader.query(contig.getSequenceName(), 0, 0, false), newOrder);
        }
        // write the unmapped reads
        writeReads( outWriter, inReader.queryUnmapped(), newOrder);

        inReader.close();
        outWriter.close();
        return 0;
    }

    private static int newOrderIndex(SAMRecord read, int oldIndex, Map<Integer, Integer> newOrder) {
        if ( oldIndex == -1 )
            return -1; // unmapped read
        else if ( ! newOrder.containsKey(oldIndex) )
            throw new PicardException("BUG: no mapping found for read " + read.format());
        else
            return newOrder.get(oldIndex);
    }

    private static void writeReads(SAMFileWriter out, SAMRecordIterator it, Map<Integer, Integer> newOrder) {
        while ( it.hasNext() ) {
            SAMRecord read = it.next();
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
    }

    private Map<Integer, Integer> buildSequenceDictionaryMap(SAMSequenceDictionary refDict, SAMSequenceDictionary readsDict) {
        Map<Integer, Integer> newOrder = new HashMap<Integer, Integer>();

        for ( SAMSequenceRecord refRec : refDict.getSequences() ) {
            SAMSequenceRecord readsRec = readsDict.getSequence(refRec.getSequenceName());
            if ( readsRec != null ) {
                System.out.printf("Mapping %s [%d] => %s [%d]%n",
                        readsRec.getSequenceName(), readsRec.getSequenceIndex(),
                        refRec.getSequenceName(), refRec.getSequenceIndex());
                newOrder.put(readsRec.getSequenceIndex(), refRec.getSequenceIndex());
            }
        }

        for ( SAMSequenceRecord readsRec : readsDict.getSequences() ) {
            if ( ! newOrder.containsKey(readsRec.getSequenceIndex()) ) {
                if ( REQUIRE_COMPLETE_DICT_CONCORDANCE )
                    throw new PicardException("No matching contig in new reference found for " + readsRec.getSequenceName());
                else
                    newOrder.put(readsRec.getSequenceIndex(), -1);
            }
        }


        return newOrder;
    }

    private static void printDictionary(String name, SAMSequenceDictionary dict) {
        System.out.println(name);
        for ( SAMSequenceRecord contig : dict.getSequences() ) {
            System.out.println(contig.getSequenceName());
        }
    }
}


