/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.aligner.maq;

import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.util.StringSortingCollectionFactory;
import edu.mit.broad.picard.util.Log;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.sam.util.SortingCollection;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.*;

import java.io.File;
import java.io.BufferedInputStream;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.nio.ByteBuffer;

/**
 * Class to write a BAM file that includes the results from a Maq .map file along with the unaligned
 * reads from the original BAM file.
 *
 * Information on the meaning of the elements of the map file is drawn from the Maq documentation
 * on this page: http://maq.sourceforge.net/maqmap_format.shtml
 */
public class MaqMapMerger {

    private final File mapFile;
    private final File sourceBamFile;
    private final File targetBamFile;
    private final boolean pairedReads;
    private final Log log = Log.getInstance(MaqMapMerger.class);
    private String commandLine = null;
    private List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>();


    /**
     * Constructor
     *
     * @param mapFile           The Maq map file to parse
     * @param sourceBamFile     The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map
     * @param targetBamFile     The file to which to write the merged
     */
    public MaqMapMerger(File mapFile, File sourceBamFile, File targetBamFile, boolean pairedReads) {
        IoUtil.assertFileIsReadable(mapFile);
        IoUtil.assertFileIsReadable(sourceBamFile);
        IoUtil.assertFileIsWritable(targetBamFile);
        this.mapFile = mapFile;
        this.sourceBamFile = sourceBamFile;
        this.targetBamFile = targetBamFile;
        this.pairedReads = pairedReads;
    }

    /**
     * Merges the alignment from the map file with the remaining records from the source BAM file.
     */
    public void mergeAlignment() {
        log.info("Processing map file: " + mapFile.getAbsolutePath());
        // Write the header
        MapFileIterator it = new MapFileIterator(getCommandLine(), this.pairedReads, false, this.mapFile);
        SAMFileHeader header = it.getHeader();
        SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, targetBamFile);

        // Write the alignments
        SortingCollection<String> readNames = writeAlignments(it, writer);

        // We're done with the map file, so close it
        it.close();
        writeUnalignedReads(writer, readNames.iterator());

        // Now close the writer
        writer.close();
    }

    
    private void writeUnalignedReads(SAMFileWriter writer, CloseableIterator<String> nameIterator) {

        int skipCount = 0;
        SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(this.sourceBamFile));
        CloseableIterator<SAMRecord> bamRecords = reader.iterator();

        String readName = nameIterator.hasNext() ? nameIterator.next() : null;
        while(bamRecords.hasNext()) {
            SAMRecord rec = bamRecords.next();
            if (rec.getReadName().equals(readName)) {
                // skip it and pull the next name off the name iterator
                readName = nameIterator.hasNext() ? nameIterator.next() : null;
                skipCount++;
            }
            else {
                writer.addAlignment(rec);
            }
        }
System.out.println("Skipped " + skipCount + " already-aligned records.");
        bamRecords.close();
        nameIterator.close();
    }

    private SortingCollection<String> writeAlignments(MapFileIterator iterator, SAMFileWriter writer) {

int wrote = 0;
        SortingCollection<String> readNames = StringSortingCollectionFactory.newCollection();
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            readNames.add(record.getReadName());
            writer.addAlignment(record);
wrote++;
        }
System.out.println("Wrote " + wrote + " alignment records.");
        return readNames;
    }

    public void setCommandLine(String commandLine) { this.commandLine = commandLine; }
    public String getCommandLine() { return this.commandLine; }
}
