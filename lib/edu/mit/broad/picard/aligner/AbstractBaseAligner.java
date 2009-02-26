/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.aligner;

import edu.mit.broad.picard.io.IoUtil;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.List;

/**
 * Abstract base class for use by <code>Aligner</code> implementations.  Provides a constructor and
 * accessors for common inputs and outputs.
 *
 * @author Kathleen Tibbetts
 */
public abstract class AbstractBaseAligner implements Aligner {

    private final Stringency stringency;        // The stringency of the alignment
    private final File readsBamFile;            // The BAM file containing the read data
    private final String outputPrefix;          // The directory and file name prefix for outputs
    private final String referenceFileDir;      // The directory where the reference file can be found
    private final int clipPoints[];             // The clip points to use
    private final Integer expectedInsertSize;   // Expected insert size; null for non-paired-end lanes
    private final Integer readsToAlign;         // The number of reads to align (all if null)
    private final boolean pairedReads;          // Whether this is a paired-end run
    private final int readLength;
    // Parameters specific to the Aligner implementation being used
    private final Map<String, String> customParametersMap;

    /**
     * Constructor that sets every parameter. 
     *
     * @param stringency            the stringency of the alignment
     * @param readsBamFile          the BAM file containing the reads
     * @param outputPrefix          the directory and filename prefix for output
     * @param referenceFileDir      the directory where the reference file is located
     * @param clipPoints            the clip points
     * @param expectedInsertSize    the expected insert size (null for non-PE lanes)
     * @param readsToAlign          the number of reads to align
     * @param customParametersMap   parameters specific to the Aligner implementation
     */
    public AbstractBaseAligner(Stringency stringency, File readsBamFile, String outputPrefix,
                               String referenceFileDir, int clipPoints[], Integer expectedInsertSize,
                               Integer readsToAlign, Map<String, String> customParametersMap,
                               boolean pairedReads, int readLength) {

        // First, a little validation
        if (clipPoints != null && clipPoints.length != 4) {
            throw new IllegalArgumentException("Length of clipPoints array argument must be 4.");
        }
        IoUtil.assertFileIsReadable(readsBamFile);

        this.stringency = stringency;
        this.readsBamFile = readsBamFile;
        this.outputPrefix = outputPrefix;
        this.referenceFileDir = referenceFileDir;
        this.clipPoints = clipPoints != null ? clipPoints : new int[4];
        this.expectedInsertSize = expectedInsertSize;
        this.readsToAlign = readsToAlign;
        this.customParametersMap = customParametersMap;
        this.pairedReads = pairedReads;
        this.readLength = readLength;
    }

    /**
     * Utility method for deleting a list of files, to be used by the
     * cleanup method of sub-classes
     *
     * @param files         the list of files to delete
     */
    protected final void deleteFiles(List<File> files) {
        for (File f : files) {
            f.delete();
        }
    }

    // Accessors
    protected final Stringency getStringency() { return stringency; }
    protected final File getReadsBamFile() { return readsBamFile; }
    protected final String getOutputPrefix() { return outputPrefix; }
    protected final String getReferenceFileDir() { return referenceFileDir; }
    protected final int[] getClipPoints() { return clipPoints; }
    protected final Integer getExpectedInsertSize() { return expectedInsertSize; }
    protected final Integer getReadsToAlign() { return readsToAlign; }
    protected final Map<String, String> getCustomParametersMap() { return customParametersMap; }
    protected final boolean isPairedReads() { return pairedReads; }
    protected final int getReadLength() { return readLength; }
}
