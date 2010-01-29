package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.sam.SAMReadViolationHistogram;

import java.io.File;
import java.util.Collection;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 2:36:16 PM
 * <p/>
 * Converts shards to SAM iterators over the specified region
 */
public abstract class SAMDataSource implements SimpleDataSource {

    /** Backing support for reads. */
    protected final Reads reads;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMDataSource.class);

    // do we take unmapped reads
    protected boolean includeUnmappedReads = true;

    /**
     * A histogram of exactly what reads were removed from the input stream and why.
     */
    protected SAMReadViolationHistogram violations = new SAMReadViolationHistogram();

    /**
     * Returns a histogram of reads that were screened out, grouped by the nature of the error.
     * @return Histogram of reads.  Will not be null.
     */
    public SAMReadViolationHistogram getViolationHistogram() {
        return violations;
    }

    /**
     * constructor, given sam files
     *
     * @param reads the list of sam files
     */
    public SAMDataSource( Reads reads ) throws SimpleDataSourceLoadException {
        this.reads = reads;

        // check the length
        if (reads.getReadsFiles().size() < 1) {
            throw new SimpleDataSourceLoadException("SAMDataSource: you must provide a list of length greater then 0");
        }
        for (File smFile : reads.getReadsFiles()) {
            if (!smFile.canRead()) {
                throw new SimpleDataSourceLoadException("SAMDataSource: Unable to load file: " + smFile.getName());
            }
        }
    }

    /**
     * Do all BAM files backing this data source have an index?  The case where hasIndex() is false
     * is supported, but only in a few extreme cases.
     * @return True if an index is present; false otherwise.
     */
    public abstract boolean hasIndex();

    /**
     * Gets the (potentially merged) SAM file header.
     *
     * @return SAM file header.
     */
    public abstract SAMFileHeader getHeader();


    /**
     * Returns Reads data structure containing information about the reads data sources placed in this pool as well as
     * information about how they are downsampled, sorted, and filtered
     * @return
     */
    public Reads getReadsInfo() { return reads; }

    /**
     * Returns readers used by this data source.
     */
    public abstract Collection<SAMFileReader> getReaders();

    /** Returns true if there are read group duplicates within the merged headers. */
    public abstract boolean hasReadGroupCollisions();

    /** Returns the read group id that should be used for the input read and RG id. */
    public abstract String getReadGroupId(final SAMFileReader reader, final String originalReadGroupId);

    /**
     *
     * @param shard the shard to get data for
     *
     * @return an iterator for that region
     */
    public abstract StingSAMIterator seek(Shard shard);

    /**
     * If we're in by-read mode, this indicates if we want
     * to see unmapped reads too.  Only seeing mapped reads
     * is much faster, but most BAM files have significant
     * unmapped read counts.
     *
     * @param seeUnMappedReads true to see unmapped reads, false otherwise
     */
    public void viewUnmappedReads( boolean seeUnMappedReads ) {
        includeUnmappedReads = seeUnMappedReads;
    }

    /**
     * Filter reads based on user-specified criteria.
     *
     * @param enableVerification Verify the order of reads.
     * @param wrappedIterator the raw data source.
     * @param downsamplingFraction whether and how much to downsample the reads themselves (not at a locus).
     * @param noValidationOfReadOrder Another trigger for the verifying iterator?  TODO: look into this.
     * @param supplementalFilters additional filters to apply to the reads.
     * @return An iterator wrapped with filters reflecting the passed-in parameters.  Will not be null.
     */
    protected StingSAMIterator applyDecoratingIterators(boolean enableVerification,
                                                        StingSAMIterator wrappedIterator,
                                                        Double downsamplingFraction,
                                                        Boolean noValidationOfReadOrder,
                                                        Collection<SamRecordFilter> supplementalFilters) {
        // NOTE: this (and other filtering) should be done before on-the-fly sorting
        //  as there is no reason to sort something that we will end of throwing away
        if (downsamplingFraction != null)
            wrappedIterator = new DownsampleIterator(wrappedIterator, downsamplingFraction);

        // unless they've said not to validate read ordering (!noValidationOfReadOrder) and we've enabled verification,
        // verify the read ordering by applying a sort order iterator
        if (!noValidationOfReadOrder && enableVerification)
            wrappedIterator = new VerifyingSamIterator(wrappedIterator);

        for( SamRecordFilter supplementalFilter: supplementalFilters )
            wrappedIterator = StingSAMIteratorAdapter.adapt(wrappedIterator.getSourceInfo(),
                                                            new FilteringIterator(wrappedIterator,supplementalFilter));

        return wrappedIterator;
    }
}


