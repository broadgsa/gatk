/*
 * Copyright (c) 2012 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.bqsr;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.recalibration.ReadCovariates;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 12/18/12
 * Time: 3:50 PM
 *
 * TODO -- merge in ReadCovariates?
 */
public final class ReadRecalibrationInfo {
    private final GATKSAMRecord read;
    private final int length;
    private final ReadCovariates covariates;
    private final boolean[] skips;
    private final byte[] baseQuals, insertionQuals, deletionQuals;
    private final double[] snpErrors, insertionErrors, deletionErrors;

    public ReadRecalibrationInfo(final GATKSAMRecord read,
                                 final ReadCovariates covariates,
                                 final boolean[] skips,
                                 final double[] snpErrors,
                                 final double[] insertionErrors,
                                 final double[] deletionErrors) {
        if ( read == null ) throw new IllegalArgumentException("read cannot be null");
        if ( covariates == null ) throw new IllegalArgumentException("covariates cannot be null");
        if ( skips == null ) throw new IllegalArgumentException("skips cannot be null");
        if ( snpErrors == null ) throw new IllegalArgumentException("snpErrors cannot be null");
        // future: may allow insertionErrors && deletionErrors to be null, so don't enforce

        this.read = read;
        this.baseQuals = read.getBaseQualities();
        this.length = baseQuals.length;
        this.covariates = covariates;
        this.skips = skips;
        this.insertionQuals = read.getExistingBaseInsertionQualities();
        this.deletionQuals = read.getExistingBaseDeletionQualities();
        this.snpErrors = snpErrors;
        this.insertionErrors = insertionErrors;
        this.deletionErrors = deletionErrors;

        if ( skips.length != length ) throw new IllegalArgumentException("skips.length " + snpErrors.length + " != length " + length);
        if ( snpErrors.length != length ) throw new IllegalArgumentException("snpErrors.length " + snpErrors.length + " != length " + length);
        if ( insertionErrors != null && insertionErrors.length != length ) throw new IllegalArgumentException("insertionErrors.length " + snpErrors.length + " != length " + length);
        if ( deletionErrors != null && deletionErrors.length != length ) throw new IllegalArgumentException("deletionErrors.length " + snpErrors.length + " != length " + length);
    }

    /**
     * Get the qual score for event type at offset
     *
     * @param eventType the type of event we want the qual for
     * @param offset the offset into this read for the qual
     * @return a valid quality score for event at offset
     */
    @Requires("validOffset(offset)")
    @Ensures("validQual(result)")
    public byte getQual(final EventType eventType, final int offset) {
        switch ( eventType ) {
            case BASE_SUBSTITUTION: return baseQuals[offset];
            // note optimization here -- if we don't have ins/del quals we just return the default byte directly
            case BASE_INSERTION: return insertionQuals == null ? GATKSAMRecord.DEFAULT_INSERTION_DELETION_QUAL : insertionQuals[offset];
            case BASE_DELETION: return deletionQuals == null ? GATKSAMRecord.DEFAULT_INSERTION_DELETION_QUAL : deletionQuals[offset];
            default: throw new IllegalStateException("Unknown event type " + eventType);
        }
    }

    /**
     * Get the error fraction for event type at offset
     *
     * The error fraction is a value between 0 and 1 that indicates how much certainty we have
     * in the error occurring at offset.  A value of 1 means that the error definitely occurs at this
     * site, a value of 0.0 means it definitely doesn't happen here.  0.5 means that half the weight
     * of the error belongs here
     *
     * @param eventType the type of event we want the qual for
     * @param offset the offset into this read for the qual
     * @return a fractional weight for an error at this offset
     */
    @Requires("validOffset(offset)")
    @Ensures("result >= 0.0 && result <= 1.0")
    public double getErrorFraction(final EventType eventType, final int offset) {
        switch ( eventType ) {
            case BASE_SUBSTITUTION: return snpErrors[offset];
            case BASE_INSERTION: return insertionErrors[offset];
            case BASE_DELETION: return deletionErrors[offset];
            default: throw new IllegalStateException("Unknown event type " + eventType);
        }
    }

    /**
     * Get the read involved in this recalibration info
     * @return a non-null GATKSAMRecord
     */
    @Ensures("result != null")
    public GATKSAMRecord getRead() {
        return read;
    }

    /**
     * Should offset in this read be skipped (because it's covered by a known variation site?)
     * @param offset a valid offset into this info
     * @return true if offset should be skipped, false otherwise
     */
    @Requires("validOffset(offset)")
    public boolean skip(final int offset) {
        return skips[offset];
    }

    /**
     * Get the ReadCovariates object carrying the mapping from offsets -> covariate key sets
     * @return a non-null ReadCovariates object
     */
    @Ensures("result != null")
    public ReadCovariates getCovariatesValues() {
        return covariates;
    }

    /**
     * Ensures an offset is valid.  Used in contracts
     * @param offset a proposed offset
     * @return true if offset is valid w.r.t. the data in this object, false otherwise
     */
    private boolean validOffset(final int offset) {
        return offset >= 0 && offset < baseQuals.length;
    }

    private boolean validQual(final byte result) {
        return result >= 0 && result <= QualityUtils.MAX_QUAL_SCORE;
    }
}
