/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.gatk.utils.NGSPlatform;

/**
 * @author ebanks
 * GATKSAMReadGroupRecord
 *
 * this class extends the samtools SAMReadGroupRecord class and caches important
 * (and oft-accessed) data that's not already cached by the SAMReadGroupRecord class
 *
 */
public class GATKSAMReadGroupRecord extends SAMReadGroupRecord {
    // the SAMReadGroupRecord data we're caching
    private String mSample = null;
    private String mPlatform = null;
    private NGSPlatform mNGSPlatform = null;

    // because some values can be null, we don't want to duplicate effort
    private boolean retrievedSample = false;
    private boolean retrievedPlatform = false;
    private boolean retrievedNGSPlatform = false;

    public GATKSAMReadGroupRecord(final String id) {
        super(id);
    }

    public GATKSAMReadGroupRecord(SAMReadGroupRecord record) {
        super(record.getReadGroupId(), record);
    }

    /**
     * Get the NGSPlatform enum telling us the platform of this read group
     *
     * This function call is caching, so subsequent calls to it are free, while
     * the first time it's called there's a bit of work to resolve the enum
     *
     * @return an NGSPlatform enum value
     */
    public NGSPlatform getNGSPlatform() {
        if ( ! retrievedNGSPlatform ) {
            mNGSPlatform = NGSPlatform.fromReadGroupPL(getPlatform());
            retrievedNGSPlatform = true;
        }

        return mNGSPlatform;
    }

    @Override
    public String toString() {
        return "GATKSAMReadGroupRecord @RG:" + getReadGroupId();
    }

    ///////////////////////////////////////////////////////////////////////////////
    // *** The following methods are overloaded to cache the appropriate data ***//
    ///////////////////////////////////////////////////////////////////////////////

    @Override
    public String getSample() {
        if ( !retrievedSample ) {
            mSample = super.getSample();
            retrievedSample = true;
        }
        return mSample;
    }

    @Override
    public void setSample(String s) {
        super.setSample(s);
        mSample = s;
        retrievedSample = true;
    }

    @Override
    public String getPlatform() {
        if ( !retrievedPlatform ) {
            mPlatform = super.getPlatform();
            retrievedPlatform = true;
        }
        return mPlatform;
    }

    @Override
    public void setPlatform(String s) {
        super.setPlatform(s);
        mPlatform = s;
        retrievedPlatform = true;
        retrievedNGSPlatform = false; // recalculate the NGSPlatform
    }
}