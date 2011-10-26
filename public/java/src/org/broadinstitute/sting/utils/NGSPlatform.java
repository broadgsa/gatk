/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

/**
 * A canonical, master list of the standard NGS platforms.  These values
 * can be obtained (efficiently) from a GATKSAMRecord object with the
 * getNGSPlatform method.
 *
 * @author Mark DePristo
 * @since 2011
 */
public enum NGSPlatform {
    ILLUMINA("ILLUMINA", "SLX", "SOLEXA"),
    SOLID("SOLID"),
    LS454("454"),
    COMPLETE_GENOMICS("COMPLETE"),
    PACBIO("PACBIO"),
    ION_TORRENT("IONTORRENT"),
    UNKNOWN("UNKNOWN");

    /**
     * Array of the prefix names in a BAM file for each of the platforms.
     */
    private final String[] BAM_PL_NAMES;

    NGSPlatform(final String... BAM_PL_NAMES) {
        for ( int i = 0; i < BAM_PL_NAMES.length; i++ )
            BAM_PL_NAMES[i] = BAM_PL_NAMES[i].toUpperCase();
        this.BAM_PL_NAMES = BAM_PL_NAMES;
    }

    /**
     * Returns a representative PL string for this platform
     * @return
     */
    public final String getDefaultPlatform() {
        return BAM_PL_NAMES[0];
    }

    /**
     * Convenience constructor -- calculates the NGSPlatfrom from a SAMRecord.
     * Note you should not use this function if you have a GATKSAMRecord -- use the
     * accessor method instead.
     *
     * @param read
     * @return an NGSPlatform object matching the PL field of the header, of UNKNOWN if there was no match
     */
    public static final NGSPlatform fromRead(SAMRecord read) {
        return fromReadGroup(read.getReadGroup());
    }

    /**
     * Returns the NGSPlatform corresponding to the PL tag in the read group
     * @param rg
     * @return an NGSPlatform object matching the PL field of the header, of UNKNOWN if there was no match
     */
    public static final NGSPlatform fromReadGroup(SAMReadGroupRecord rg) {
        return fromReadGroupPL(rg.getPlatform());
    }

    /**
     * Returns the NGSPlatform corresponding to the PL tag in the read group
     * @param plFromRG -- the PL field (or equivalent) in a ReadGroup object
     * @return an NGSPlatform object matching the PL field of the header, of UNKNOWN if there was no match
     */
    public static final NGSPlatform fromReadGroupPL(final String plFromRG) {
        if ( plFromRG == null ) return UNKNOWN;

        // todo -- algorithm could be implemented more efficiently, as the list of all
        // todo -- names is known upfront, so a decision tree could be used to identify
        // todo -- a prefix common to PL
        final String pl = plFromRG.toUpperCase();
        for ( final NGSPlatform ngsPlatform : NGSPlatform.values() ) {
            for ( final String bamPLName : ngsPlatform.BAM_PL_NAMES ) {
                if ( pl.contains(bamPLName) )
                    return ngsPlatform;
            }
        }

        return UNKNOWN;
    }
}
