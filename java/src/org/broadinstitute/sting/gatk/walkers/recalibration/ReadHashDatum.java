package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import edu.mit.broad.picard.illumina.parser.IlluminaUtil;

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
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 14, 2009
 */

public class ReadHashDatum {
    public String readGroup;
    public String platform;
    public byte[] quals;
    public byte[] bases;
    public boolean isNegStrand;
    public int mappingQuality;
    public int length;
    public Integer tile;
    private static boolean warnUserNullReadGroup = false; // Has the walker warned the user about null read groups yet?

    public ReadHashDatum(String _readGroup, String _platform, byte[] _quals, byte[] _bases, boolean _isNegStrand,
                         int _mappingQuality, int _length, Integer _tile) {
        readGroup = _readGroup;
        platform = _platform;
        quals = _quals;
        bases = _bases;
        isNegStrand = _isNegStrand;
        mappingQuality = _mappingQuality;
        length = _length;
        tile = _tile;
    }

    public ReadHashDatum(ReadHashDatum that) {
        this.readGroup = that.readGroup;
        this.platform = that.platform;
        this.quals = that.quals;
        this.bases = that.bases;
        this.isNegStrand = that.isNegStrand;
        this.mappingQuality = that.mappingQuality;
        this.length = that.length;
        this.tile = that.tile;
    }

    /**
     * Pulls the important bits out of a SAMRecord and stuffs them into a leaner Object to be used by the Covariates and cached by the Recalibration walkers
     * @param read The SAMRecord whose data we want to get at
     * @param RAC The collection of arguments that are shared between CovariateCounterWalker and TableRecalibrationWalker
     * @return A new ReadHashDatum holding all the important parts of the read
     */
    public static ReadHashDatum parseSAMRecord( final SAMRecord read, final RecalibrationArgumentCollection RAC ) {

        // Lots of lines just pulling things out of the SAMRecord and stuffing into local variables, many edge cases to worry about
        byte[] quals = read.getBaseQualities();

        // Check if we need to use the original quality scores instead
        if ( RAC.USE_ORIGINAL_QUALS && read.getAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG) != null ) {
            Object obj = read.getAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG);
            if ( obj instanceof String )
                quals = QualityUtils.fastqToPhred((String)obj);
            else {
                throw new RuntimeException(String.format("Value encoded by %s in %s isn't a string!", RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG, read.getReadName()));
            }
        }

        final byte[] bases = read.getReadBases(); // BUGBUG: DinucCovariate is relying on this method returning the same byte for bases 'a' and 'A'.
        final boolean isNegStrand = read.getReadNegativeStrandFlag();
        final SAMReadGroupRecord readGroup = read.getReadGroup();
        String readGroupId;
        String platform;

        // If there are no read groups we have to default to something, and that something could be specified by the user using command line arguments
        if( readGroup == null ) {
            if( !warnUserNullReadGroup && RAC.FORCE_READ_GROUP == null ) {
                Utils.warnUser("The input .bam file contains reads with no read group. " +
                                "Defaulting to read group ID = " + RAC.DEFAULT_READ_GROUP + " and platform = " + RAC.DEFAULT_PLATFORM + ". " +
                                "First observed at read with name = " + read.getReadName() );
                warnUserNullReadGroup = true;
            }
            // There is no readGroup so defaulting to these values
            readGroupId = RAC.DEFAULT_READ_GROUP;
            platform = RAC.DEFAULT_PLATFORM;
        } else {
            readGroupId = readGroup.getReadGroupId();
            platform = readGroup.getPlatform();
        }

        final int mappingQuality = read.getMappingQuality();
        final int length = bases.length;

        if( RAC.FORCE_READ_GROUP != null ) { // Collapse all the read groups into a single common String provided by the user
            readGroupId = RAC.FORCE_READ_GROUP;
        }

        if( RAC.FORCE_PLATFORM != null ) {
            platform = RAC.FORCE_PLATFORM;
        }

        Integer tile = IlluminaUtil.getTileFromReadName(read.getReadName());

        return new ReadHashDatum( readGroupId, platform, quals, bases, isNegStrand, mappingQuality, length, tile );

    }
}
