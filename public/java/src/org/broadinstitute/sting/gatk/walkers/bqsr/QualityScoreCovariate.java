package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;

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
 * Date: Nov 3, 2009
 *
 * The Reported Quality Score covariate.
 */

public class QualityScoreCovariate implements RequiredCovariate {

    private byte defaultMismatchesQuality;                                      // walker parameter. Must be > 0 to be used, otherwise we use the quality from the read.
    private byte defaultInsertionsQuality;                                      // walker parameter. Must be > 0 to be used, otherwise we use the quality from the read.
    private byte  defaultDeletionsQuality;                                      // walker parameter. Must be > 0 to be used, otherwise we use the quality from the read.
        
    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        defaultMismatchesQuality = RAC.MISMATCHES_DEFAULT_QUALITY;
        defaultInsertionsQuality = RAC.INSERTIONS_DEFAULT_QUALITY;
         defaultDeletionsQuality = RAC.DELETIONS_DEFAULT_QUALITY; 
    }

    @Override
    public CovariateValues getValues(final GATKSAMRecord read) {
        int readLength = read.getReadLength();
        
        Byte [] mismatches = new Byte[readLength];
        Byte [] insertions = new Byte[readLength];
        Byte []  deletions = new Byte[readLength];
        
        byte [] baseQualities = read.getBaseQualities();

        if (defaultMismatchesQuality >= 0)
            Arrays.fill(mismatches, defaultMismatchesQuality);                  // if the user decides to override the base qualities in the read, use the flat value
        else {
            for (int i=0; i<baseQualities.length; i++)
                mismatches[i] = baseQualities[i];
        }

        Arrays.fill(insertions, defaultInsertionsQuality);                      // Some day in the future when base insertion and base deletion quals exist the samtools API will
        Arrays.fill( deletions, defaultDeletionsQuality);                       // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat value (parameter)

        return new CovariateValues(mismatches, insertions, deletions);
    }

}
