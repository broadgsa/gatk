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

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;

/**
 * Allows the user to steadily accumulate information about what
 * components go into a SAM file writer, ultimately using this
 * information to create a SAM file writer on demand.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMFileReaderBuilder {
    /**
     * To which file should output be written?
     */
    private File samFile = null;

    /**
     * What compression level should be used when building this file?
     */
    private SAMFileReader.ValidationStringency validationStringency = null;

    /**
     * Sets the handle of the sam file to which data should be written.
     * @param samFile The SAM file into which data should flow.
     */
    public void setSAMFile( File samFile ) {
        this.samFile = samFile;
    }

    /**
     * Sets the validation stringency to apply when reading this sam file.
     * @param validationStringency Stringency to apply.  Must not be null.
     */
    public void setValidationStringency( SAMFileReader.ValidationStringency validationStringency ) {
        this.validationStringency = validationStringency;
    }

    /**
     * Create the SAM writer, given the constituent parts accrued.
     * @return Newly minted SAM file writer.
     */
    public SAMFileReader build() {
        if( samFile == null )
            throw new ReviewedStingException( "Filename for output sam file must be supplied.");
        if( validationStringency == null )
            throw new ReviewedStingException( "Header for output sam file must be supplied.");

        SAMFileReader reader = new SAMFileReader( samFile );
        reader.setValidationStringency( validationStringency );

        return reader;
    }
}
