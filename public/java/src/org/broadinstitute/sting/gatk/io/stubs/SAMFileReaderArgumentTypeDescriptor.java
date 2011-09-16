/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk.io.stubs;

import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.commandline.ArgumentMatches;
import org.broadinstitute.sting.commandline.ArgumentSource;
import org.broadinstitute.sting.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.sting.commandline.ParsingEngine;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.SAMFileReaderBuilder;

import java.io.File;
import java.lang.reflect.Type;

/**
 * Describe how to parse SAMFileReaders.
 */
public class SAMFileReaderArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * The engine into which output stubs should be fed.
     */
    private GenomeAnalysisEngine engine;

    /**
     * Create a new SAMFileReader argument, notifying the given engine when that argument has been created.
     * @param engine
     */
    public SAMFileReaderArgumentTypeDescriptor( GenomeAnalysisEngine engine ) {
        this.engine = engine;
    }

    @Override
    public boolean supports( Class type ) {
        return SAMFileReader.class.isAssignableFrom(type);
    }

    @Override
    public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches ) {
        SAMFileReaderBuilder builder = new SAMFileReaderBuilder();

        String readerFileName = getArgumentValue( createDefaultArgumentDefinition(source), matches );

        if( readerFileName == null )
            throw new UserException.CommandLineException("SAM file compression was supplied, but no associated writer was supplied with it.");

        builder.setSAMFile(new File(readerFileName));

        // WARNING: Skipping required side-effect because stub is impossible to generate.
        engine.addInput(source, builder);

        // MASSIVE KLUDGE!  SAMFileReader is tricky to implement and we don't yet have a stub.  Return null, then
        // let the output tracker load it in.
        // TODO: Add a stub for SAMFileReader.
        return null;
    }
}
