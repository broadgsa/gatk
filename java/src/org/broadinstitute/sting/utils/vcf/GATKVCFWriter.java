/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.utils.vcf;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broad.tribble.vcf.*;

import java.io.*;


/**
 * @author ebanks
 *         <p/>
 *         Class GATKVCFWriter
 *         <p/>
 *         GATK-specific version of the VCF Writer
 */
public class GATKVCFWriter extends VCFWriter implements VCFGenotypeWriter {

    public GATKVCFWriter(File writeTo) {
        super(writeTo);
    }

    public GATKVCFWriter(OutputStream writeTo) {
        super(writeTo);
    }

    public void writeHeader(VCFHeader header) {
        // TODO -- put the command-line generating code for the header right here
        super.writeHeader(header);
    }

    public void append(File file) {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line = reader.readLine();
            while ( line != null ) {
                if ( !VCFHeaderLine.isHeaderLine(line) ) {
                    mWriter.write(line);
                    mWriter.write("\n");
                }
                line = reader.readLine();
            }

            reader.close();
        } catch (IOException e) {
            throw new StingException("Error reading file " + file + " in GATKVCFWriter: ", e);
        }
    }
}
