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

package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.TribbleRMDTrackBuilder;
import org.broadinstitute.sting.utils.StingException;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFFormatHeaderLine;
import org.broad.tribble.vcf.VCFCodec;
import org.broad.tribble.source.BasicFeatureSource;
import org.broad.tribble.index.Index;

import java.io.*;
import java.util.HashSet;
import java.util.Iterator;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFGenotypeWriterAdapter
 *         <p/>
 *         Adapt the VCF writter to the genotype output system
 */
public class VCFGenotypeWriterAdapter extends VCFWriter implements VCFGenotypeWriter {

    // allowed genotype format strings
    private HashSet<String> allowedGenotypeFormatKeys = null;

    public VCFGenotypeWriterAdapter(File writeTo) {
        super(writeTo);
    }

    public VCFGenotypeWriterAdapter(OutputStream writeTo) {
        super(writeTo);
    }

    public void addCall(VariantContext vc, byte ref) {
        vc = VariantContextUtils.purgeUnallowedGenotypeAttributes(vc, allowedGenotypeFormatKeys);
        add(vc, ref);
    }

    public void writeHeader(VCFHeader header) {
        for ( VCFHeaderLine field : header.getMetaData() ) {
            if ( field instanceof VCFFormatHeaderLine ) {
                if ( allowedGenotypeFormatKeys == null )
                    allowedGenotypeFormatKeys = new HashSet<String>();
                allowedGenotypeFormatKeys.add(((VCFFormatHeaderLine)field).getName());
            }
        }
        super.writeHeader(header);
    }

    public void append(File file) {
        // if we don't need to restrict the FORMAT fields, then read blindly
        if ( allowedGenotypeFormatKeys == null )
            blindAppend(file);
        else
            smartAppend(file);
    }

    private void blindAppend(File file) {
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
            throw new StingException("Error reading file " + file + " in VCFGenotypeWriter: ", e);
        }
    }

    private void smartAppend(File file) {
        try {
            VCFCodec codec = new VCFCodec();
            Index index = TribbleRMDTrackBuilder.loadIndex(file, codec, false);
            BasicFeatureSource<VariantContext> vcfReader = new BasicFeatureSource(file.getAbsolutePath(),index,codec);
            Iterator<VariantContext> iterator = vcfReader.iterator();
            while ( iterator.hasNext() ) {
                VariantContext vc = iterator.next();
                vc = VariantContextUtils.purgeUnallowedGenotypeAttributes(vc, allowedGenotypeFormatKeys);
                add(vc, vc.getReferenceBaseForIndel());
            }
            vcfReader.close();
        } catch (FileNotFoundException e) {
            throw new StingException("Error reading file " + file + " in VCFGenotypeWriter: ", e);
        } catch (IOException e) {
            throw new StingException("Error reading file " + file + " in VCFGenotypeWriter: ", e);
        }



    }
}
