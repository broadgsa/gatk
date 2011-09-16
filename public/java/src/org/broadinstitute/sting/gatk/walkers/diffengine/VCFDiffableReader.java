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

package org.broadinstitute.sting.gatk.walkers.diffengine;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;


/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 1:09 PM
 *
 * Class implementing diffnode reader for VCF
 */
public class VCFDiffableReader implements DiffableReader {
    @Override
    public String getName() { return "VCF"; }

    @Override
    public DiffElement readFromFile(File file, int maxElementsToRead) {
        DiffNode root = DiffNode.rooted(file.getName());
        try {
            // read the version line from the file
            LineReader lineReader = new AsciiLineReader(new FileInputStream(file));
            final String version = lineReader.readLine();
            root.add("VERSION", version);
            lineReader.close();

            lineReader = new AsciiLineReader(new FileInputStream(file));
            VCFCodec vcfCodec = new VCFCodec();

            // must be read as state is stored in reader itself
            VCFHeader header = (VCFHeader)vcfCodec.readHeader(lineReader);
            for ( VCFHeaderLine headerLine : header.getMetaData() ) {
                String key = headerLine.getKey();
                if ( headerLine instanceof VCFNamedHeaderLine )
                    key += "_" + ((VCFNamedHeaderLine) headerLine).getName();
                root.add(key, headerLine.toString());
            }

            String line = lineReader.readLine();
            int count = 0, nRecordsAtPos = 1;
            String prevName = "";
            while ( line != null ) {
                if ( count++ > maxElementsToRead && maxElementsToRead != -1)
                    break;

                VariantContext vc = (VariantContext)vcfCodec.decode(line);
                String name = vc.getChr() + ":" + vc.getStart();
                if ( name.equals(prevName) ) {
                    name += "_" + ++nRecordsAtPos;
                } else {
                    prevName = name;
                }
                DiffNode vcRoot = DiffNode.empty(name, root);

                // add fields
                vcRoot.add("CHROM", vc.getChr());
                vcRoot.add("POS", vc.getStart());
                vcRoot.add("ID", vc.hasID() ? vc.getID() : VCFConstants.MISSING_VALUE_v4);
                vcRoot.add("REF", vc.getReference());
                vcRoot.add("ALT", vc.getAlternateAlleles());
                vcRoot.add("QUAL", vc.hasNegLog10PError() ? vc.getNegLog10PError() * 10 : VCFConstants.MISSING_VALUE_v4);
                vcRoot.add("FILTER", vc.getFilters());

                // add info fields
                for (Map.Entry<String, Object> attribute : vc.getAttributes().entrySet()) {
                    if ( ! attribute.getKey().startsWith("_") && ! attribute.getKey().equals(VariantContext.ID_KEY))
                        vcRoot.add(attribute.getKey(), attribute.getValue());
                }

                for (Genotype g : vc.getGenotypes().values() ) {
                    DiffNode gRoot = DiffNode.empty(g.getSampleName(), vcRoot);
                    gRoot.add("GT", g.getGenotypeString());
                    gRoot.add("GQ", g.hasNegLog10PError() ? g.getNegLog10PError() * 10 : VCFConstants.MISSING_VALUE_v4 );

                    for (Map.Entry<String, Object> attribute : g.getAttributes().entrySet()) {
                        if ( ! attribute.getKey().startsWith("_") )
                            gRoot.add(attribute.getKey(), attribute.getValue());
                    }

                    vcRoot.add(gRoot);
                }

                root.add(vcRoot);
                line = lineReader.readLine();
            }

            lineReader.close();
        } catch ( IOException e ) {
            return null;
        }

        return root.getBinding();
    }

    @Override
    public boolean canRead(File file) {
        return AbstractVCFCodec.canDecodeFile(file, VCFCodec.VCF4_MAGIC_HEADER);
    }
}
