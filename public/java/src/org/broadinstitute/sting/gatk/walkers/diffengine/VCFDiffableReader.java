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

import org.apache.log4j.Logger;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;
import java.util.Arrays;
import java.util.Iterator;
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
    private static Logger logger = Logger.getLogger(VCFDiffableReader.class);

    @Override
    public String getName() { return "VCF"; }

    @Override
    public DiffElement readFromFile(File file, int maxElementsToRead) {
        DiffNode root = DiffNode.rooted(file.getName());
        try {
            // read the version line from the file
            BufferedReader br = new BufferedReader(new FileReader(file));
            final String version = br.readLine();
            root.add("VERSION", version);
            br.close();

            // must be read as state is stored in reader itself
            FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new VCFCodec(), false);
            VCFHeader header = (VCFHeader)reader.getHeader();
            for ( VCFHeaderLine headerLine : header.getMetaData() ) {
                String key = headerLine.getKey();
                if ( headerLine instanceof VCFIDHeaderLine)
                    key += "_" + ((VCFIDHeaderLine) headerLine).getID();
                if ( root.hasElement(key) )
                    logger.warn("Skipping duplicate header line: file=" + file + " line=" + headerLine.toString());
                else
                    root.add(key, headerLine.toString());
            }

            int count = 0, nRecordsAtPos = 1;
            String prevName = "";
            Iterator<VariantContext> it = reader.iterator();
            while ( it.hasNext() ) {
                if ( count++ > maxElementsToRead && maxElementsToRead != -1)
                    break;

                VariantContext vc = it.next();
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
                vcRoot.add("ID", vc.getID());
                vcRoot.add("REF", vc.getReference());
                vcRoot.add("ALT", vc.getAlternateAlleles());
                vcRoot.add("QUAL", vc.hasLog10PError() ? vc.getLog10PError() * -10 : VCFConstants.MISSING_VALUE_v4);
                vcRoot.add("FILTER", vc.getFilters());

                // add info fields
                for (Map.Entry<String, Object> attribute : vc.getAttributes().entrySet()) {
                    if ( ! attribute.getKey().startsWith("_") )
                        vcRoot.add(attribute.getKey(), attribute.getValue());
                }

                for (Genotype g : vc.getGenotypes() ) {
                    DiffNode gRoot = DiffNode.empty(g.getSampleName(), vcRoot);
                    gRoot.add("GT", g.getGenotypeString());
                    if ( g.hasGQ() ) gRoot.add("GQ", g.getGQ() );
                    if ( g.hasDP() ) gRoot.add("DP", g.getDP() );
                    if ( g.hasAD() ) gRoot.add("AD", Utils.join(",", g.getAD()));
                    if ( g.hasPL() ) gRoot.add("PL", Utils.join(",", g.getPL()));

                    for (Map.Entry<String, Object> attribute : g.getExtendedAttributes().entrySet()) {
                        if ( ! attribute.getKey().startsWith("_") )
                            gRoot.add(attribute.getKey(), attribute.getValue());
                    }

                    vcRoot.add(gRoot);
                }

                root.add(vcRoot);
            }

            reader.close();
        } catch ( IOException e ) {
            return null;
        }

        return root.getBinding();
    }

    @Override
    public boolean canRead(File file) {
        return AbstractVCFCodec.canDecodeFile(file.getPath(), VCFCodec.VCF4_MAGIC_HEADER);
    }
}
