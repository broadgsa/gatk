/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.utils.codecs.bcf2;

import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.writer.Options;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriterFactory;

import java.io.*;
import java.util.*;

/**
 * Testing BCF2
 *
 * @author Mark DePristo
 * @since 2012
 */
public class BCF2TestWalker extends RodWalker<Integer, Integer> {
    /**
     * Variants from this VCF file are used by this tool as input.
     * The file must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;

    @Argument(doc="keep variants", required=false)
    public boolean keepVariants = false;

    @Argument(doc="quiet", required=false)
    public boolean quiet = false;

    @Argument(doc="dontIndexOnTheFly", required=false)
    public boolean dontIndexOnTheFly = false;

    @Output(doc="File to which results should be written",required=true)
    protected File bcfFile;

    private final List<VariantContext> vcs = new ArrayList<VariantContext>();
    protected VariantContextWriter writer;

    @Override
    public void initialize() {
        final Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), Collections.singletonList(variants));
        final VCFHeader header = VCFUtils.withUpdatedContigs(vcfRods.values().iterator().next(), getToolkit());
        try {
            EnumSet<Options> options = EnumSet.of(Options.FORCE_BCF);
            if ( !dontIndexOnTheFly ) options.add(Options.INDEX_ON_THE_FLY);
            writer = VariantContextWriterFactory.create(bcfFile, new FileOutputStream(bcfFile), getToolkit().getMasterSequenceDictionary(), options);
            writer.writeHeader(header);
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(bcfFile, e);
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        for ( VariantContext vc : tracker.getValues(variants, context.getLocation())) {
            writer.add(vc);
            if ( keepVariants ) vcs.add(vc);
        }

        return 1;
    }

    //
    // default reduce -- doesn't do anything at all
    //
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer counter, Integer sum) { return counter + sum; }

    public void onTraversalDone(Integer sum) {
        try {
            writer.close();
            logger.info("Closed writer");

            // read in the BCF records
            BCF2Codec codec = new BCF2Codec();
            PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(bcfFile));
            FeatureCodecHeader header = codec.readHeader(pbs);
            pbs.close();

            pbs = new PositionalBufferedStream(new FileInputStream(bcfFile));
            pbs.skip(header.getHeaderEnd());
            Iterator<VariantContext> it = vcs.iterator();
            while ( ! pbs.isDone() ) {
                if ( keepVariants ) {
                    VariantContext expected = it.next();
                    if ( ! quiet )
                        System.out.printf("vcf = %s %d %s%n", expected.getChr(), expected.getStart(), expected);
                }
                VariantContext bcfRaw = codec.decode(pbs);
                VariantContext bcf = new VariantContextBuilder(bcfRaw).source("variant").make();
                if ( ! quiet ) {
                    System.out.printf("bcf = %s %d %s%n", bcf.getChr(), bcf.getStart(), bcf.toString());
                    System.out.printf("--------------------------------------------------%n");
                }
            }

        } catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(bcfFile, "bad user!");
        }
    }
}
