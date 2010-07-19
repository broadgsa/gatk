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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.vcf.*;
import org.broad.tribble.util.ParsingUtils;
import org.broad.tribble.FeatureCodec;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.features.vcf4.VCF4Codec;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.commandline.Argument;

import java.util.*;
import java.io.*;

import com.sun.xml.internal.ws.wsdl.parser.ParserUtil;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Apr 13, 2010
 */
public class VCF4ReaderTestWalker extends RodWalker<VCFRecord,Long> {
    @Argument(shortName="MR", doc="", required=false)
    int maxRecords = -1;
    @Argument(shortName="vcf", doc="", required=true)
    File vcfFile = null;
    @Argument(shortName="Parse", doc="", required=true)
    ParsingStatus splitFile = ParsingStatus.NONE;
    @Argument(shortName="DontValidate", doc="", required=false)
    boolean DontValidate = false;

    @Argument(shortName="USE_VCF3", doc="", required=false)
    boolean USE_VCF3 = false;


    public enum ParsingStatus { NONE, SPLIT_LINES, VARIANTS, GENOTYPES }

    public void initialize() {
    }

    public VCFRecord map(RefMetaDataTracker tracker, ReferenceContext context, AlignmentContext alicon) {
        return null;
    }

    public Long reduce(VCFRecord con, Long num) {
        if ( con == null ) {
            return num;
        }

        return 1 + num;
    }

    public Long reduceInit() {
        return 0l;
    }

    String[] parts = new String[10000];
    public void onTraversalDone(Long num){
        VCF4Codec vcf4codec = new VCF4Codec();
        //VCF4Codec.parseGenotypesToo = splitFile == ParsingStatus.GENOTYPES;
        VCF4Codec.validate = ! DontValidate;

        VCFCodec  vcf3codec = new VCFCodec();

        FeatureCodec codec = USE_VCF3 ? vcf3codec : vcf4codec;

        try {
            AsciiLineReader lineReader = new AsciiLineReader(new FileInputStream(vcfFile));
            VCFHeader header = (VCFHeader)codec.readHeader(lineReader);
            out.printf("Read %d header lines%n", header.getMetaData().size()+1);

            // a counter of the number of lines we've read
            int lineNumber = header.getMetaData().size()+1;
            while (true) {
                String line = lineReader.readLine();

                if ( line == null )
                    break;

                lineNumber++;
                if ( lineNumber >= maxRecords && maxRecords != -1 ) {
                    return;
                }

                if ( line.charAt(0) == '#' ) 
                    continue;

                Object vc = null;
                if ( splitFile == ParsingStatus.NONE ) {

                }
                else if ( splitFile == ParsingStatus.SPLIT_LINES ) {
                    // todo -- look at header and determine number of elements that need to be parsed.  Should be static per file
                    int nParts = ParsingUtils.split(line, parts, '\t');
                } else {
                    vc = codec.decode(line);
                    if ( USE_VCF3 ) {
                        VCFRecord rec = (VCFRecord)vc;
                        GenomeLoc loc = GenomeLocParser.createGenomeLoc(rec.getChr(), rec.getStart());
                        ReferenceContext ref = new ReferenceContext(loc, (byte)rec.getReference().charAt(0));
                        vc = VariantContextAdaptors.toVariantContext("X", vc, ref);
                    }
                }
                
                if ( lineNumber % 10000 == 0 ) {
                    System.out.printf("%10d: %s%n", lineNumber, line.subSequence(0, 50));
                    System.out.printf("%10d: %s%n", lineNumber, vc);
                }
            }
        } catch ( FileNotFoundException e ) {
            throw new StingException(e.getMessage());
        } catch ( IOException e ) {
            throw new StingException(e.getMessage());
        }
    }
}