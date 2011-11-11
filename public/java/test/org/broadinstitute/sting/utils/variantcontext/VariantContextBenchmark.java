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

package org.broadinstitute.sting.utils.variantcontext;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import com.google.caliper.runner.CaliperMain;
import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;

import java.io.*;
import java.util.*;

/**
 * Caliper microbenchmark of parsing a VCF file
 */
public class VariantContextBenchmark extends SimpleBenchmark {
    @Param({"/Users/depristo/Desktop/broadLocal/localData/ALL.chr20.merged_beagle_mach.20101123.snps_indels_svs.genotypes.vcf"})
    String vcfFile;

    @Param({"1000"})
    int linesToRead; // set automatically by framework

    @Param({"100"})
    int nSamplesToTake; // set automatically by framework

    @Param({"READ", "READ_SUBSET"})
    Operation operation; // set automatically by framework

    @Param({"OF_GENOTYPES", "OF_SAMPLES"})
    SubContextOp subContextOp; // set automatically by framework

    private String INPUT_STRING;

    public enum Operation {
        READ,
        READ_SUBSET
    }

    public enum SubContextOp {
        OF_GENOTYPES,
        OF_SAMPLES
    }

    @Override protected void setUp() {
        // read it into a String so that we don't try to benchmark IO issues
        try {
            FileInputStream s = new FileInputStream(new File(vcfFile));
            AsciiLineReader lineReader = new AsciiLineReader(s);
            int counter = 0;
            StringBuffer sb = new StringBuffer();
            while (counter++ < linesToRead ) {
                String line = lineReader.readLine();
                if ( line == null )
                    break;
                sb.append(line + "\n");
            }
            s.close();
            INPUT_STRING = sb.toString();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void parseGenotypes(VCFCodec codec, Operation op, SubContextOp subop ) {
        try {
            InputStream is = new ByteArrayInputStream(INPUT_STRING.getBytes());
            AsciiLineReader lineReader = new AsciiLineReader(is);
            codec.readHeader(lineReader);

            int counter = 0;
            List<String> samples = null;
            while (counter++ < linesToRead ) {
                String line = lineReader.readLine();
                if ( line == null )
                    break;

                VariantContext vc = (VariantContext)codec.decode(line);
                if ( samples == null ) {
                    samples = new ArrayList<String>(vc.getSampleNames()).subList(0, nSamplesToTake);
                }

                if ( op == Operation.READ_SUBSET)
                    processOneVC(vc, samples, subop);
            }
        } catch (Exception e) {
            System.out.println("Benchmarking run failure because of " + e.getMessage());
        }
    }

    public void timeMe(int rep) {
        for ( int i = 0; i < rep; i++ )
            parseGenotypes(new VCFCodec(), operation, subContextOp);
    }

    public static void main(String[] args) {
        CaliperMain.main(VariantContextBenchmark.class, args);
    }

    private static final void processOneVC(VariantContext vc, List<String> samples, SubContextOp subop) {
        VariantContext sub;

        switch ( subop ) {
            case OF_GENOTYPES:
                sub = vc.subContextFromGenotypes(vc.getGenotypes(samples).values(), vc.getAlleles());
                break;
            case OF_SAMPLES:
                sub = vc.subContextFromSamples(samples, vc.getAlleles());
                break;
            default:
                throw new RuntimeException("Unexpected op: " + subop);
        }

        sub.getNSamples();
    }
}
