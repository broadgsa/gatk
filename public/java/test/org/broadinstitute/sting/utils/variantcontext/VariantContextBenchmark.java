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
import net.sf.picard.reference.ReferenceSequenceFile;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;

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

    @Param({"10"})
    int dupsToMerge; // set automatically by framework

    @Param
    Operation operation; // set automatically by framework

    private String INPUT_STRING;

    public enum Operation {
        READ,
        SUBSET_TO_SAMPLES,
        GET_TYPE,
        GET_ID,
        GET_GENOTYPES,
        GET_ATTRIBUTE_STRING,
        GET_ATTRIBUTE_INT,
        GET_N_SAMPLES,
        GET_GENOTYPES_FOR_SAMPLES,
        GET_GENOTYPES_IN_ORDER_OF_NAME,
        CALC_GENOTYPE_COUNTS,
        MERGE
    }

    private GenomeLocParser b37GenomeLocParser;

    @Override protected void setUp() {
        // TODO -- update for new tribble interface
//        try {
//            ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(BaseTest.b37KGReference));
//            b37GenomeLocParser = new GenomeLocParser(seq);
//        } catch ( FileNotFoundException e) {
//            throw new RuntimeException(e);
//        }
//
//        // read it into a String so that we don't try to benchmark IO issues
//        try {
//            FileInputStream s = new FileInputStream(new File(vcfFile));
//            AsciiLineReader lineReader = new AsciiLineReader(s);
//            int counter = 0;
//            StringBuffer sb = new StringBuffer();
//            while (counter++ < linesToRead ) {
//                String line = lineReader.readLine();
//                if ( line == null )
//                    break;
//                sb.append(line + "\n");
//            }
//            s.close();
//            INPUT_STRING = sb.toString();
//        } catch (IOException e) {
//            throw new RuntimeException(e);
//        }
    }

    private interface FunctionToBenchmark<T extends Feature> {
        public void run(T vc);
    }

    private <T extends Feature> void runBenchmark(FeatureCodec<T> codec, FunctionToBenchmark<T> func) {
        // TODO -- update for new Tribble interface
//        try {
//            InputStream is = new ByteArrayInputStream(INPUT_STRING.getBytes());
//            AsciiLineReader lineReader = new AsciiLineReader(is);
//            codec.readHeader(lineReader);
//
//            int counter = 0;
//            while (counter++ < linesToRead ) {
//                String line = lineReader.readLine();
//                if ( line == null )
//                    break;
//
//                T vc = codec.decode(line);
//                func.run(vc);
//            }
//        } catch (Exception e) {
//            System.out.println("Benchmarking run failure because of " + e.getMessage());
//        }
    }

    public void timeV14(int rep) {
        for ( int i = 0; i < rep; i++ ) {
            FunctionToBenchmark<VariantContext> func = getV14FunctionToBenchmark();
            FeatureCodec<VariantContext> codec = new VCFCodec();
            runBenchmark(codec, func);
        }
    }

    public FunctionToBenchmark<VariantContext> getV14FunctionToBenchmark() {
        switch ( operation ) {
            case READ:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        ; // empty operation
                    }
                };
            case SUBSET_TO_SAMPLES:
                return new FunctionToBenchmark<VariantContext>() {
                    Set<String> samples;
                    public void run(final VariantContext vc) {
                        if ( samples == null )
                            samples = new HashSet<String>(new ArrayList<String>(vc.getSampleNames()).subList(0, nSamplesToTake));
                        VariantContext sub = vc.subContextFromSamples(samples, true);
                        sub.getNSamples();
                    }
                };
            case GET_TYPE:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        vc.getType();
                    }
                };
            case GET_ID:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        vc.getID();
                    }
                };
            case GET_GENOTYPES:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        vc.getGenotypes().size();
                    }
                };

            case GET_GENOTYPES_FOR_SAMPLES:
                return new FunctionToBenchmark<VariantContext>() {
                    Set<String> samples;
                    public void run(final VariantContext vc) {
                        if ( samples == null )
                            samples = new HashSet<String>(new ArrayList<String>(vc.getSampleNames()).subList(0, nSamplesToTake));
                        vc.getGenotypes(samples).size();
                    }
                };

            case GET_ATTRIBUTE_STRING:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        vc.getAttribute("AN", null);
                    }
                };

            case GET_ATTRIBUTE_INT:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        vc.getAttributeAsInt("AC", 0);
                    }
                };

            case GET_N_SAMPLES:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        vc.getNSamples();
                    }
                };

            case GET_GENOTYPES_IN_ORDER_OF_NAME:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        ; // TODO - TEST IS BROKEN
//                        int n = 0;
//                        for ( final Genotype g: vc.getGenotypesOrderedByName() ) n++;
                    }
                };

            case CALC_GENOTYPE_COUNTS:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        vc.getHetCount();
                    }
                };

            case MERGE:
                return new FunctionToBenchmark<VariantContext>() {
                    public void run(final VariantContext vc) {
                        List<VariantContext> toMerge = new ArrayList<VariantContext>();

                        for ( int i = 0; i < dupsToMerge; i++ ) {
                            GenotypesContext gc = GenotypesContext.create(vc.getNSamples());
                            for ( final Genotype g : vc.getGenotypes() ) {
                                gc.add(new GenotypeBuilder(g).name(g.getSampleName()+"_"+i).make());
                            }
                            toMerge.add(new VariantContextBuilder(vc).genotypes(gc).make());
                        }

                        VariantContextUtils.simpleMerge(b37GenomeLocParser, toMerge, null,
                                VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                                VariantContextUtils.GenotypeMergeType.UNSORTED,
                                true, false, "set", false, true);
                    }
                };

            default: throw new IllegalArgumentException("Unexpected operation " + operation);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // V13
    //
    // In order to use this, you must move the v13 version from archive and uncomment
    //
    // git mv private/archive/java/src/org/broadinstitute/sting/utils/variantcontext/v13 public/java/test/org/broadinstitute/sting/utils/variantcontext/v13
    //
    // --------------------------------------------------------------------------------

//    public void timeV13(int rep) {
//        for ( int i = 0; i < rep; i++ ) {
//            FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext> func = getV13FunctionToBenchmark();
//            FeatureCodec<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext> codec = new org.broadinstitute.sting.utils.variantcontext.v13.VCFCodec();
//            runBenchmark(codec, func);
//        }
//    }
//
//    public FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext> getV13FunctionToBenchmark() {
//        switch ( operation ) {
//            case READ:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        ; // empty operation
//                    }
//                };
//            case SUBSET_TO_SAMPLES:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    List<String> samples;
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        if ( samples == null )
//                            samples = new ArrayList<String>(vc.getSampleNames()).subList(0, nSamplesToTake);
//                        org.broadinstitute.sting.utils.variantcontext.v13.VariantContext sub = vc.subContextFromGenotypes(vc.getGenotypes(samples).values());
//                        sub.getNSamples();
//                    }
//                };
//
//            case GET_TYPE:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        vc.getType();
//                    }
//                };
//            case GET_ID:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        vc.getID();
//                    }
//                };
//            case GET_GENOTYPES:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        vc.getGenotypes().size();
//                    }
//                };
//
//            case GET_GENOTYPES_FOR_SAMPLES:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    Set<String> samples;
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        if ( samples == null )
//                            samples = new HashSet<String>(new ArrayList<String>(vc.getSampleNames()).subList(0, nSamplesToTake));
//                        vc.getGenotypes(samples).size();
//                    }
//                };
//
//            case GET_ATTRIBUTE_STRING:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        vc.getAttribute("AN", null);
//                    }
//                };
//
//            case GET_ATTRIBUTE_INT:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        vc.getAttributeAsInt("AC", 0);
//                    }
//                };
//
//            case GET_N_SAMPLES:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        vc.getNSamples();
//                    }
//                };
//
//            case GET_GENOTYPES_IN_ORDER_OF_NAME:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        ; // TODO - TEST IS BROKEN
//                        //vc.getGenotypesOrderedByName();
//                    }
//                };
//
//            case CALC_GENOTYPE_COUNTS:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        vc.getHetCount();
//                    }
//                };
//
//            case MERGE:
//                return new FunctionToBenchmark<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>() {
//                    public void run(final org.broadinstitute.sting.utils.variantcontext.v13.VariantContext vc) {
//                        List<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext> toMerge = new ArrayList<org.broadinstitute.sting.utils.variantcontext.v13.VariantContext>();
//
//                        for ( int i = 0; i < dupsToMerge; i++ ) {
//                            Map<String, org.broadinstitute.sting.utils.variantcontext.v13.Genotype> gc = new HashMap<String, org.broadinstitute.sting.utils.variantcontext.v13.Genotype>();
//                            for ( final org.broadinstitute.sting.utils.variantcontext.v13.Genotype g : vc.getGenotypes().values() ) {
//                                String name = g.getSampleName()+"_"+i;
//                                gc.put(name, new org.broadinstitute.sting.utils.variantcontext.v13.Genotype(name,
//                                        g.getAlleles(), g.getLog10PError(), g.getFilters(), g.getAttributes(), g.isPhased(), g.getLikelihoods().getAsVector()));
//                                toMerge.add(org.broadinstitute.sting.utils.variantcontext.v13.VariantContext.modifyGenotypes(vc, gc));
//                            }
//                        }
//
//                        org.broadinstitute.sting.utils.variantcontext.v13.VariantContextUtils.simpleMerge(b37GenomeLocParser,
//                                toMerge, null,
//                                org.broadinstitute.sting.utils.variantcontext.v13.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
//                                org.broadinstitute.sting.utils.variantcontext.v13.VariantContextUtils.GenotypeMergeType.UNSORTED,
//                                true, false, "set", false, true);
//                    }
//                };
//
//            default: throw new IllegalArgumentException("Unexpected operation " + operation);
//        }
//    }

    public static void main(String[] args) {
        CaliperMain.main(VariantContextBenchmark.class, args);
    }
}
