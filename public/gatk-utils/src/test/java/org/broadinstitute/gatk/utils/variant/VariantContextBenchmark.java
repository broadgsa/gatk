/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.variant;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFCodec;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

    private <T extends Feature> void runBenchmark(FeatureCodec codec, FunctionToBenchmark<T> func) {
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
            final VCFCodec codec = new VCFCodec();
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
                            samples = new HashSet<>(new ArrayList<>(vc.getSampleNames()).subList(0, nSamplesToTake));
                        VariantContext sub = vc.subContextFromSamples(samples);
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
                            samples = new HashSet<>(new ArrayList<>(vc.getSampleNames()).subList(0, nSamplesToTake));
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
                        List<VariantContext> toMerge = new ArrayList<>();

                        for ( int i = 0; i < dupsToMerge; i++ ) {
                            GenotypesContext gc = GenotypesContext.create(vc.getNSamples());
                            for ( final Genotype g : vc.getGenotypes() ) {
                                gc.add(new GenotypeBuilder(g).name(g.getSampleName()+"_"+i).make());
                            }
                            toMerge.add(new VariantContextBuilder(vc).genotypes(gc).make());
                        }

                        GATKVariantContextUtils.simpleMerge(toMerge, null,
                                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                                GATKVariantContextUtils.GenotypeMergeType.UNSORTED,
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
//            FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext> func = getV13FunctionToBenchmark();
//            FeatureCodec<htsjdk.variant.variantcontext.v13.VariantContext> codec = new htsjdk.variant.variantcontext.v13.VCFCodec();
//            runBenchmark(codec, func);
//        }
//    }
//
//    public FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext> getV13FunctionToBenchmark() {
//        switch ( operation ) {
//            case READ:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        ; // empty operation
//                    }
//                };
//            case SUBSET_TO_SAMPLES:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    List<String> samples;
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        if ( samples == null )
//                            samples = new ArrayList<String>(vc.getSampleNames()).subList(0, nSamplesToTake);
//                        htsjdk.variant.variantcontext.v13.VariantContext sub = vc.subContextFromGenotypes(vc.getGenotypes(samples).values());
//                        sub.getNSamples();
//                    }
//                };
//
//            case GET_TYPE:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        vc.getType();
//                    }
//                };
//            case GET_ID:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        vc.getID();
//                    }
//                };
//            case GET_GENOTYPES:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        vc.getGenotypes().size();
//                    }
//                };
//
//            case GET_GENOTYPES_FOR_SAMPLES:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    Set<String> samples;
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        if ( samples == null )
//                            samples = new HashSet<String>(new ArrayList<String>(vc.getSampleNames()).subList(0, nSamplesToTake));
//                        vc.getGenotypes(samples).size();
//                    }
//                };
//
//            case GET_ATTRIBUTE_STRING:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        vc.getExtendedAttribute("AN", null);
//                    }
//                };
//
//            case GET_ATTRIBUTE_INT:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        vc.getAttributeAsInt("AC", 0);
//                    }
//                };
//
//            case GET_N_SAMPLES:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        vc.getNSamples();
//                    }
//                };
//
//            case GET_GENOTYPES_IN_ORDER_OF_NAME:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        ; // TODO - TEST IS BROKEN
//                        //vc.getGenotypesOrderedByName();
//                    }
//                };
//
//            case CALC_GENOTYPE_COUNTS:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        vc.getHetCount();
//                    }
//                };
//
//            case MERGE:
//                return new FunctionToBenchmark<htsjdk.variant.variantcontext.v13.VariantContext>() {
//                    public void run(final htsjdk.variant.variantcontext.v13.VariantContext vc) {
//                        List<htsjdk.variant.variantcontext.v13.VariantContext> toMerge = new ArrayList<htsjdk.variant.variantcontext.v13.VariantContext>();
//
//                        for ( int i = 0; i < dupsToMerge; i++ ) {
//                            Map<String, htsjdk.variant.variantcontext.v13.Genotype> gc = new HashMap<String, htsjdk.variant.variantcontext.v13.Genotype>();
//                            for ( final htsjdk.variant.variantcontext.v13.Genotype g : vc.getGenotypes().values() ) {
//                                String name = g.getSampleName()+"_"+i;
//                                gc.put(name, new htsjdk.variant.variantcontext.v13.Genotype(name,
//                                        g.getAlleles(), g.getLog10PError(), g.getFilters(), g.getAttributes(), g.isPhased(), g.getLikelihoods().getAsVector()));
//                                toMerge.add(htsjdk.variant.variantcontext.v13.VariantContext.modifyGenotypes(vc, gc));
//                            }
//                        }
//
//                        htsjdk.variant.variantcontext.v13.VariantContextUtils.simpleMerge(b37GenomeLocParser,
//                                toMerge, null,
//                                htsjdk.variant.variantcontext.v13.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
//                                htsjdk.variant.variantcontext.v13.VariantContextUtils.GenotypeMergeType.UNSORTED,
//                                true, false, "set", false, true, false);
//                    }
//                };
//
//            default: throw new IllegalArgumentException("Unexpected operation " + operation);
//        }
//    }

    public static void main(String[] args) {
        com.google.caliper.Runner.main(VariantContextBenchmark.class, args);
    }
}
