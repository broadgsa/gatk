/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.vcf.VCFRecord;
import org.broad.tribble.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;

import java.io.PrintStream;
import java.util.*;

/**
 * Produces an input file to Beagle imputation engine, listing genotype likelihoods for each sample in input VCF file
 */
@Requires(value={},referenceMetaData=@RMD(name="vcf",type= VCFRecord.class))
public class ProduceBeagleInputWalker extends RodWalker<Integer, Integer> {

    @Argument(fullName = "beagle_file", shortName = "beagle", doc = "File to print BEAGLE-specific data for use with imputation", required = true)
    public PrintStream beagleWriter = null;

    final TreeSet<String> samples = new TreeSet<String>();

    public void initialize() {

        final List<ReferenceOrderedDataSource> dataSources = this.getToolkit().getRodDataSources();
        for ( final ReferenceOrderedDataSource source : dataSources ) {
            final RMDTrack rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(VCFRecord.class) ) {
                final VCFReader reader = new VCFReader(rod.getFile());
                final Set<String> vcfSamples = reader.getHeader().getGenotypeSamples();
                samples.addAll(vcfSamples);
                reader.close();
            }
        }

        beagleWriter.print("marker alleleA alleleB");
        for ( String sample : samples )
            beagleWriter.print(String.format(" %s %s %s", sample, sample, sample));

        beagleWriter.println();
    }
    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

         if( tracker != null ) {
            EnumSet<VariantContext.Type> vc = EnumSet.of(VariantContext.Type.SNP);
            GenomeLoc loc = context.getLocation();
            VariantContext vc_eval;


            try {
                vc_eval = tracker.getVariantContext(ref,"vcf", vc, loc, true);
            } catch (java.util.NoSuchElementException e) {
                return 0;
            }

            // output marker ID to Beagle input file
            beagleWriter.print(String.format("%s ",vc_eval.getLocation().toString()));

            for (Allele allele: vc_eval.getAlleles()) {
                // TODO -- check whether this is really needed by Beagle
                String bglPrintString;
                if (allele.isNoCall() || allele.isNull())
                    bglPrintString = "0";
                else
                    bglPrintString = allele.toString().substring(0,1);  // get rid of * in case of reference allele

                beagleWriter.print(String.format("%s ", bglPrintString));
            }

            if ( !vc_eval.hasGenotypes() )
                return null;

            Map<String, Genotype> genotypes = vc_eval.getGenotypes();
            for ( String sample : samples ) {
                // use sample as key into genotypes structure
                Genotype genotype = genotypes.get(sample);
                if (genotype.isCalled() && genotype.hasAttribute(VCFGenotypeRecord.GENOTYPE_LIKELIHOODS_KEY)) {
                    String[] glArray = genotype.getAttributeAsString(VCFGenotypeRecord.GENOTYPE_LIKELIHOODS_KEY).split(",");

                    for (String gl : glArray) {
                        Double d_gl = Math.pow(10, Double.valueOf(gl));
                        beagleWriter.print(String.format("%5.2f ",d_gl));
                    }
                }
                else
                    beagleWriter.print("0 0 0 "); // write 0 likelihoods for uncalled genotypes.

            }

            beagleWriter.println();

        }
        return 1;
        
    }

    public Integer reduceInit() {
        return 0; // Nothing to do here
    }

    public Integer reduce( Integer value, Integer sum ) {
        return 0; // Nothing to do here
    }

    public void onTraversalDone( Integer sum ) {

   }
    
}
