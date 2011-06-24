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

package org.broadinstitute.sting.playground.gatk.walkers.variantutils;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;

import java.util.*;

/**
 * Converts variants from the Complete Genomics VAR format to VCF format.
 */
@Requires(value={},referenceMetaData=@RMD(name=CGVarToVCF.INPUT_ROD_NAME, type=VariantContext.class))
@Reference(window=@Window(start=-40,stop=400))
public class CGVarToVCF extends RodWalker<Integer, Integer> {

    public static final String INPUT_ROD_NAME = "variant";

    @Output(doc="File to which variants should be written", required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="sample", shortName="sample", doc="The sample name represented by the variant rod", required=true)
    protected String sampleName = null;

    public void initialize() {
        HashSet<String> samples = new HashSet<String>(1);
        samples.add(sampleName);
        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), samples));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> contexts = tracker.getVariantContexts(ref, INPUT_ROD_NAME, null, ref.getLocus(), true, false);

        // for now, we don't support the mixed type
        if ( contexts.size() == 0 || contexts.size() > 2 )
            return 0;

        Iterator<VariantContext> iter = contexts.iterator();
        if ( contexts.size() == 1 ) {
            writeHet(iter.next(), ref.getBase());
        } else {
            VariantContext vc1 = iter.next();
            VariantContext vc2 = iter.next();
            if ( vc1.getType().equals(vc2.getType()) )
                writeHom(vc1, ref.getBase());
        }

        return 0;
    }

    private void writeHet(VariantContext vc, byte ref) {
        List<Allele> alleles = new ArrayList<Allele>(vc.getAlleles());
        Genotype g = new Genotype(sampleName, alleles);
        write(vc, ref, g);
    }

    private void writeHom(VariantContext vc, byte ref) {
        List<Allele> alleles = new ArrayList<Allele>(2);
        alleles.add(vc.getAlternateAllele(0));
        alleles.add(vc.getAlternateAllele(0));
        Genotype g = new Genotype(sampleName, alleles);
        write(vc, ref, g);
    }

    private void write(VariantContext vc, byte ref, Genotype g) {
        HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>(1);
        genotypes.put(sampleName, g);
        vc = VariantContext.modifyGenotypes(vc, genotypes);
        if ( vc.isSNP() )
            vc = VariantContext.modifyLocation(vc, vc.getChr(), vc.getStart()+1, vc.getStart()+1);        
        vcfWriter.add(vc, ref);
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
