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
package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.phasing.WriteVCF;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

/**
 * At each locus in the input data set, prints the reference in VCF format.
 */

@By(DataSource.REFERENCE)
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES})

public class PrintReferenceVariantsWalker extends LocusWalker<Integer, Integer> {
    @Output(doc = "File to which reference variants should be written", required = true)
    protected VCFWriter writer = null;

    private final static String REFERENCE = "REFERENCE";
    private final static double REF_NEG_LOG_10_P_ERROR = 9.9;

    private String REF_FILE_NAME = null;

    public void initialize() {
        this.REF_FILE_NAME = getToolkit().getArguments().referenceFile.getName();

        initializeVcfWriter();
    }

    private void initializeVcfWriter() {
        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", REF_FILE_NAME));

        Set<String> samples = new TreeSet<String>();
        samples.add(REFERENCE);
        writer.writeHeader(new VCFHeader(hInfo, samples));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (ref == null)
            return 0;

        GenomeLoc refLoc = ref.getLocus();

        Allele refAllele = Allele.create(ref.getBase(), true); // create a single-base allele
        Set<Allele> alleles = new HashSet<Allele>();
        alleles.add(refAllele);

        Map<String, Genotype> genotypes = new HashMap<String, Genotype>();
        boolean isPhased = true; // trivially true for a haploid genotype
        Genotype haploidRefGt = new Genotype(REFERENCE, new LinkedList<Allele>(alleles), VCFConstants.MAX_GENOTYPE_QUAL, new HashSet<String>(), new HashMap<String, Object>(), isPhased);
        genotypes.put(REFERENCE, haploidRefGt);

        // Ensure that the genotype refers to alleles of length 1 (by using refLoc.getStart() as the stop position):
        VariantContext vc = new VariantContext(REF_FILE_NAME, refLoc.getContig(), refLoc.getStart(), refLoc.getStart(), alleles, genotypes, REF_NEG_LOG_10_P_ERROR, new HashSet<String>(), new HashMap<String, Object>());

        WriteVCF.writeVCF(vc, writer, logger);

        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer addIn, Integer sum) {
        return sum + addIn;
    }

    public void onTraversalDone(Integer result) {
        System.out.println("Processed " + result + " sites.");
    }
}