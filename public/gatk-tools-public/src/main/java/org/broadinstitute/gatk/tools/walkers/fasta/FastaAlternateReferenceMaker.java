/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.List;
import java.util.Set;


/**
 * Generates an alternative reference sequence over the specified interval.
 *
 * <p>
 * Given variant tracks, it replaces the reference bases at variation sites with the bases supplied by the ROD(s).
 * Additionally, allows for one or more "snpmask" VCFs to set overlapping bases to 'N'.
 *
 * The output format can be partially controlled using the provided command-line arguments.
 * Specify intervals with the usual -L argument to output only the reference bases within your intervals.
 * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
 * separate fasta sequence (named numerically in order).
 *
 * Several important notes:
 * 1) if there are multiple variants that start at a site, it chooses one of them randomly.
 * 2) when there are overlapping indels (but with different start positions) only the first will be chosen.
 * 3) this tool works only for SNPs and for simple indels (but not for things like complex substitutions).
 * Reference bases for each interval will be output as a separate fasta sequence (named numerically in order).
 *
 * <h3>Input</h3>
 * <p>
 * The reference, requested intervals, and any number of variant rod files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A fasta file representing the requested intervals.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T FastaAlternateReferenceMaker \
 *   -o output.fasta \
 *   -L input.intervals \
 *   --variant input.vcf \
 *   [--snpmask mask.vcf]
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_REFUTILS, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-1,stop=50))
@Requires(value={DataSource.REFERENCE})
public class FastaAlternateReferenceMaker extends FastaReferenceMaker {

    /**
     * Variants from this input file are used by this tool to construct an alternate reference.
     */
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * Snps from this file are used as a mask (inserting N's in the sequence) when constructing the alternate reference
     * (regardless of whether they overlap a variant site).
     */
    @Input(fullName="snpmask", shortName = "snpmask", doc="SNP mask VCF file", required=false)
    protected RodBinding<VariantContext> snpmask;

    /**
     * This option will generate an error if the specified sample does not exist in the VCF.
     * Non-diploid (or non-called) genotypes are ignored.
     */
    @Argument(fullName="use_IUPAC_sample", shortName="IUPAC", doc = "If specified, heterozygous SNP sites will be output using IUPAC ambiguity codes given the genotypes for this sample", required=false)
    private String iupacSample = null;

    private int deletionBasesRemaining = 0;

    @Override
    public void initialize() {
        super.initialize();
        if ( iupacSample != null ) {
            final List<String> rodName = Arrays.asList(variantCollection.variants.getName());
            final Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName);
            if ( !samples.contains(iupacSample) )
                throw new UserException.BadInput("the IUPAC sample specified is not present in the provided VCF file");
        }
    }

    @Override
    public Pair<GenomeLoc, String> map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {

        if (deletionBasesRemaining > 0) {
            deletionBasesRemaining--;
            return new Pair<>(context.getLocation(), "");
        }

        final String refBase = String.valueOf((char)ref.getBase());

        // Check to see if we have a called snp
        for ( final VariantContext vc : tracker.getValues(variantCollection.variants, ref.getLocus()) ) {
            if ( vc.isFiltered() )
                continue;

            if ( vc.isSimpleDeletion()) {
                deletionBasesRemaining = vc.getReference().length() - 1;
                // delete the next n bases, not this one
                return new Pair<>(context.getLocation(), refBase);
            } else if ( vc.isSimpleInsertion()) {
                return new Pair<>(context.getLocation(), vc.getAlternateAllele(0).toString());
            } else if (vc.isSNP()) {
                final String base = (iupacSample != null) ? getIUPACbase(vc.getGenotype(iupacSample), refBase) : vc.getAlternateAllele(0).toString();
                return new Pair<>(context.getLocation(), base);
            }
        }

        // if we don't have a called site, and we have a mask at this site, mask it
        for ( final VariantContext vc : tracker.getValues(snpmask) ) {
            if ( vc.isSNP()) {
                return new Pair<>(context.getLocation(), "N");
            }
        }

        // if we got here then we're just ref
        return new Pair<>(context.getLocation(), refBase);
    }

    /**
     * Returns the IUPAC encoding for the given genotype or the reference base if not possible
     *
     * @param genotype  the genotype to encode
     * @param ref       the reference base
     * @return non-null, non-empty String
     */
    private String getIUPACbase(final Genotype genotype, final String ref) {
        if ( genotype == null )
            throw new IllegalStateException("The genotype is null for sample " + iupacSample);

        if ( !genotype.isHet() )
            return genotype.isHom() ? genotype.getAllele(0).getBaseString() : ref;

        final byte allele1 = genotype.getAllele(0).getBases()[0];
        final byte allele2 = genotype.getAllele(1).getBases()[0];
        return new String(new byte[] {BaseUtils.basesToIUPAC(allele1, allele2)});
    }
}