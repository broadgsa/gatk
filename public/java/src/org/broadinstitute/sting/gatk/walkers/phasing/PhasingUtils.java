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

package org.broadinstitute.sting.gatk.walkers.phasing;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * [Short one sentence description of this walker]
 * <p/>
 * <p>
 * [Functionality of this walker]
 * </p>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 *  </pre>
 *
 * @author Your Name
 * @since Date created
 */
class PhasingUtils {
    static VariantContext mergeIntoMNP(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile, AlleleMergeRule alleleMergeRule) {
        if (!mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2))
            return null;

        // Check that it's logically possible to merge the VCs:
        if (!allSamplesAreMergeable(vc1, vc2))
            return null;

        // Check if there's a "point" in merging the VCs (e.g., annotations could be changed)
        if (!alleleMergeRule.allelesShouldBeMerged(vc1, vc2))
            return null;

        return reallyMergeIntoMNP(vc1, vc2, referenceFile);
    }

    static VariantContext reallyMergeIntoMNP(VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile) {
        int startInter = vc1.getEnd() + 1;
        int endInter = vc2.getStart() - 1;
        byte[] intermediateBases = null;
        if (startInter <= endInter) {
            intermediateBases = referenceFile.getSubsequenceAt(vc1.getChr(), startInter, endInter).getBases();
            StringUtil.toUpperCase(intermediateBases);
        }
        MergedAllelesData mergeData = new MergedAllelesData(intermediateBases, vc1, vc2); // ensures that the reference allele is added

        GenotypesContext mergedGenotypes = GenotypesContext.create();
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            List<Allele> mergedAllelesForSample = new LinkedList<Allele>();

            /* NOTE: Since merged alleles are added to mergedAllelesForSample in the SAME order as in the input VC records,
               we preserve phase information (if any) relative to whatever precedes vc1:
             */
            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next(); // this is OK, since allSamplesAreMergeable()

                Allele mergedAllele = mergeData.ensureMergedAllele(all1, all2);
                mergedAllelesForSample.add(mergedAllele);
            }

            double mergedGQ = Math.max(gt1.getLog10PError(), gt2.getLog10PError());

            Map<String, Object> mergedGtAttribs = new HashMap<String, Object>();
            PhaseAndQuality phaseQual = calcPhaseForMergedGenotypes(gt1, gt2);
            if (phaseQual.PQ != null)
                mergedGtAttribs.put(ReadBackedPhasingWalker.PQ_KEY, phaseQual.PQ);

            Genotype mergedGt = new GenotypeBuilder(gt1.getSampleName(), mergedAllelesForSample).log10PError(mergedGQ).attributes(mergedGtAttribs).phased(phaseQual.isPhased).make();
            mergedGenotypes.add(mergedGt);
        }

        String mergedName = mergeVariantContextNames(vc1.getSource(), vc2.getSource());
        double mergedLog10PError = Math.min(vc1.getLog10PError(), vc2.getLog10PError());
        Set<String> mergedFilters = new HashSet<String>(); // Since vc1 and vc2 were unfiltered, the merged record remains unfiltered
        Map<String, Object> mergedAttribs = mergeVariantContextAttributes(vc1, vc2);

        // ids
        List<String> mergedIDs = new ArrayList<String>();
        if ( vc1.hasID() ) mergedIDs.add(vc1.getID());
        if ( vc2.hasID() ) mergedIDs.add(vc2.getID());
        String mergedID = mergedIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(VCFConstants.ID_FIELD_SEPARATOR, mergedIDs);

        VariantContextBuilder mergedBuilder = new VariantContextBuilder(mergedName, vc1.getChr(), vc1.getStart(), vc2.getEnd(), mergeData.getAllMergedAlleles()).id(mergedID).genotypes(mergedGenotypes).log10PError(mergedLog10PError).filters(mergedFilters).attributes(mergedAttribs);
        VariantContextUtils.calculateChromosomeCounts(mergedBuilder, true);
        return mergedBuilder.make();
    }

    static String mergeVariantContextNames(String name1, String name2) {
        return name1 + "_" + name2;
    }

    static Map<String, Object> mergeVariantContextAttributes(VariantContext vc1, VariantContext vc2) {
        Map<String, Object> mergedAttribs = new HashMap<String, Object>();

        List<VariantContext> vcList = new LinkedList<VariantContext>();
        vcList.add(vc1);
        vcList.add(vc2);

        String[] MERGE_OR_ATTRIBS = {VCFConstants.DBSNP_KEY};
        for (String orAttrib : MERGE_OR_ATTRIBS) {
            boolean attribVal = false;
            for (VariantContext vc : vcList) {
                attribVal = vc.getAttributeAsBoolean(orAttrib, false);
                if (attribVal) // already true, so no reason to continue:
                    break;
            }
            mergedAttribs.put(orAttrib, attribVal);
        }

        return mergedAttribs;
    }

    static boolean mergeIntoMNPvalidationCheck(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = VariantContextUtils.getLocation(genomeLocParser, vc1);
        GenomeLoc loc2 = VariantContextUtils.getLocation(genomeLocParser, vc2);

        if (!loc1.onSameContig(loc2))
            throw new ReviewedStingException("Can only merge vc1, vc2 if on the same chromosome");

        if (!loc1.isBefore(loc2))
            throw new ReviewedStingException("Can only merge if vc1 is BEFORE vc2");

        if (vc1.isFiltered() || vc2.isFiltered())
            return false;

        if (!vc1.getSampleNames().equals(vc2.getSampleNames())) // vc1, vc2 refer to different sample sets
            return false;

        if (!allGenotypesAreUnfilteredAndCalled(vc1) || !allGenotypesAreUnfilteredAndCalled(vc2))
            return false;

        return true;
    }

    static boolean allGenotypesAreUnfilteredAndCalled(VariantContext vc) {
        for (final Genotype gt : vc.getGenotypes()) {
            if (gt.isNoCall() || gt.isFiltered())
                return false;
        }

        return true;
    }

    static boolean allSamplesAreMergeable(VariantContext vc1, VariantContext vc2) {
        // Check that each sample's genotype in vc2 is uniquely appendable onto its genotype in vc1:
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            if (!alleleSegregationIsKnown(gt1, gt2)) // can merge if: phased, or if either is a hom
                return false;
        }

        return true;
    }

    static boolean alleleSegregationIsKnown(Genotype gt1, Genotype gt2) {
        if (gt1.getPloidy() != gt2.getPloidy())
            return false;

        /* If gt2 is phased or hom, then could even be MERGED with gt1 [This is standard].

           HOWEVER, EVEN if this is not the case, but gt1.isHom(),
           it is trivially known that each of gt2's alleles segregate with the single allele type present in gt1.
         */
        return (gt2.isPhased() || gt2.isHom() || gt1.isHom());
    }

    static PhaseAndQuality calcPhaseForMergedGenotypes(Genotype gt1, Genotype gt2) {
        if (gt2.isPhased() || gt2.isHom())
            return new PhaseAndQuality(gt1); // maintain the phase of gt1

        if (!gt1.isHom())
            throw new ReviewedStingException("alleleSegregationIsKnown(gt1, gt2) implies: gt2.genotypesArePhased() || gt2.isHom() || gt1.isHom()");

        /* We're dealing with: gt1.isHom(), gt2.isHet(), !gt2.genotypesArePhased(); so, the merged (het) Genotype is not phased relative to the previous Genotype

           For example, if we're merging the third Genotype with the second one:
           0/1
           1|1
           0/1

           Then, we want to output:
           0/1
           1/2
         */
        return new PhaseAndQuality(gt2); // maintain the phase of gt2 [since !gt2.genotypesArePhased()]
    }

    static boolean someSampleHasDoubleNonReferenceAllele(VariantContext vc1, VariantContext vc2) {
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next(); // this is OK, since allSamplesAreMergeable()

                if (all1.isNonReference() && all2.isNonReference()) // corresponding alleles are alternate
                    return true;
            }
        }

        return false;
    }

    static boolean doubleAllelesSegregatePerfectlyAmongSamples(VariantContext vc1, VariantContext vc2) {
        // Check that Alleles at vc1 and at vc2 always segregate together in all samples (including reference):
        Map<Allele, Allele> allele1ToAllele2 = new HashMap<Allele, Allele>();
        Map<Allele, Allele> allele2ToAllele1 = new HashMap<Allele, Allele>();

        // Note the segregation of the alleles for the reference genome:
        allele1ToAllele2.put(vc1.getReference(), vc2.getReference());
        allele2ToAllele1.put(vc2.getReference(), vc1.getReference());

        // Note the segregation of the alleles for each sample (and check that it is consistent with the reference and all previous samples).
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next();

                Allele all1To2 = allele1ToAllele2.get(all1);
                if (all1To2 == null)
                    allele1ToAllele2.put(all1, all2);
                else if (!all1To2.equals(all2)) // all1 segregates with two different alleles at site 2
                    return false;

                Allele all2To1 = allele2ToAllele1.get(all2);
                if (all2To1 == null)
                    allele2ToAllele1.put(all2, all1);
                else if (!all2To1.equals(all1)) // all2 segregates with two different alleles at site 1
                    return false;
            }
        }

        return true;
    }

    abstract static class AlleleMergeRule {
        // vc1, vc2 are ONLY passed to allelesShouldBeMerged() if mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2) AND allSamplesAreMergeable(vc1, vc2):
        abstract public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2);

        public String toString() {
            return "all samples are mergeable";
        }
    }

    static class AlleleOneAndTwo {
        private Allele all1;
        private Allele all2;

        public AlleleOneAndTwo(Allele all1, Allele all2) {
            this.all1 = all1;
            this.all2 = all2;
        }

        public int hashCode() {
            return all1.hashCode() + all2.hashCode();
        }

        public boolean equals(Object other) {
            if (!(other instanceof AlleleOneAndTwo))
                return false;

            AlleleOneAndTwo otherAot = (AlleleOneAndTwo) other;
            return (this.all1.equals(otherAot.all1) && this.all2.equals(otherAot.all2));
        }
    }

    static class MergedAllelesData {
        private Map<AlleleOneAndTwo, Allele> mergedAlleles;
        private byte[] intermediateBases;
        private int intermediateLength;

        public MergedAllelesData(byte[] intermediateBases, VariantContext vc1, VariantContext vc2) {
            this.mergedAlleles = new HashMap<AlleleOneAndTwo, Allele>(); // implemented equals() and hashCode() for AlleleOneAndTwo
            this.intermediateBases = intermediateBases;
            this.intermediateLength = this.intermediateBases != null ? this.intermediateBases.length : 0;

            this.ensureMergedAllele(vc1.getReference(), vc2.getReference(), true);
        }

        public Allele ensureMergedAllele(Allele all1, Allele all2) {
            return ensureMergedAllele(all1, all2, false); // false <-> since even if all1+all2 = reference, it was already created in the constructor
        }

        private Allele ensureMergedAllele(Allele all1, Allele all2, boolean creatingReferenceForFirstTime) {
            AlleleOneAndTwo all12 = new AlleleOneAndTwo(all1, all2);
            Allele mergedAllele = mergedAlleles.get(all12);

            if (mergedAllele == null) {
                byte[] bases1 = all1.getBases();
                byte[] bases2 = all2.getBases();

                byte[] mergedBases = new byte[bases1.length + intermediateLength + bases2.length];
                System.arraycopy(bases1, 0, mergedBases, 0, bases1.length);
                if (intermediateBases != null)
                    System.arraycopy(intermediateBases, 0, mergedBases, bases1.length, intermediateLength);
                System.arraycopy(bases2, 0, mergedBases, bases1.length + intermediateLength, bases2.length);

                mergedAllele = Allele.create(mergedBases, creatingReferenceForFirstTime);
                mergedAlleles.put(all12, mergedAllele);
            }

            return mergedAllele;
        }

        public Set<Allele> getAllMergedAlleles() {
            return new HashSet<Allele>(mergedAlleles.values());
        }
    }

    static class PhaseAndQuality {
        public boolean isPhased;
        public Double PQ = null;

        public PhaseAndQuality(Genotype gt) {
            this.isPhased = gt.isPhased();
            if (this.isPhased) {
                this.PQ = gt.getAttributeAsDouble(ReadBackedPhasingWalker.PQ_KEY, -1);
                if ( this.PQ == -1 ) this.PQ = null;
            }
        }
    }
}
