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

package org.broadinstitute.sting.utils.codecs.vcf;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * All of the gross allele clipping and padding routines in one place
 *
 * Having attempted to understand / fix / document this code myself
 * I can only conclude that this entire approach needs to be rethought.  This
 * code just doesn't work robustly with symbolic alleles, with multiple alleles,
 * requires a special "reference base for indels" stored in the VariantContext
 * whose correctness isn't enforced, and overall has strange special cases
 * all over the place.
 *
 * The reason this code is so complex is due to symbolics and multi-alleleic
 * variation, which frequently occur when combining variants from multiple
 * VCF files.
 *
 * TODO rethink this class, make it clean, and make it easy to create, mix, and write out alleles
 *
 * @author Mark DePristo
 * @since 6/12
 */
public final class VCFAlleleClipper {
    // TODO
    // TODO
    // TODO
    // TODO I have honestly no idea what this code is really supposed to do
    // TODO computeForwardClipping is called in two places:
    // TODO
    // TODO clipAlleles() where all alleles, including ref, are passed in
    // TODO createVariantContextWithTrimmedAlleles() where only the alt alleles are passed in
    // TODO
    // TODO The problem is that the code allows us to clip
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    @Ensures("result == 0 || result == 1")
    public static int computeForwardClipping(final List<Allele> unclippedAlleles,
                                             final byte ref0) {
        boolean clipping = true;
        int symbolicAlleleCount = 0;

        for ( final Allele a : unclippedAlleles ) {
            if ( a.isSymbolic() ) {
                symbolicAlleleCount++;
                continue;
            }

            if ( a.length() < 1 || (a.getBases()[0] != ref0) ) {
                clipping = false;
                break;
            }
        }

        // don't clip if all alleles are symbolic
        return (clipping && symbolicAlleleCount != unclippedAlleles.size()) ? 1 : 0;
    }

    public static int computeReverseClipping(final List<Allele> unclippedAlleles,
                                             final byte[] ref,
                                             final int forwardClipping,
                                             final boolean allowFullClip) {
        int clipping = 0;
        boolean stillClipping = true;

        while ( stillClipping ) {
            for ( final Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() )
                    continue;

                // we need to ensure that we don't reverse clip out all of the bases from an allele because we then will have the wrong
                // position set for the VariantContext (although it's okay to forward clip it all out, because the position will be fine).
                if ( a.length() - clipping == 0 )
                    return clipping - (allowFullClip ? 0 : 1);

                if ( a.length() - clipping <= forwardClipping || a.length() - forwardClipping == 0 ) {
                    stillClipping = false;
                }
                else if ( ref.length == clipping ) {
                    if ( allowFullClip )
                        stillClipping = false;
                    else
                        return -1;
                }
                else if ( a.getBases()[a.length()-clipping-1] != ref[ref.length-clipping-1] ) {
                    stillClipping = false;
                }
            }
            if ( stillClipping )
                clipping++;
        }

        return clipping;
    }

    /**
     * Are the alleles describing a polymorphism substitution one base for another?
     *
     * @param alleles a list of alleles, must not be empty
     * @return Return true if the length of any allele in alleles isn't 1
     */
    @Requires("!alleles.isEmpty()")
    private static boolean isSingleNucleotideEvent(final List<Allele> alleles) {
        for ( final Allele a : alleles ) {
            if ( a.length() != 1 )
                return false;
        }
        return true;
    }

    /**
     * clip the alleles, based on the reference, returning a ClippedAlleles object describing what happened
     *
     * The ClippedAlleles object contains the implied stop position of the alleles, given the provided start
     * position, after clipping.  It also contains the list of alleles, in the same order as the provided
     * unclipped ones, that are the fully clipped version of the input alleles.  If an error occurs
     * during this option the getError() function returns a string describing the problem (for use in parsers).
     *
     * The basic operation are:
     *
     * single allele
     *      => stop == start and clipped == unclipped
     * any number of single nucleotide events
     *      => stop == start and clipped == unclipped
     * two alleles, second being symbolic
     *      => stop == start and clipped == unclipped
     *      Note in this case that the STOP should be computed by other means (from END in VCF, for example)
     *      Note that if there's more than two alleles and the second is a symbolic the code produces an error
     * Any other case:
     *      The alleles are trimmed of any sequence shared at the end of the alleles.  If N bases
     *      are common then the alleles will all be at least N bases shorter.
     *      The stop position returned is the start position + the length of the
     *      reverse trimmed only reference allele - 1.
     *      If the alleles all share a single common starting sequence (just one base is considered)
     *      then the alleles have this leading common base removed as well.
     *
     * TODO This code is gross and brittle and needs to be rethought from scratch
     *
     * @param start the unadjusted start position (pre-clipping)
     * @param ref the reference string
     * @param unclippedAlleles the list of unclipped alleles, including the reference allele
     * @return the new reference end position of this event
     */
    @Requires({"position > 0", "ref != null && ref.length() > 0", "!unclippedAlleles.isEmpty()"})
    @Ensures("result != null")
    public static ClippedAlleles clipAlleles(final int start,
                                             final String ref,
                                             final List<Allele> unclippedAlleles) {
        // no variation or single nucleotide events are by definition fully clipped
        if ( unclippedAlleles.size() == 1 || isSingleNucleotideEvent(unclippedAlleles) )
            return new ClippedAlleles(start, unclippedAlleles);

        // symbolic alleles aren't unclipped by default, and there can be only two (better than highlander style)
        if ( unclippedAlleles.size() > 1 ) {
            if ( unclippedAlleles.get(1).isSymbolic() ) {
                if ( unclippedAlleles.size() > 2 )
                    return new ClippedAlleles("Detected multiple alleles at a site, including a symbolic one, which is currently unsupported");
                // we do not unclip symbolic alleles
                return new ClippedAlleles(start, unclippedAlleles);
            }
        }

        // we've got to sort out the clipping by looking at the alleles themselves
        final int forwardClipping = computeForwardClipping(unclippedAlleles, (byte)ref.charAt(0));
        final int reverseClipping = computeReverseClipping(unclippedAlleles, ref.getBytes(), forwardClipping, false);

        if ( reverseClipping == -1 )
            return new ClippedAlleles("computeReverseClipping failed due to bad alleles");

        final List<Allele> clippedAlleles = new ArrayList<Allele>(unclippedAlleles.size());
        for ( final Allele a : unclippedAlleles ) {
            if ( a.isSymbolic() ) {
                clippedAlleles.add(a);
            } else {
                final byte[] allele = Arrays.copyOfRange(a.getBases(), forwardClipping, a.getBases().length - reverseClipping);
                if ( !Allele.acceptableAlleleBases(allele) )
                    return new ClippedAlleles("Unparsable vcf record with bad allele [" + allele + "]");
                clippedAlleles.add(Allele.create(allele, a.isReference()));
            }
        }

        // the new reference length
        final int refLength = ref.length() - reverseClipping;
        final int stop = start + Math.max(refLength - 1, 0);
        return new ClippedAlleles(stop, clippedAlleles);
    }

    /**
     * Returns true if the alleles in inputVC should have reference bases added for padding
     *
     * We need to pad a VC with a common base if the length of the reference allele is
     * less than the length of the VariantContext. This happens because the position of
     * e.g. an indel is always one before the actual event (as per VCF convention).
     *
     * @param inputVC the VC to evaluate, cannot be null
     * @return true if
     */
    public static boolean needsPadding(final VariantContext inputVC) {
        // biallelic sites with only symbolic never need padding
        if ( inputVC.isBiallelic() && inputVC.getAlternateAllele(0).isSymbolic() )
            return false;

        final int recordLength = inputVC.getEnd() - inputVC.getStart() + 1;
        final int referenceLength = inputVC.getReference().length();

        if ( referenceLength == recordLength )
            return false;
        else if ( referenceLength == recordLength - 1 )
            return true;
        else if ( !inputVC.hasSymbolicAlleles() )
            throw new IllegalArgumentException("Badly formed variant context at location " + String.valueOf(inputVC.getStart()) +
                    " in contig " + inputVC.getChr() + ". Reference length must be at most one base shorter than location size");
        else
            return false;
    }

    public static Allele padAllele(final VariantContext vc, final Allele allele) {
        assert needsPadding(vc);

        if ( allele.isSymbolic() )
            return allele;
        else {
            // get bases for current allele and create a new one with trimmed bases
            final StringBuilder sb = new StringBuilder();
            sb.append((char)vc.getReferenceBaseForIndel().byteValue());
            sb.append(allele.getDisplayString());
            final String newBases = sb.toString();
            return Allele.create(newBases, allele.isReference());
        }
    }

    public static VariantContext createVariantContextWithPaddedAlleles(VariantContext inputVC) {
        final boolean padVC = needsPadding(inputVC);

        // nothing to do if we don't need to pad bases
        if ( padVC ) {
            if ( !inputVC.hasReferenceBaseForIndel() )
                throw new ReviewedStingException("Badly formed variant context at location " + inputVC.getChr() + ":" + inputVC.getStart() + "; no padded reference base is available.");

            final ArrayList<Allele> alleles = new ArrayList<Allele>(inputVC.getNAlleles());
            final Map<Allele, Allele> unpaddedToPadded = new HashMap<Allele, Allele>(inputVC.getNAlleles());

            for (final Allele a : inputVC.getAlleles()) {
                final Allele padded = padAllele(inputVC, a);
                alleles.add(padded);
                unpaddedToPadded.put(a, padded);
            }

            // now we can recreate new genotypes with trimmed alleles
            GenotypesContext genotypes = GenotypesContext.create(inputVC.getNSamples());
            for (final Genotype g : inputVC.getGenotypes() ) {
                final List<Allele> newGenotypeAlleles = new ArrayList<Allele>(g.getAlleles().size());
                for (final Allele a : g.getAlleles()) {
                    newGenotypeAlleles.add( a.isCalled() ? unpaddedToPadded.get(a) : Allele.NO_CALL);
                }
                genotypes.add(new GenotypeBuilder(g).alleles(newGenotypeAlleles).make());

            }

            return new VariantContextBuilder(inputVC).alleles(alleles).genotypes(genotypes).make();
        }
        else
            return inputVC;

    }

    public static VariantContext reverseTrimAlleles( final VariantContext inputVC ) {
        // see if we need to trim common reference base from all alleles

        final int trimExtent = computeReverseClipping(inputVC.getAlleles(), inputVC.getReference().getDisplayString().getBytes(), 0, true);
        if ( trimExtent <= 0 || inputVC.getAlleles().size() <= 1 )
            return inputVC;

        final List<Allele> alleles = new ArrayList<Allele>();
        final GenotypesContext genotypes = GenotypesContext.create();
        final Map<Allele, Allele> originalToTrimmedAlleleMap = new HashMap<Allele, Allele>();

        for (final Allele a : inputVC.getAlleles()) {
            if (a.isSymbolic()) {
                alleles.add(a);
                originalToTrimmedAlleleMap.put(a, a);
            } else {
                // get bases for current allele and create a new one with trimmed bases
                final byte[] newBases = Arrays.copyOfRange(a.getBases(), 0, a.length()-trimExtent);
                final Allele trimmedAllele = Allele.create(newBases, a.isReference());
                alleles.add(trimmedAllele);
                originalToTrimmedAlleleMap.put(a, trimmedAllele);
            }
        }

        // now we can recreate new genotypes with trimmed alleles
        for ( final Genotype genotype : inputVC.getGenotypes() ) {
            final List<Allele> originalAlleles = genotype.getAlleles();
            final List<Allele> trimmedAlleles = new ArrayList<Allele>();
            for ( final Allele a : originalAlleles ) {
                if ( a.isCalled() )
                    trimmedAlleles.add(originalToTrimmedAlleleMap.get(a));
                else
                    trimmedAlleles.add(Allele.NO_CALL);
            }
            genotypes.add(new GenotypeBuilder(genotype).alleles(trimmedAlleles).make());
        }

        return new VariantContextBuilder(inputVC).stop(inputVC.getStart() + alleles.get(0).length() + (inputVC.isMixed() ? -1 : 0)).alleles(alleles).genotypes(genotypes).make();
    }

    public static class ClippedAlleles {
        final int stop;
        final List<Allele> clippedAlleles;
        final String error;

        public ClippedAlleles(final int stop, final List<Allele> clippedAlleles) {
            this.stop = stop;
            this.clippedAlleles = clippedAlleles;
            this.error = null;
        }

        public ClippedAlleles(final String error) {
            this.stop = -1;
            this.clippedAlleles = null;
            this.error = error;
        }

        public String getError() {
            return error;
        }

        public int getStop() {
            return stop;
        }

        public List<Allele> getClippedAlleles() {
            return clippedAlleles;
        }
    }
}
