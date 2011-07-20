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
package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

public class SNPallelePair extends AllelePair {

    public SNPallelePair(Genotype gt) {
        super(gt);

        if (getTopAllele().getBases().length != 1)
            throw new ReviewedStingException("LOGICAL ERROR: SNPallelePair may not contain non-SNP site!");
        if (getBottomAllele().getBases().length != 1)
            throw new ReviewedStingException("LOGICAL ERROR: SNPallelePair may not contain non-SNP site!");
    }

    public byte getTopBase() {
        byte[] topBases = getTopAllele().getBases();
        return getSingleBase(topBases);
    }

    public byte getBottomBase() {
        byte[] bottomBases = getBottomAllele().getBases();
        return getSingleBase(bottomBases);
    }

    public boolean matchesTopBase(byte base) {
        return BaseUtils.basesAreEqual(base, getTopBase());
    }

    public byte getOtherBase(byte base) {
        byte topBase = getTopBase();
        byte botBase = getBottomBase();

        if (BaseUtils.basesAreEqual(base, topBase))
            return botBase;
        else if (BaseUtils.basesAreEqual(base, botBase))
            return topBase;
        else
            throw new ReviewedStingException("LOGICAL ERROR: base MUST match either TOP or BOTTOM!");
    }

    public static byte getSingleBase(byte[] bases) {
        return bases[0];
    }

    public static byte getSingleBase(Allele all) {
        return getSingleBase(all.getBases());
    }
}