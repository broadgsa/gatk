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

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.util.ArrayList;
import java.util.List;

public class AllelePair {
    private Allele top;
    private Allele bottom;

    public AllelePair(Genotype gt) {
        if (gt.getPloidy() != 2)
            throw new ReviewedStingException("AllelePair must have ploidy of 2!");

        this.top = gt.getAllele(0);
        this.bottom = gt.getAllele(1);
    }

    public Allele getTopAllele() {
        return top;
    }

    public Allele getBottomAllele() {
        return bottom;
    }

    public void swapAlleles() {
        Allele tmp = top;
        top = bottom;
        bottom = tmp;
    }

    public List<Allele> getAllelesAsList() {
        List<Allele> allList = new ArrayList<Allele>(2);
        allList.add(0, top);
        allList.add(1, bottom);
        return allList;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Top:\t" + top.getBaseString() + "\n");
        sb.append("Bot:\t" + bottom.getBaseString() + "\n");
        return sb.toString();
    }
}