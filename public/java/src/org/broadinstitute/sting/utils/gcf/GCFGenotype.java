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

package org.broadinstitute.sting.utils.gcf;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * GATK binary VCF record
 *
 * @author Your Name
 * @since Date created
 */
public class GCFGenotype {
    private byte gq;
    private int gt;
    private int dp;
    private int ad[];
    private byte[] pl;

    // todo -- what to do about phasing?  Perhaps we shouldn't support it
    // todo -- is the FL field generic or just a flag?  Should we even support per sample filtering?

    public GCFGenotype(final GCFHeaderBuilder GCFHeaderBuilder, final List<Allele> allAlleles, Genotype genotype) {
        gq = GCF.qualToByte(genotype.getPhredScaledQual());
        gt = encodeAlleles(genotype.getAlleles(), allAlleles);

        dp = genotype.getAttributeAsInt("DP", 0);

        int nAlleles = allAlleles.size();
        ad = new int[nAlleles];

        int npls = nAllelesToNPls(nAlleles);
        pl = new byte[npls];
    }

    private int nAllelesToNPls( int nAlleles ) {
        return nAlleles*(nAlleles+1) / 2;
    }

    public GCFGenotype(GCF GCF, DataInputStream inputStream) throws IOException {
        int gqInt = inputStream.readUnsignedByte();
        gq = (byte)gqInt;
        gt = inputStream.readInt();
        dp = inputStream.readInt();
        ad = GCF.readIntArray(inputStream, GCF.getNAlleles());
        pl = GCF.readByteArray(inputStream, nAllelesToNPls(GCF.getNAlleles()));
    }

    // 2 alleles => 1 + 8 + 8 + 3 => 20
    protected int sizeInBytes() {
        return 1 // gq
                + 4 * 2 // gt + dp
                + 4 * ad.length // ad
                + 1 * pl.length; // pl
    }

    public Genotype decode(final String sampleName, final GCFHeader header, GCF GCF, List<Allele> alleleIndex) {
        final List<Allele> alleles = decodeAlleles(gt, alleleIndex);
        final double negLog10PError = gq / 10.0;
        final Set<String> filters = Collections.emptySet();
        final Map<String, Object> attributes = new HashMap<String, Object>();
        attributes.put("DP", dp);
        attributes.put("AD", ad);
        attributes.put("PL", pl);

        return new Genotype(sampleName, alleles, negLog10PError, filters, attributes, false);
    }

    private static int encodeAlleles(List<Allele> gtList, List<Allele> allAlleles) {
        final int nAlleles = gtList.size();
        if ( nAlleles  > 4 )
            throw new IllegalArgumentException("encodeAlleles doesn't support more than 4 alt alleles, but I saw " + gtList);

        int gtInt = 0;
        for ( int i = 0; i < nAlleles ; i++ ) {
            final int bitOffset = i * 8;
            final int allelei = getAlleleIndex(gtList.get(i), allAlleles);
            final int gti = (allelei + 1) << bitOffset;
            gtInt = gtInt | gti;
        }

        return gtInt;
    }

    private static int getAlleleIndex(Allele q, List<Allele> allAlleles) {
        if ( q.isNoCall() )
            return 254;
        for ( int i = 0; i < allAlleles.size(); i++ )
            if ( q.equals(allAlleles.get(i)) )
                return i;
        throw new IllegalStateException("getAlleleIndex passed allele not in map! allele " + q + " allAlleles " + allAlleles);
    }

    private static List<Allele> decodeAlleles(int gtInt, List<Allele> alleleIndex) {
        List<Allele> alleles = new ArrayList<Allele>(4);

        for ( int i = 0; i < 32; i += 8 ) {
            final int gi = (gtInt & (0x000000FF << i)) >> i;
            if ( gi != 0 ) {
                final int allelei = gi - 1;
                alleles.add( allelei == 254 ? Allele.NO_CALL : alleleIndex.get(allelei) );
            } else {
                break;
            }
        }

        return alleles;
    }

    public int write(DataOutputStream outputStream) throws IOException {
        int startSize = outputStream.size();
        outputStream.writeByte(gq);
        outputStream.writeInt(gt);
        outputStream.writeInt(dp);
        GCF.writeIntArray(ad, outputStream, false);
        GCF.writeByteArray(pl, outputStream, false);
        return outputStream.size() - startSize;
    }
}
