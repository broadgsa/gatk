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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
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
public class GCFHeader {
    final protected static Logger logger = Logger.getLogger(GCFHeader.class);

    private static byte[] MAGIC_HEADER = "GVCF0.1\1".getBytes();
    final List<Allele> alleles;
    final List<String> strings;
    final List<String> samples;
    final List<Set<String>> filters;

    public GCFHeader(final Map<Allele, Integer> allelesIn, final Map<String, Integer> stringIn, final Map<String, Integer> samplesIn) {
        this.alleles = linearize(allelesIn);
        this.strings = linearize(stringIn);
        this.samples = linearize(samplesIn);
        this.filters = null; // not used with this constructor
    }

    public GCFHeader(DataInputStream inputStream) throws IOException {
        byte[] headerTest = new byte[MAGIC_HEADER.length];
        inputStream.read(headerTest);
        if ( ! Arrays.equals(headerTest, MAGIC_HEADER) ) {
            throw new UserException("Could not read GVCF file.  MAGIC_HEADER missing.  Saw " + headerTest);
        } else {
            alleles = stringsToAlleles(readStrings(inputStream));
            strings = readStrings(inputStream);
            samples = readStrings(inputStream);
            logger.info(String.format("Allele map of %d elements", alleles.size()));
            logger.info(String.format("String map of %d elements", strings.size()));
            logger.info(String.format("Sample map of %d elements", samples.size()));
            filters = initializeFilterCache();
        }
    }

    public int write(final DataOutputStream outputStream) throws IOException {
        int startBytes = outputStream.size();
        outputStream.write(MAGIC_HEADER);
        write(outputStream, allelesToStrings(alleles));
        write(outputStream, strings);
        write(outputStream, samples);
        return outputStream.size() - startBytes;
    }

    public void write(DataOutputStream outputStream, List<String> l) throws IOException {
        outputStream.writeInt(l.size());
        for ( String elt : l ) outputStream.writeUTF(elt);
    }

    private List<String> allelesToStrings(List<Allele> alleles) {
        List<String> strings = new ArrayList<String>(alleles.size());
        for ( Allele allele : alleles ) strings.add(allele.toString());
        return strings;
    }

    private List<Set<String>> initializeFilterCache() {
        // required to allow offset -> set lookup
        List<Set<String>> l = new ArrayList<Set<String>>(strings.size());
        for ( int i = 0; i < strings.size(); i++ ) l.add(null);
        return l;
    }

    private static List<Allele> stringsToAlleles(final List<String> strings) {
        final List<Allele> alleles = new ArrayList<Allele>(strings.size());
        for ( String string : strings ) {
            boolean isRef = string.endsWith("*");
            if ( isRef ) string = string.substring(0, string.length() - 1);
            alleles.add(Allele.create(string, isRef));
        }
        return alleles;
    }

    private static List<String> readStrings(final DataInputStream inputStream) throws IOException {
        final int nStrings = inputStream.readInt();

        final List<String> strings = new ArrayList<String>(nStrings);
        for ( int i = 0; i < nStrings; i++ ) {
            strings.add(inputStream.readUTF());
        }

        return strings;
    }

    private static <T> List<T> linearize(final Map<T, Integer> map) {
        final ArrayList<T> l = new ArrayList<T>(map.size());
        for ( int i = 0; i < map.size(); i++ ) l.add(null);
        for ( final Map.Entry<T, Integer> elt : map.entrySet() )
            l.set(elt.getValue(), elt.getKey());
        return l;
    }

    public String getSample(final int offset) { return samples.get(offset); }
    public String getString(final int offset) { return strings.get(offset); }
    public Allele getAllele(final int offset) { return alleles.get(offset); }
    public List<Allele> getAlleles(final int[] offsets) {
        final List<Allele> alleles = new ArrayList<Allele>(offsets.length);
        for ( int i : offsets ) alleles.add(getAllele(i));
        return alleles;
    }

    public Set<String> getFilters(final int offset) {
        Set<String> cached = filters.get(offset);

        if ( cached != null )
            return cached;
        else {
            final String filterString = getString(offset);
            if ( filterString.equals(VCFConstants.UNFILTERED) )
                return null; // UNFILTERED records are represented by null
            else {
                Set<String> set = VCFCodec.parseFilters(null, -1, filterString);
                filters.set(offset, set); // remember the result
                return set;
            }
        }
    }
}
