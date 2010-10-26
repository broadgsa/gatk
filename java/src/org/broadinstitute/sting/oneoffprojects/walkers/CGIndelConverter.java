/*
 * Copyright (c) 2010.
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

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;


/**
 * Converts CG indels to VCF format
 */
public class CGIndelConverter extends RefWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter writer = null;

    @Argument(fullName="raw_cg_rod_file", shortName="cg", doc="file with raw cg calls", required=true)
    private File CG_FILE = null;

    @Argument(fullName="sample_name", shortName="sample", doc="sample represented by the calls", required=true)
    private String sample = null;

    private HashMap<String, IndelCall> calls = new HashMap<String, IndelCall>();

    private static class IndelCall {
        public String chr, pos, bases, type, id;
        public int count;

        public IndelCall(StringTokenizer st) {
            chr = st.nextToken();
            pos = st.nextToken();
            bases = st.nextToken();
            type = st.nextToken();
            if ( st.hasMoreTokens() )
                id = st.nextToken();
            count = 1;
        }
    }

    public Integer reduceInit() { return 0; }

    public void initialize() {
        // set up the vcf writer
        Set<String> samples = new HashSet<String>();
        samples.add(sample);
        writer.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), samples));

        // read in the raw calls file
        try {
            AsciiLineReader lineReader = new AsciiLineReader(new FileInputStream(CG_FILE));
            String line = lineReader.readLine();
            while ( line != null ) {
                StringTokenizer st = new StringTokenizer(line);
                IndelCall call = new IndelCall(st);
                String key = call.chr + ":" + call.pos;
                if ( calls.containsKey(key) )
                    calls.get(key).count++;
                else
                    calls.put(key, call);

                line = lineReader.readLine();
            }
        }
        catch (IOException e ) {
            throw new ReviewedStingException(e.getMessage());
        }
    }

    /**
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        String pos = ref.getLocus().getContig() + ":" + ref.getLocus().getStart();
        if ( !calls.containsKey(pos) ) {
            System.out.println("We don't have any calls for position " + pos + " but somehow triggered here");
            return 0;
        }

        IndelCall call = calls.get(pos);

        Set<Allele> alleles = new HashSet<Allele>();
        Map<String, Genotype> genotypes = new HashMap<String, Genotype>();

        Allele refAllele, altAllele;
        if ( call.type.equals("del") ) {
            refAllele = Allele.create(call.bases, true);
            altAllele = Allele.create(Allele.NULL_ALLELE_STRING, false);
        } else if ( call.type.equals("ins") ) {
            refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
            altAllele = Allele.create(call.bases, false);        
        } else {
            System.out.println("Weird type seen for position " + pos + " -> " + call.type);
            return 0;
        }

        alleles.add(refAllele);
        alleles.add(altAllele);
        long end = ref.getLocus().getStart() + (call.type.equals("del") ? call.bases.length() : 0);

        List<Allele> myAlleles = new ArrayList<Allele>();
        myAlleles.add(altAllele);
        if ( call.count == 1 )
            myAlleles.add(refAllele);
        else
            myAlleles.add(altAllele);
        genotypes.put(sample, new Genotype(sample, myAlleles, -1));

        Map<String, Object> attrs = new HashMap<String, Object>();
        if ( call.id != null ) {
            String rsID = call.id.substring(call.id.lastIndexOf(":")+1);
            attrs.put(VariantContext.ID_KEY, rsID);
        }

        VariantContext vc = new VariantContext("CG",
                ref.getLocus().getContig(),
                ref.getLocus().getStart(),
                end,
                alleles,
                genotypes,
                -1,
                null,
                attrs);

        writer.add(vc, ref.getBase());

        return 1;
    }

    /**
     * Increment the number of rods processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of rods processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {}
}