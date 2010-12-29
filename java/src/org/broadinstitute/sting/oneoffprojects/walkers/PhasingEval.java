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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits specific fields as dictated by the user from one or more VCF files.
 */
@Requires(value={})
public class PhasingEval extends RodWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(doc="sample to emit", required = false)
    protected String sample = null;

    @Argument(doc="only include physical phased results", required = false)
    protected boolean requirePQ = false;

    @Argument(doc="Analysis to perform", required = true)
    protected Analysis analysis;

    public enum Analysis {
        HAPLOTYPE_SIZES,
        PHASING_BY_AC
    }

    private class Haplotype {
        GenomeLoc start, last;
        List<Genotype> genotypes;

        public Haplotype(GenomeLoc start) {
            this.start = start;
            this.last = null;
            this.genotypes = new ArrayList<Genotype>();
        }
    }

    private class PhasingByAC {
        int myAC = 0;
        int myAN = 0;
        int nHets = 0;
        int nHetsPhased = 0;

        public PhasingByAC(int myAC, int myAN) {
            this.myAC = myAC;
            this.myAN = myAN;
        }
    }

    Map<String, Haplotype> haplotypes = new Hashtable<String, Haplotype>();
    List<PhasingByAC> phasingByACs = new ArrayList<PhasingByAC>();

    public void initialize() {
        if ( analysis == Analysis.HAPLOTYPE_SIZES ) {
            out.println(Utils.join("\t", Arrays.asList("sample", "haplotype.length", "n.genotypes", "start", "stop")));
        }

        Set<String> samples = SampleUtils.getSampleList(VCFUtils.getVCFHeadersFromRods(getToolkit(), null));
        int AN = 2 * samples.size();
        for ( int i = 0; i <= AN; i++ ) {
            phasingByACs.add(new PhasingByAC(i, AN));
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref, context.getLocation());
        for ( VariantContext vc : vcs) {
            if ( sample != null )
                vc = vc.subContextFromGenotypes(vc.getGenotype(sample));

            if ( analysis == Analysis.HAPLOTYPE_SIZES) {
                GenomeLoc loc = getToolkit().getGenomeLocParser().parseGenomeLoc(vc.getChr(), vc.getStart(), vc.getEnd());
                for ( Genotype g : vc.getGenotypes().values() ) {
                    if ( ! haplotypes.containsKey(g.getSampleName()) )
                        haplotypes.put(g.getSampleName(), new Haplotype(loc));

                    Haplotype h = haplotypes.get(g.getSampleName());

                    if ( g.isPhased() ) {
                        h.genotypes.add(g);
                        h.last = loc;
                    } else {
                        if ( ! requirePQ || isPhysicallyPhased(h.genotypes) )
                            out.printf("%s %d %d %s %s%n", g.getSampleName(),
                                    h.last != null ? h.start.distance(h.last) : 0,
                                    h.genotypes.size(), h.start, h.last);
                        haplotypes.put(g.getSampleName(), new Haplotype(loc));
                    }
                }
            } else if ( analysis == Analysis.PHASING_BY_AC ) {
                int homref = vc.getHomRefCount();
                int homalt = vc.getHomVarCount();
                int het = vc.getHetCount();
                int ac = 2 * homalt + het;
                //int an = 2 * (homref + homalt + het);
                PhasingByAC data = phasingByACs.get(ac);
                data.nHets += het > 0 ? 1 : 0;
                data.nHetsPhased += isPhysicallyPhased(vc.getGenotypes().values()) ? 1 : 0;
            }
        }

        return 1;
    }

    private boolean isPhysicallyPhased(Collection<Genotype> genotypes) {
        for ( Genotype g : genotypes ) {
            if ( g.isHet() && g.hasAttribute("PQ") )
                return true;
        }

        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        if ( analysis == Analysis.PHASING_BY_AC ) {
            out.println(Utils.join("\t", Arrays.asList("ac", "nhets", "nhetphased")));
            for ( PhasingByAC pac : phasingByACs ) {
                out.printf("%d\t%d\t%d%n", pac.myAC, pac.nHets, pac.nHetsPhased);
            }
        }
    }
}
