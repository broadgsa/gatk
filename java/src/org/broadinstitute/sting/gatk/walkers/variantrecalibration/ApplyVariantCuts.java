/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Applies cuts to the input vcf file (by adding filter lines) to achieve the desired novel FDR levels which were specified during VariantRecalibration
 *
 * @author rpoplin
 * @since Jun 2, 2010
 *
 * @help.summary Applies cuts to the input vcf file (by adding filter lines) to achieve the desired novel FDR levels which were specified during VariantRecalibration
 */

public class ApplyVariantCuts extends RodWalker<Integer, Integer> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="tranchesFile", shortName="tf", doc="The input tranches file describing where to cut the data", required=true)
    private String TRANCHE_FILENAME = "optimizer.dat.tranches";
    @Argument(fullName="outputVCFFile", shortName="outputVCF", doc="The output filtered VCF file", required=true)
    private String OUTPUT_FILENAME = "optimizer.vcf";
    @Argument(fullName="fdr_filter_level", shortName="fdr_filter_level", doc="The FDR level at which to start filtering.", required=false)
    private double FDR_FILTER_LEVEL = 0.0;


    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VCFWriter vcfWriter;
    final ExpandingArrayList<Double> qCuts = new ExpandingArrayList<Double>();
    final ExpandingArrayList<String> filterName = new ExpandingArrayList<String>();

    public static class Tranche {
        public double fdr, pCut, novelTiTv;
        public int numNovel;
        public String name;

        public Tranche(double fdr, double pCut, double novelTiTv, int numNovel, String name) {
            this.fdr = fdr;
            this.pCut = pCut;
            this.novelTiTv = novelTiTv;
            this.numNovel = numNovel;
            this.name = name;
        }

        public Tranche(final String line) {
            final String[] vals = line.split(",");
            this.fdr = Double.parseDouble(vals[0]);
            this.novelTiTv = Double.parseDouble(vals[1]);
            this.pCut = Double.parseDouble(vals[2]);
            this.numNovel = Integer.parseInt(vals[3]);
            this.name = vals[4];
        }

        public String toString() {
            return String.format("[Tranche %s cut = %.3f with %d novels @ %.2f]", name, pCut, numNovel, novelTiTv);
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public static List<Tranche> readTraches(File f) {
        boolean firstLine = true;
        List<Tranche> tranches = new ArrayList<Tranche>();

        try {
            for( final String line : new XReadLines(f) ) {
                if( ! firstLine ) {
                    tranches.add(new Tranche(line));
                }
                firstLine = false;
            }

            return tranches;
        } catch( FileNotFoundException e ) {
            throw new StingException("Can not find input file: " + f);
        }
    }

    public void initialize() {

        // todo -- ryan, it's always best to use a read data structure, I need to read these in.
        // todo -- I would have updated your code but there's no integration test to protect me from unexpected effects
        boolean firstLine = true;
        try {
            for( final String line : new XReadLines(new File( TRANCHE_FILENAME )) ) {
                if( !firstLine ) {
                    final String[] vals = line.split(",");
                    if(Double.parseDouble(vals[0]) >= FDR_FILTER_LEVEL) {
                        qCuts.add(Double.parseDouble(vals[2]));
                        filterName.add(vals[4]);
                    }
                }
                firstLine = false;
            }
        } catch( FileNotFoundException e ) {
            throw new StingException("Can not find input file: " + TRANCHE_FILENAME);
        }

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFInfoHeaderLine("OQ", 1, VCFHeaderLineType.Float, "The original variant quality score"));
        hInfo.add(new VCFHeaderLine("source", "VariantOptimizer"));
        vcfWriter = new VCFWriter( new File(OUTPUT_FILENAME) );
        final TreeSet<String> samples = new TreeSet<String>();
        samples.addAll(SampleUtils.getSampleListWithVCFHeader(getToolkit(), null));
        
        for( int iii = 1; iii < filterName.size(); iii++ ) {
            hInfo.add(new VCFFilterHeaderLine(filterName.get(iii), String.format("FDR tranche level at qual " + qCuts.get(iii))));
        }
        hInfo.add(new VCFFilterHeaderLine(filterName.get(0)+"+", String.format("FDR tranche level at qual > " + qCuts.get(filterName.size()-1))));

        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return 1;
        }

        for( VariantContext vc : tracker.getAllVariantContexts(ref, null, context.getLocation(), false, false) ) {
            if( vc != null && !vc.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) && vc.isSNP() ) {
                String filterString = null;
                if( !vc.isFiltered() ) {
                    final double qual = vc.getPhredScaledQual();
                    for( int tranche = qCuts.size() - 1; tranche >= 0; tranche-- ) {
                        if( qual >= qCuts.get(tranche) ) {
                            if(tranche == qCuts.size() - 1) {
                                filterString = VCFConstants.PASSES_FILTERS_v4;
                            } else {
                                filterString = filterName.get(tranche);
                            }
                            break;
                        }
                    }
                    if( filterString == null )
                        filterString = filterName.get(0)+"+";

                    if ( !filterString.equals(VCFConstants.PASSES_FILTERS_v4) ) {
                        Set<String> filters = new HashSet<String>();
                        filters.add(filterString);
                        vc = new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), vc.getGenotypes(), vc.getNegLog10PError(), filters, vc.getAttributes());
                    }
                }
                vcfWriter.add( vc, ref.getBase() );
            }

        }

        return 1; // This value isn't used for anything
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer reduceInit() {
        return 1;
    }

    public Integer reduce( final Integer mapValue, final Integer reduceSum ) {
        return 1;
    }

    public void onTraversalDone( Integer reduceSum ) {
        vcfWriter.close();
    }
}

