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
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
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
    // Inputs
    /////////////////////////////
    @Input(fullName="tranches_file", shortName="tranchesFile", doc="The input tranches file describing where to cut the data", required=true)
    private File TRANCHES_FILE;

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Output( doc="The output filtered VCF file", required=true)
    private VCFWriter vcfWriter = null;

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="fdr_filter_level", shortName="fdr_filter_level", doc="The FDR level at which to start filtering.", required=false)
    private double FDR_FILTER_LEVEL = 0.0;

    /////////////////////////////
    // Debug Arguments
    /////////////////////////////
    @Hidden
    @Argument(fullName = "NO_HEADER", shortName = "NO_HEADER", doc = "Don't output the usual VCF header tag with the command line. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.", required = false)
    protected Boolean NO_VCF_HEADER_LINE = false;
    
    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
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
            throw new UserException.CouldNotCreateOutputFile(f, e);
        }
    }

    public void initialize() {

        // todo -- ryan, it's always best to use a data structure, I need to read these in too.
        boolean firstLine = true;
        try {
            for( final String line : new XReadLines( TRANCHES_FILE ) ) {
                if( !firstLine ) {
                    final String[] vals = line.split(",");
                    double FDR = Double.parseDouble(vals[0]);
                    double TsTv = Double.parseDouble(vals[1]);
                    double pCut = Double.parseDouble(vals[2]);
                    String name = vals[4];
                    //String statusMsg = "Excluding, below FDR level";
                    if (FDR >= FDR_FILTER_LEVEL) {
                        qCuts.add(pCut);
                        filterName.add(name);
                        //statusMsg = "Keeping, above FDR threshold";
                    }
                    logger.info(String.format("Tranche %s with %.2f FDR, TsTv %.2f and pCut %.2f, threshold %.2f",
                            name, FDR, TsTv, pCut, FDR_FILTER_LEVEL));
                }
                firstLine = false;
            }
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(TRANCHES_FILE, e);
        }

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        if( !NO_VCF_HEADER_LINE ) {
            hInfo.add(new VCFHeaderLine("ApplyVariantCuts", "\"" + CommandLineUtils.createApproximateCommandLineArgumentString(getToolkit(), this) + "\""));
        }
        final TreeSet<String> samples = new TreeSet<String>();
        samples.addAll(SampleUtils.getUniqueSamplesFromRods(getToolkit()));

        if( filterName.size() >= 2 ) {
            for( int iii = 0; iii < filterName.size() - 1; iii++ ) {
                hInfo.add(new VCFFilterHeaderLine(filterName.get(iii), String.format("FDR tranche level at qual: " + qCuts.get(iii) + " <= x < " + qCuts.get(iii+1))));
            }
        }
        if( filterName.size() >= 1 ) {
            hInfo.add(new VCFFilterHeaderLine(filterName.get(0)+"+", String.format("FDR tranche level at qual < " + qCuts.get(0))));
        }

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
                    if( filterString == null ) {
                        filterString = filterName.get(0)+"+";
                    }

                    if ( !filterString.equals(VCFConstants.PASSES_FILTERS_v4) ) {
                        Set<String> filters = new HashSet<String>();
                        filters.add(filterString);
                        vc = VariantContext.modifyFilters(vc, filters);
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
    }
}

