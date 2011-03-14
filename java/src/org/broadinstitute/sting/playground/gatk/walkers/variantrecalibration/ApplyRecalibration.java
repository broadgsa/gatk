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

package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.lang.Double;
import java.util.*;

/**
 * Applies cuts to the input vcf file (by adding filter lines) to achieve the desired novel FDR levels which were specified during VariantRecalibration
 *
 * @author rpoplin
 * @since Mar 14, 2011
 *
 * @help.summary Applies cuts to the input vcf file (by adding filter lines) to achieve the desired novel FDR levels which were specified during VariantRecalibration
 */

public class ApplyRecalibration extends RodWalker<Integer, Integer> {


    /////////////////////////////
    // Inputs
    /////////////////////////////
    @Output(fullName="recal_file", shortName="recalFile", doc="The output recal file used by ApplyRecalibration", required=true)
    private File RECAL_FILE;
    @Input(fullName="tranches_file", shortName="tranchesFile", doc="The input tranches file describing where to cut the data", required=true)
    private File TRANCHES_FILE;

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Output( doc="The output filtered, recalibrated VCF file", required=true)
    private VCFWriter vcfWriter = null;

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="fdr_filter_level", shortName="fdr_filter_level", doc="The FDR level at which to start filtering.", required=false)
    private double FDR_FILTER_LEVEL = 1.0;
    
    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    final private List<Tranche> tranches = new ArrayList<Tranche>();
    final private Set<String> inputNames = new HashSet<String>();
    final private NestedHashMap lodMap = new NestedHashMap();

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        for ( Tranche t : Tranche.readTranches(TRANCHES_FILE) ) {
            if ( t.fdr >= FDR_FILTER_LEVEL) {
                tranches.add(t);
                //statusMsg = "Keeping, above FDR threshold";
            }
            logger.info(String.format("Read tranche " + t));
        }
        Collections.reverse(tranches); // this algorithm wants the tranches ordered from worst to best

        for( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if( d.getName().startsWith("input") ) {
                inputNames.add(d.getName());
                logger.info("Found input variant track with name " + d.getName());
            } else {
                logger.info("Not evaluating ROD binding " + d.getName());
            }
        }

        if( inputNames.size() == 0 ) {
            throw new UserException.BadInput( "No input variant tracks found. Input variant binding names must begin with 'input'." );
        }

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit(), inputNames));
        final TreeSet<String> samples = new TreeSet<String>();
        samples.addAll(SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames));

        if( tranches.size() >= 2 ) {
            for( int iii = 0; iii < tranches.size() - 1; iii++ ) {
                Tranche t = tranches.get(iii);
                hInfo.add(new VCFFilterHeaderLine(t.name, String.format("FDR tranche level at VSQ Lod: " + t.minVQSLod + " <= x < " + tranches.get(iii+1).minVQSLod)));
            }
        }
        if( tranches.size() >= 1 ) {
            hInfo.add(new VCFFilterHeaderLine(tranches.get(0).name + "+", String.format("FDR tranche level at VQS Lod < " + tranches.get(0).minVQSLod)));
        } else {
            throw new UserException("No tranches were found in the file or were above the FDR Filter level " + FDR_FILTER_LEVEL);
        }

        logger.info("Keeping all variants in tranche " + tranches.get(tranches.size()-1));

        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);

        try {
            logger.info("Reading in recalibration table...");
            for ( final String line : new XReadLines( RECAL_FILE ) ) {
                final String[] vals = line.split(",");
                lodMap.put( Double.parseDouble(vals[2]), vals[0], Integer.parseInt(vals[1]));
            }
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(RECAL_FILE, e);
        }

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

        for( VariantContext vc : tracker.getVariantContexts(ref, inputNames, null, context.getLocation(), true, false) ) {
            if( vc != null ) {
                String filterString = null;
                final Double lod = (Double) lodMap.get( ref.getLocus().getContig(), ref.getLocus().getStart() );
                if( !vc.isFiltered() ) {
                    try {
                        for( int i = tranches.size() - 1; i >= 0; i-- ) {
                            Tranche tranche = tranches.get(i);
                            if( lod >= tranche.minVQSLod ) {
                                if (i == tranches.size() - 1) {
                                    filterString = VCFConstants.PASSES_FILTERS_v4;
                                } else {
                                    filterString = tranche.name;
                                }
                                break;
                            }
                        }

                        if( filterString == null ) {
                            filterString = tranches.get(0).name+"+";
                        }

                        if ( !filterString.equals(VCFConstants.PASSES_FILTERS_v4) ) {
                            Set<String> filters = new HashSet<String>();
                            filters.add(filterString);
                            vc = VariantContext.modifyFilters(vc, filters);
                        }
                    } catch ( Exception e ) {
                        throw new UserException.MalformedFile(vc.getSource(), "Invalid value for VQS key " + ContrastiveRecalibrator.VQS_LOD_KEY + " at variant " + vc, e);
                    }
                }

                final Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());

                if(lod != null) {
                    attrs.put(ContrastiveRecalibrator.VQS_LOD_KEY, String.format("%.4f", lod));
                }
                VariantContext newVC = VariantContext.modifyPErrorFiltersAndAttributes(vc, vc.getNegLog10PError(), new HashSet<String>(), attrs);
                vcfWriter.add( newVC, ref.getBase() );
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
        return 1; // This value isn't used for anything
    }

    public Integer reduce( final Integer mapValue, final Integer reduceSum ) {
        return 1; // This value isn't used for anything
    }

    public void onTraversalDone( Integer reduceSum ) {
    }
}

