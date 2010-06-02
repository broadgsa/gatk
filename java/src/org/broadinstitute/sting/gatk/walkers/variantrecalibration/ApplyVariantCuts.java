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

import Jama.Matrix;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
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

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VCFWriter vcfWriter;
    private final ArrayList<String> ALLOWED_FORMAT_FIELDS = new ArrayList<String>();
    final ExpandingArrayList<Double> qCuts = new ExpandingArrayList<Double>();
    final ExpandingArrayList<String> filterName = new ExpandingArrayList<String>();

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

        boolean firstLine = true;
        try {
            for( final String line : new XReadLines(new File( TRANCHE_FILENAME )) ) {
                if( !firstLine ) {
                    final String[] vals = line.split(",");
                    qCuts.add(Double.parseDouble(vals[2]));
                    filterName.add(vals[4]);
                }
                firstLine = false;
            }
        } catch( FileNotFoundException e ) {
            throw new StingException("Can not find input file: " + TRANCHE_FILENAME);
        }

        ALLOWED_FORMAT_FIELDS.add(VCFGenotypeRecord.GENOTYPE_KEY); // copied from VariantsToVCF
        ALLOWED_FORMAT_FIELDS.add(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);
        ALLOWED_FORMAT_FIELDS.add(VCFGenotypeRecord.DEPTH_KEY);
        ALLOWED_FORMAT_FIELDS.add(VCFGenotypeRecord.GENOTYPE_LIKELIHOODS_KEY);

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFInfoHeaderLine("OQ", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "The original variant quality score"));
        hInfo.add(new VCFHeaderLine("source", "VariantOptimizer"));
        vcfWriter = new VCFWriter( new File(OUTPUT_FILENAME) );
        final TreeSet<String> samples = new TreeSet<String>();
        final List<ReferenceOrderedDataSource> dataSources = this.getToolkit().getRodDataSources();
        for( final ReferenceOrderedDataSource source : dataSources ) {
            final RMDTrack rod = source.getReferenceOrderedData();
            if( rod.getRecordType().equals(VCFRecord.class) ) {
                final VCFReader reader = new VCFReader(rod.getFile());
                final Set<String> vcfSamples = reader.getHeader().getGenotypeSamples();
                samples.addAll(vcfSamples);
                reader.close();
            }
        }

        // BUGBUG : why doesn't this work?  VCFFilterHeaderLine is protected?
        //for( int iii = 1; iii < filterName.size(); iii++ ) {
        //    hInfo.add(new VCFFilterHeaderLine(filterName.get(iii)), String.format("FDR tranche level at qual " + qCuts.get(iii)));
        //}
        //hInfo.add(new VCFFilterHeaderLine(filterName.get(0)+"+"), String.format("FDR tranche level at qual > " + qCuts.get(filterName.size()-1)));

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

        for( final VariantContext vc : tracker.getAllVariantContexts(ref, null, context.getLocation(), false, false) ) {
            if( vc != null && !vc.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) && vc.isSNP() ) {
                final VCFRecord vcf = VariantContextAdaptors.toVCF(vc, ref.getBase(), ALLOWED_FORMAT_FIELDS, false, false);
                if( !vc.isFiltered() ) {
                    final double qual = vc.getPhredScaledQual();
                    boolean setFilter = false;
                    for( int tranche = qCuts.size() - 1; tranche >= 0; tranche-- ) {
                        if( qual >= qCuts.get(tranche) ) {
                            if(tranche == qCuts.size() - 1) {
                                vcf.setFilterString(VCFRecord.PASSES_FILTERS);
                                setFilter = true;
                            } else {
                                vcf.setFilterString(filterName.get(tranche));
                                setFilter = true;
                            }
                            break;
                        }
                    }
                    if( !setFilter ) {
                        vcf.setFilterString(filterName.get(0)+"+");
                    }
                }
                vcfWriter.addRecord( vcf );
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

    public void onTraversalDone( ExpandingArrayList<VariantDatum> reduceSum ) {
        vcfWriter.close();
    }
}

