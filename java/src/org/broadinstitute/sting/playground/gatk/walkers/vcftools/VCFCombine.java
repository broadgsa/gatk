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

package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.util.*;

/**
 * Combines VCF records from different sources; supports both full merges and set unions.
 * Merge: combines multiple records into a single one; if sample names overlap then they are uniquified.
 * Union: assumes each rod represents the same set of samples (although this is not enforced); using the
 *   priority list (if provided), emits a single record instance at every position represented in the rods.
 */
@Requires(value={})
public class VCFCombine extends RodWalker<VCFRecord, VCFWriter> {
    // the types of combinations we currently allow
    public enum COMBINATION_TYPE {
        UNION, MERGE
    }

    @Argument(fullName="vcf_output_file", shortName="O", doc="VCF file to write results", required=false)
    protected File OUTPUT_FILE = null;

    @Argument(fullName="combination_type", shortName="type", doc="combination type; currently UNION and MERGE are supported", required=true)
    protected COMBINATION_TYPE COMBO_TYPE;

    @Argument(fullName="rod_priority_list", shortName="priority", doc="For the UNION combination type: a comma-separated string describing the priority ordering for the rods as far as which record gets emitted; a complete priority list MUST be provided", required=false)
    protected String PRIORITY_STRING = null;

    @Argument(fullName="annotateUnion", shortName="A", doc="For the UNION combination type: if provided, the output union VCF will contain venn information about where each call came from", required=false)
    protected boolean annotateUnion = false;

    private VCFWriter vcfWriter = null;
    private String[] priority = null;

    // a map of rod name to uniquified sample name
    private HashMap<Pair<String, String>, String> rodNamesToSampleNames;


    public void initialize() {

        Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
        metaData.add(new VCFHeaderLine("source", "VCF"));

        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));

        if ( OUTPUT_FILE != null )
            vcfWriter = new VCFWriter(OUTPUT_FILE);
        else
            vcfWriter = new VCFWriter(out);

        if ( PRIORITY_STRING != null )
            priority = PRIORITY_STRING.split(",");

        Set<String> samples;
        switch (COMBO_TYPE ) {
            case MERGE:
                samples = new TreeSet<String>();
                rodNamesToSampleNames = new HashMap<Pair<String, String>, String>();
                // get the list of all sample names from the various input rods (they need to be uniquified in case there's overlap)
                SampleUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, rodNamesToSampleNames);
                break;
            case UNION:
                if ( annotateUnion ) {
                    validateAnnotateUnionArguments(priority);
                }

                samples = SampleUtils.getUniqueSamplesFromRods(getToolkit());
                break;
            default:
                throw new IllegalArgumentException("Unsupported combination type: " + COMBO_TYPE);
        }

        vcfWriter.writeHeader(new VCFHeader(hInfo, samples));
    }

    private void validateAnnotateUnionArguments(String[] priority) {
        Set<RMDTrack> rods = VCFUtils.getRodVCFs(getToolkit());
        if ( rods.size() != priority.length ) {
            throw new StingException("A complete priority list must be provided when annotateUnion is provided");
        }
        if ( priority.length != 2 ) {
            throw new StingException("When annotateUnion is provided only 2 VCF files can be merged");
        }

        for ( String p : priority ) {
            boolean good = false;
            for ( RMDTrack data : rods ) {
                if ( p.equals(data.getName()) )
                    good = true;
            }
            if ( ! good ) throw new StingException("Priority item not provided as a ROD! " + p);
        }
    }

    public VCFRecord map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        // get all of the vcf rods at this locus
        Map<VCFRecord, String> vcfRods = new LinkedHashMap<VCFRecord,String>();
        Iterator<GATKFeature> rods = tracker.getAllRods().iterator();
        while (rods.hasNext()) {
            GATKFeature feat = rods.next();
            Object rod = feat.getUnderlyingObject();
            if ( rod instanceof VCFRecord )
                vcfRods.put((VCFRecord)rod,feat.getName());
        }

        if ( vcfRods.size() == 0 )
            return null;

        VCFRecord record = null;
        switch (COMBO_TYPE ) {
            case MERGE:
                record = VCFUtils.mergeRecords(vcfRods, rodNamesToSampleNames);
                break;
            case UNION:
                record = vcfUnion(vcfRods);
                break;
            default:
                break;
        }

        return record;
    }

    public VCFWriter reduceInit() {
        return vcfWriter;
    }

    private VCFRecord vcfUnion(Map<VCFRecord, String> rods) {
        if ( priority == null )
            return rods.keySet().iterator().next();

        if ( annotateUnion ) {
            Map<String, VCFRecord> rodMap = new HashMap<String, VCFRecord>();
            for ( VCFRecord vcf : rods.keySet() ) {
                rodMap.put(rods.get(vcf),vcf);
            }

            String priority1 = priority[0];
            String priority2 = priority[1];
            VCFRecord vcf1 = rodMap.containsKey(priority1) ? rodMap.get(priority1) : null;
            VCFRecord vcf2 = rodMap.containsKey(priority2) ? rodMap.get(priority2) : null;

            // for simplicity, we are setting set and call for vcf1
            String set = priority1;
            VCFRecord call = vcf1;

            if ( vcf1 == null ) {
                if ( vcf2 == null )
                    //return null;
                    throw new StingException("BUG: VCF1 and VCF2 are both null!");
                else {
                    set = priority2;
                    call = vcf2;
                }
            } else if ( vcf1.isFiltered() ) {
                if ( vcf2 != null ) {
                    if ( vcf2.isFiltered() ) {
                        set = "filteredInBoth";
                    } else {
                        set = priority2 + "-filteredInOther";
                        call = vcf2;
                    }
                }
            } else { // good call
                if ( vcf2 != null ) {
                    if ( vcf2.isFiltered() )
                        set = priority1 + "-filteredInOther";
                    else
                        set = "Intersection";
                }
            }

            call.addInfoField("set", set);
            return call;
        } else {
            for ( String rodname : priority ) {
                for ( VCFRecord rod : rods.keySet() ) {
                    if ( rods.get(rod).equals(rodname) )
                        return rod;
                }
            }
        }

        return null;
    }

    public VCFWriter reduce(VCFRecord record, VCFWriter writer) {
        if ( record != null )
            writer.addRecord(record);
        return writer;
    }

    public void onTraversalDone(VCFWriter writer) {
        if ( writer != null )
            writer.close();
    }
}