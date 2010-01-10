package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.io.File;

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

    @Argument(fullName="rod_priority_list", shortName="priority", doc="For the UNION combination type: a comma-separated string describing the priority ordering for the rods as far as which record gets emitted; if not specified, then the first rod seen will be used", required=false)
    protected String PRIORITY_STRING = null;


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

        Set<String> samples;
        switch (COMBO_TYPE ) {
            case MERGE:
                samples = new TreeSet<String>();
                rodNamesToSampleNames = new HashMap<Pair<String, String>, String>();
                // get the list of all sample names from the various input rods (they need to be uniquified in case there's overlap)
                SampleUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, rodNamesToSampleNames);
                break;
            case UNION:
                samples = SampleUtils.getUniqueSamplesFromRods(getToolkit());
                break;
            default:
                throw new IllegalArgumentException("Unsupported combination type: " + COMBO_TYPE);
        }

        vcfWriter.writeHeader(new VCFHeader(hInfo, samples));

        if ( PRIORITY_STRING != null )
            priority = PRIORITY_STRING.split(",");
    }

    public VCFRecord map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        // get all of the vcf rods at this locus
        ArrayList<RodVCF> vcfRods = new ArrayList<RodVCF>();
        Iterator<ReferenceOrderedDatum> rods = tracker.getAllRods().iterator();
        while (rods.hasNext()) {
            ReferenceOrderedDatum rod = rods.next();
            if ( rod instanceof RodVCF )
                vcfRods.add((RodVCF)rod);
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

    private VCFRecord vcfUnion(ArrayList<RodVCF> rods) {
        if ( priority == null )
            return rods.get(0).mCurrentRecord;

        for ( String rodname : priority ) {
            for ( RodVCF rod : rods ) {
                if ( rod.getName().equals(rodname) )
                    return rod.mCurrentRecord;
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