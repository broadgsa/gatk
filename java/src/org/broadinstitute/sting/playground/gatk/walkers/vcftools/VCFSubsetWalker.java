package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.io.File;

/**
 * Extracts subsets of a VCF file like one or more samples, all or only variant loci, all or filtered loci.
 */
public class VCFSubsetWalker extends RefWalker<ArrayList<VCFRecord>, VCFWriter> {
    @Argument(fullName="sample", shortName="SN", doc="Sample to include (or nothing to specify all samples)", required=false)
    private HashSet<String> SAMPLES;

    @Argument(fullName="vcfsubset", shortName="O", doc="File to write VCF subset to", required=false)
    private File VPATH;

    @Argument(fullName="includeNonVariants", shortName="INV", doc="Include non-variant loci", required=false)
    private boolean INCLUDE_NON_VARIANTS = false;

    @Argument(fullName="includeFiltered", shortName="IF", doc="Include filtered loci", required=false)
    private boolean INCLUDE_FILTERED = false;

    private VCFHeader vheader = null;
    private VCFWriter vwriter = null;

    public void initializeWriter() {

        Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
        metaData.add(new VCFHeaderLine("source", "VariantsToVCF"));
        metaData.add(new VCFHeaderLine("reference", this.getToolkit().getArguments().referenceFile.getAbsolutePath()));

        Set<String> additionalColumns = new HashSet<String>();
        additionalColumns.add("FORMAT");
        additionalColumns.addAll(SAMPLES);

        vheader = new VCFHeader(metaData, additionalColumns);
        if (VPATH != null) {
            vwriter = new VCFWriter(vheader, VPATH);
        }
    }

    public ArrayList<VCFRecord> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ArrayList<VCFRecord> records = new ArrayList<VCFRecord>();

        for (ReferenceOrderedDatum rod : tracker.getAllRods()) {
            if (rod instanceof RodVCF) {
                RodVCF vcfrod = (RodVCF) rod;
                VCFRecord record = vcfrod.mCurrentRecord;

                if (SAMPLES == null) {
                    SAMPLES = new HashSet<String>();
                    SAMPLES.addAll(vcfrod.getHeader().getGenotypeSamples());
                }

                if (VPATH != null && vwriter == null) {
                    initializeWriter();
                }

                //out.println(record.toStringEncoding(vcfrod.getHeader()));

                records.add(record);
            }
        }
        
        return records;
    }

    public VCFWriter reduceInit() {
        return vwriter;
    }

    private VCFRecord subsetRecord(VCFRecord record) {
        ArrayList<VCFGenotypeRecord> genotypeRecords = new ArrayList<VCFGenotypeRecord>();
        for (int i = 0; i < record.getGenotypes().size(); i++) {
            VCFGenotypeRecord gr = (VCFGenotypeRecord)record.getGenotypes().get(i);

            //if (gr.getSampleName().equalsIgnoreCase(SAMPLE)) {
            if (SAMPLES.contains(gr.getSampleName())) {
                //out.println(gr.getSampleName() + " " + gr.toGenotypeString(record.getAlternateAlleles()));
                genotypeRecords.add(gr);
            }
        }

        VCFRecord subset = new VCFRecord(record.getReferenceBase(),
                                         record.getLocation().getContig(),
                                         (int) record.getLocation().getStart(),
                                         record.getID(),
                                         record.getAlternateAlleles(),
                                         record.getQual(),
                                         record.getFilterString(),
                                         record.getInfoValues(),
                                         record.getGenotypeFormatString(),
                                         genotypeRecords);

        return subset;
    }

    public VCFWriter reduce(ArrayList<VCFRecord> records, VCFWriter writer) {
        for (VCFRecord record : records) {
            VCFRecord subset = subsetRecord(record);

            boolean isVariant = false;
            for (VCFGenotypeEncoding ge : ((VCFGenotypeRecord)subset.getGenotypes().get(0)).getAlleles()) {
                if (record.getReferenceBase() != ge.getBases().charAt(0)) {
                    isVariant = true;
                }
            }

            //if (isVariant && !subset.isFiltered()) {
            if ((isVariant || INCLUDE_NON_VARIANTS) && (!subset.isFiltered() || INCLUDE_FILTERED)) {
                if (writer != null) {
                    writer.addRecord(subset);
                } else {
                    out.println(subset.toStringEncoding(vheader));
                }
            }
        }

        return writer;
    }

    public void onTraversalDone(VCFWriter writer) {
        if (writer != null) {
            writer.close();
        }
    }
}
