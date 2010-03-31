package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * Extracts subsets of a VCF file like one or more samples, all or only variant loci, all or filtered loci.
 */
public class VCFSubsetWalker extends RodWalker<ArrayList<VCFRecord>, VCFWriter> {
    @Argument(fullName="sample", shortName="SN", doc="Sample to include (or nothing to specify all samples)", required=false)
    private HashSet<String> SAMPLES;

    @Argument(fullName="vcfsubset", shortName="O", doc="File to write VCF subset to", required=false)
    private File VPATH = null;

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
            vwriter = new VCFWriter(VPATH);
            vwriter.writeHeader(vheader);
        }
    }

    public ArrayList<VCFRecord> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ArrayList<VCFRecord> records = new ArrayList<VCFRecord>();

        if (tracker != null) {
            for (GATKFeature feature : tracker.getAllRods()) {
                Object rod = feature.getUnderlyingObject();
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
        }
        
        return records;
    }

    public VCFWriter reduceInit() {
        return vwriter;
    }

    private VCFRecord subsetRecord(VCFRecord record) {
        ArrayList<VCFGenotypeRecord> genotypeRecords = new ArrayList<VCFGenotypeRecord>();
        HashSet<VCFGenotypeEncoding> genotypeEncodingSet = new HashSet<VCFGenotypeEncoding>();

        for ( VCFGenotypeRecord gr : record.getVCFGenotypeRecords() ) {
            if (SAMPLES.contains(gr.getSampleName())) {
                genotypeRecords.add(gr);

                for (VCFGenotypeEncoding allele : gr.getAlleles()) {
                    if (!allele.getBases().equalsIgnoreCase(record.getReference())) {
                        genotypeEncodingSet.add(allele);
                    }
                }
            }
        }

        ArrayList<VCFGenotypeEncoding> genotypeEncodings = new ArrayList<VCFGenotypeEncoding>();
        for (VCFGenotypeEncoding allele : genotypeEncodingSet) {
            genotypeEncodings.add(allele);
        }

        VCFRecord subset = new VCFRecord(record.getReference(),
                                         record.getLocation().getContig(),
                                         (int) record.getLocation().getStart(),
                                         record.getID(),
                                         genotypeEncodings,
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
            for ( VCFGenotypeEncoding ge : subset.getVCFGenotypeRecords().get(0).getAlleles() ) {
                if (!record.getReference().equals(ge.getBases())) {
                    isVariant = true;
                }
            }

            //if (isVariant && !subset.isFiltered()) {
            if ((isVariant || INCLUDE_NON_VARIANTS) && (!subset.isFiltered() || INCLUDE_FILTERED)) {
                if (vwriter != null) {
                    vwriter.addRecord(subset);
                } else {
                    out.println(subset.toStringEncoding(vheader));
                }
            }
        }

        return writer;
    }

    public void onTraversalDone(VCFWriter writer) {
        if (vwriter != null) {
            vwriter.close();
        }
    }
}
