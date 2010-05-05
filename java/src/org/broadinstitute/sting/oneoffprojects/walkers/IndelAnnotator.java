package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeatureIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.util.*;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Apr 21, 2010
 */
public class IndelAnnotator extends RodWalker<Integer,Long>{
    @Argument(fullName="refseq", shortName="refseq",
			doc="Name of RefSeq transcript annotation file. If specified, indels will be annotated with GENOMIC/UTR/INTRON/CODING and with the gene name", required=true)
	String RefseqFileName = null;

    private static String annGenomic = "GENOMIC";
    private static String annIntron = "INTRON";
    private static String annUTR = "UTR";
    private static String annCoding = "CODING";
    private static String annUnknown = "UNKNOWN";

    private SeekableRODIterator refseqIterator;
    private VCFWriter vcfWriter;

    private String getAnnotationString(RODRecordList ann) {
        if ( ann == null ) return annGenomic;
        else {
            StringBuilder b = new StringBuilder();

            if ( rodRefSeq.isExon(ann) ) {
                if ( rodRefSeq.isCodingExon(ann) ) b.append(annCoding); // both exon and coding = coding exon sequence
                else b.append(annUTR); // exon but not coding = UTR
            } else {
                if ( rodRefSeq.isCoding(ann) ) b.append(annIntron); // not in exon, but within the coding region = intron
                else b.append(annUnknown); // we have no idea what this is. this may actually happen when we have a fully non-coding exon...
            }
            b.append('\t');
            b.append(((Transcript)ann.get(0).getUnderlyingObject()).getGeneName()); // there is at least one transcript in the list, guaranteed
            return b.toString();
        }
    }

    public void initialize() {
        if ( RefseqFileName != null ) {
            ReferenceOrderedData<rodRefSeq> refseq = new ReferenceOrderedData<rodRefSeq>("refseq",
                    new java.io.File(RefseqFileName), rodRefSeq.class);

            refseqIterator = new SeekableRODIterator(new GATKFeatureIterator(refseq.iterator()));
            logger.info("Using RefSeq annotations from "+RefseqFileName);
        }

        if ( refseqIterator == null ) logger.info("No annotations available");

    Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "IndelAnnotator"));
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));
        HashSet<VCFInfoHeaderLine> anno = new HashSet<VCFInfoHeaderLine>();
        anno.add(new VCFInfoHeaderLine("type",1,VCFInfoHeaderLine.INFO_TYPE.String,"Genomic interpretation (according to RefSeq)"));
        hInfo.addAll(anno);

        vcfWriter = new VCFWriter(out);
        VCFHeader vcfHeader = new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit()));
        vcfWriter.writeHeader(vcfHeader);
    }

    public Long reduceInit() {
        return 0l;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext con) {
        if ( tracker == null )
            return 0;

        List<Object> rods = tracker.getReferenceMetaData("variant");
        // ignore places where we don't have a variant
        if ( rods.size() == 0 )
            return 0;

        Object variant = rods.get(0);

        if ( variant instanceof VCFRecord) {
            RODRecordList annotationList = (refseqIterator == null ? null : refseqIterator.seekForward(ref.getLocus()));
            String annotationString = (refseqIterator == null ? "" : getAnnotationString(annotationList));
            annotationString = annotationString.split("\\s+")[0];
            ((VCFRecord) variant).addInfoField("type",annotationString);
            vcfWriter.addRecord((VCFRecord) variant);
        } else {
            throw new StingException("This one-off walker only deals with VCF files.");
        }

        return 1;
    }

    public Long reduce(Integer i, Long j) {
        return i + j;
    }

    public void onTraversalDone(Long l) {
        return;
    }

}