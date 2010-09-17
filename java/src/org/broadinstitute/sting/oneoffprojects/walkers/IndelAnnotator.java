package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.FeatureSource;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqCodec;
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqFeature;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.FeatureToGATKFeatureIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class IndelAnnotator extends RodWalker<Integer,Long> {
    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="refseq", shortName="refseq", doc="Name of RefSeq transcript annotation file", required=true)
	String RefseqFileName = null;

    private SeekableRODIterator refseqIterator;

    public void initialize() {
        if ( RefseqFileName != null ) {
            RMDTrackBuilder builder = new RMDTrackBuilder();
            FeatureSource refseq = builder.createFeatureReader(RefSeqCodec.class,new File(RefseqFileName)).first;

            try {
                refseqIterator = new SeekableRODIterator(new FeatureToGATKFeatureIterator(refseq.iterator(),"refseq"));
            } catch (IOException e) {
                throw new UserException.CouldNotReadInputFile(RefseqFileName, e);
            }

            logger.info("Using RefSeq annotations from " + RefseqFileName);
        }

        if ( refseqIterator == null ) logger.info("No annotations available");

        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "IndelAnnotator"));
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));

        HashSet<VCFInfoHeaderLine> anno = new HashSet<VCFInfoHeaderLine>();
        anno.add(new VCFInfoHeaderLine("cDNAchange", 1, VCFHeaderLineType.String, "cDNAchange"));
        anno.add(new VCFInfoHeaderLine("classification", 1, VCFHeaderLineType.String, "classification"));
        anno.add(new VCFInfoHeaderLine("codonchange", 1, VCFHeaderLineType.String, "codonchange"));
        anno.add(new VCFInfoHeaderLine("gene", 1, VCFHeaderLineType.String, "gene"));
        anno.add(new VCFInfoHeaderLine("genomechange", 1, VCFHeaderLineType.String, "genomechange"));
        anno.add(new VCFInfoHeaderLine("proteinchange", 1, VCFHeaderLineType.String, "proteinchange"));
        anno.add(new VCFInfoHeaderLine("strand", 1, VCFHeaderLineType.String, "strand"));
        anno.add(new VCFInfoHeaderLine("transcript", 1, VCFHeaderLineType.String, "transcript"));
        anno.add(new VCFInfoHeaderLine("type", 1, VCFHeaderLineType.String, "type"));
        hInfo.addAll(anno);

        VCFHeader vcfHeader = new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit()));
        vcfWriter.writeHeader(vcfHeader);
    }

    public Long reduceInit() {
        return 0l;
    }

    private TreeMap<String, Object> getAnnotationMap(VariantContext vc, VariantContext dbsnp, RODRecordList ann) {
        TreeMap<String, Object> anns = new TreeMap<String, Object>();
        anns.put("gene", "---");
        anns.put("type", "IGR");
        anns.put("transcript", "---");
        anns.put("strand", "+");
        anns.put("proteinchange", "---");
        anns.put("genomechange", "---");
        anns.put("codonchange", "---");
        anns.put("cDNAchange", "---");

        if (dbsnp != null) {
            anns.put("ID", dbsnp.getAttribute("ID"));
        }

        if (vc.isIndel()) {
            anns.put("classification", vc.isInsertion() ? "INS" : "DEL");
        }

        if ( ann != null ) {
            TreeMap<Integer, TreeMap<String, Object>> deleteriousnessRankedAnnotations = new TreeMap<Integer, TreeMap<String, Object>>();

            for (int transcriptIndex = 0; transcriptIndex < ann.size(); transcriptIndex++) {
                Transcript t = (Transcript) ann.get(transcriptIndex).getUnderlyingObject();
                TreeMap<String, Object> plausibleAnnotations = new TreeMap<String, Object>();
                Integer rank = 0;

                plausibleAnnotations.put("gene", t.getGeneName());
                plausibleAnnotations.put("transcript", t.getTranscriptId());
                plausibleAnnotations.put("strand", t.getStrand() == -1 ? "-" : "+");

                if ( RefSeqFeature.isExon(ann) ) {
                    if ( RefSeqFeature.isCodingExon(ann) ) {
                        //b.append(annCoding); // both exon and coding = coding exon sequence
                        if (vc.getIndelLengths().get(0) % 3 == 0) {
                            plausibleAnnotations.put("type", "Non-frameshift");
                            rank = 4;
                        } else {
                            plausibleAnnotations.put("type", "Frameshift");
                            rank = 0;
                        }
                    } else {
                        //b.append(annUTR); // exon but not coding = UTR
                        if (t.getStrand() == 1) {
                            plausibleAnnotations.put("type", "5'-UTR");
                            rank = 2;
                        } else {
                            plausibleAnnotations.put("type", "3'-UTR");
                            rank = 3;
                        }
                    }
                } else {
                    if ( RefSeqFeature.isCoding(ann) ) {
                        //b.append(annIntron); // not in exon, but within the coding region = intron
                        GenomeLoc ig = GenomeLocParser.createGenomeLoc(vc.getChr(), vc.getStart(), vc.getEnd());
                        GenomeLoc cl = t.getCodingLocation();
                        GenomeLoc g = t.getLocation();

                        boolean spliceSiteDisruption = false;

                        for (GenomeLoc exon : t.getExons()) {
                            GenomeLoc expandedExon = GenomeLocParser.createGenomeLoc(exon.getContig(), exon.getStart() - 6, exon.getStop() + 6);

                            if (ig.overlapsP(expandedExon)) {
                                spliceSiteDisruption = true;
                            }
                        }

                        if (spliceSiteDisruption) {
                            plausibleAnnotations.put("type", "SpliceSiteDisruption");
                            rank = 1;
                        } else {
                            plausibleAnnotations.put("type", "Intron");
                            rank = 5;
                        }
                    } else {
                        //b.append(annUnknown); // we have no idea what this is. this may actually happen when we have a fully non-coding exon...
                        plausibleAnnotations.put("type", "Unknown");
                        rank = 6;
                    }
                }

                deleteriousnessRankedAnnotations.put(rank, plausibleAnnotations);
            }

            anns.putAll(deleteriousnessRankedAnnotations.firstEntry().getValue());
        }
        

        return anns;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext con) {
        if ( tracker == null ) { return 0; }

        VariantContext vc = tracker.getVariantContext(ref, "variant", null, con.getLocation(), true);
        if ( vc == null ) { return 0; }

        Collection<VariantContext> dbsnps = tracker.getVariantContexts(ref, "dbsnp", null, con.getLocation(), true, true);
        VariantContext dbsnp = null;

        if (dbsnps != null && dbsnps.size() > 0) {
        ArrayList<VariantContext> dbsnpsarray = new ArrayList<VariantContext>(dbsnps);
            dbsnp = dbsnpsarray.get(0);
        }
        
        RODRecordList annotationList = refseqIterator.seekForward(ref.getLocus());
        TreeMap<String, Object> annotationMap = getAnnotationMap(vc, dbsnp, annotationList);

        Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());
        attrs.putAll(annotationMap);

        vc = VariantContext.modifyAttributes(vc, attrs);
        vcfWriter.add(vc, ref.getBase());

        return 1;
    }

    public Long reduce(Integer i, Long j) {
        return i + j;
    }
}