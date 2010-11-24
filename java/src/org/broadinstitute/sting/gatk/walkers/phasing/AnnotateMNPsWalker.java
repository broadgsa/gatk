/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator.AminoAcid;
import org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator.AminoAcidTable;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

import static org.broadinstitute.sting.utils.vcf.VCFUtils.getVCFHeadersFromRods;


/**
 * Walks along all variant ROD loci, and dynamically annotates alleles at MNP records.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = AnnotateMNPsWalker.REFSEQ_ROD_NAME, type = AnnotatorInputTableFeature.class), @RMD(name = AnnotateMNPsWalker.VARIANT_ROD_NAME, type = ReferenceOrderedDatum.class)})

public class AnnotateMNPsWalker extends RodWalker<Integer, Integer> {

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter writer = null;
    private ManualSortingVCFWriter sortingWriter = null;

    private LinkedList<String> rodNames = null;
    private GenomeLocParser locParser = null;
    private TreeMap<GenomeLoc, Set<GenomeLoc>> MNPstartToStops = null; // Must be TreeMap sorted by START sites!

    public final static String REFSEQ_ROD_NAME = "refseq";
    public final static String VARIANT_ROD_NAME = "variant";

    private LocusToFeatures locusToRefSeqFeatures = null;


    protected final static String MNP_ANNOTATION_KEY_PREFIX = "MNP.refseq.";

    protected final static String REFSEQ_NAME = "name";
    protected final static String REFSEQ_NAME2 = "name2";

    protected final static String REFSEQ_POSITION_TYPE = "positionType";
    protected final static String REFSEQ_CDS = "CDS";

    protected final static String REFSEQ_STRAND = "transcriptStrand";
    protected final static String REFSEQ_POS_STRAND = "+";
    protected final static String REFSEQ_NEG_STRAND = "-";

    protected final static String REFSEQ_CODON_COORD = "codonCoord";
    protected final static String REFSEQ_CODING_FRAME = "frame";

    protected final static String REFSEQ_REF_CODON = "referenceCodon";
    protected final static String REFSEQ_REF_AA = "referenceAA";

    protected final static String REFSEQ_ALT_BASE = "haplotypeAlternate";

    protected final static String REFSEQ_VARIANT_CODON = "variantCodon";
    protected final static String REFSEQ_VARIANT_AA = "variantAA";
    protected final static String REFSEQ_CHANGES_AA = "changesAA";
    protected final static String REFSEQ_FUNCTIONAL_CLASS = "functionalClass";
    protected final static String REFSEQ_PROTEIN_COORD_DESCRIPTION = "proteinCoordStr";

    protected final static String REFSEQ_CODING_ANNOTATIONS = "codingVariants";

    public void initialize() {
        rodNames = new LinkedList<String>();
        rodNames.add(VARIANT_ROD_NAME);

        locParser = getToolkit().getGenomeLocParser();
        MNPstartToStops = new TreeMap<GenomeLoc, Set<GenomeLoc>>(); // sorted by start sites

        initializeVcfWriter();

        locusToRefSeqFeatures = new LocusToFeatures();
    }

    private void initializeVcfWriter() {
        sortingWriter = new ManualSortingVCFWriter(writer);
        writer = sortingWriter;

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), rodNames);
        writer.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(rodNameToHeader.get(rodNames.get(0)).getGenotypeSamples())));
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * For each site of interest, annotate it if it's a MNP.
     *
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return count of MNPs observed
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        int numMNPsObserved = 0;
        GenomeLoc curLocus = ref.getLocus();
        clearOldLocusFeatures(curLocus);

        boolean requireStartHere = false; // see EVERY site of the MNP
        boolean takeFirstOnly = false; // take as many entries as the VCF file has
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly)) {
            if (vc.getType() == VariantContext.Type.MNP) {
                GenomeLoc MNPloc = VariantContextUtils.getLocation(locParser, vc);
                boolean atStartOfMNP = curLocus.getStart() == MNPloc.getStart();
                boolean atEndOfMNP = curLocus.getStart() == MNPloc.getStop();
                logger.debug("Observed MNP at " + MNPloc);

                if (isChrM(vc)) {
                    if (atStartOfMNP) {
                        logger.warn("Skipping mitochondrial MNP at " + MNPloc + " due to complexity of coding table [need to know if first codon, etc.]...");
                        writeVCF(vc);
                    }
                    continue;
                }

                GenomeLoc stopLoc = locParser.createGenomeLoc(curLocus.getContig(), MNPloc.getStop());
                final List<Object> refSeqRODs = tracker.getReferenceMetaData(REFSEQ_ROD_NAME);
                for (Object refSeqObject : refSeqRODs) {
                    AnnotatorInputTableFeature refSeqAnnotation = (AnnotatorInputTableFeature) refSeqObject;
                    locusToRefSeqFeatures.putLocusFeatures(curLocus, refSeqAnnotation, stopLoc);
                }

                if (atStartOfMNP) { // MNP is starting here, so register that we're waiting for it
                    Set<GenomeLoc> stopLocs = MNPstartToStops.get(curLocus);
                    if (stopLocs == null) {
                        stopLocs = new HashSet<GenomeLoc>();
                        MNPstartToStops.put(curLocus, stopLocs);
                    }
                    stopLocs.add(stopLoc);
                }

                if (atEndOfMNP) {
                    numMNPsObserved++; // only count a MNP at its stop site
                    logger.debug("Observed end of MNP at " + curLocus);
                    logger.debug("Current list of per-locus features\n" + locusToRefSeqFeatures);

                    Map<String, Object> MNPannotations = annotateMNP(vc);
                    MNPannotations.putAll(vc.getAttributes());
                    vc = VariantContext.modifyAttributes(vc, MNPannotations);
                    writeVCF(vc);

                    GenomeLoc startLoc = locParser.createGenomeLoc(curLocus.getContig(), MNPloc.getStart());
                    Set<GenomeLoc> stopLocs = MNPstartToStops.get(startLoc);
                    if (stopLocs != null) { // otherwise, just removed stopLocs due to another MNP that has the same (start, stop)
                        stopLocs.remove(stopLoc);
                        if (stopLocs.isEmpty()) // no longer waiting for startLoc
                            MNPstartToStops.remove(startLoc);
                    }
                }
            }
            else {
                writeVCF(vc);
            }
        }

        Integer mostUpstreamWritableLoc = null;
        if (!MNPstartToStops.isEmpty()) {
            GenomeLoc waitingForLoc = MNPstartToStops.entrySet().iterator().next().getKey();
            mostUpstreamWritableLoc = waitingForLoc.getStart() - 1;
        }
        sortingWriter.setmostUpstreamWritableLocus(mostUpstreamWritableLoc);

        return numMNPsObserved;
    }

    private static boolean isChrM(final VariantContext vc) {
        return vc.getChr().equals("chrM") || vc.getChr().equals("MT");
    }

    private Map<String, Object> annotateMNP(VariantContext vc) {
        Map<String, Object> annotations = new HashMap<String, Object>();

        RefSeqNameToFeatures nameToPositionalFeatures = new RefSeqNameToFeatures(vc);
        MNPannotationKeyBuilder kb = new MNPannotationKeyBuilder(nameToPositionalFeatures);

        for (Map.Entry<String, RefSeqFeatureList> nameToFeatureEntry : nameToPositionalFeatures.entrySet()) {
            String featureName = nameToFeatureEntry.getKey();
            RefSeqFeatureList feature = nameToFeatureEntry.getValue();
            CodonAnnotationsForAltAlleles codonAnnotationsForAlleles = new CodonAnnotationsForAltAlleles(vc, feature);

            annotations.put(kb.getKey(REFSEQ_CODING_ANNOTATIONS), codonAnnotationsForAlleles.toString());
            annotations.put(kb.getKey(REFSEQ_NAME), featureName);
            annotations.put(kb.getKey(REFSEQ_NAME2), feature.name2);
            annotations.put(kb.getKey(REFSEQ_POSITION_TYPE), REFSEQ_CDS);
            annotations.put(kb.getKey(REFSEQ_STRAND), (feature.positiveStrand ? REFSEQ_POS_STRAND : REFSEQ_NEG_STRAND));
            annotations.put(kb.getKey(REFSEQ_CODON_COORD), feature.getCodonCoordString());

            kb.incrementFeatureIndex();
        }

        return annotations;
    }

    private static class MNPannotationKeyBuilder {
        private int featureIndex;
        private boolean multipleEntries;

        public MNPannotationKeyBuilder(RefSeqNameToFeatures nameToPositionalFeatures) {
            this.featureIndex = 1;
            this.multipleEntries = nameToPositionalFeatures.nameToFeatures.size() > 1;
        }

        public void incrementFeatureIndex() {
            featureIndex++;
        }

        public String getKey(String type) {
            String annotationKey = MNP_ANNOTATION_KEY_PREFIX + type;
            if (multipleEntries)
                annotationKey += "_" + featureIndex;
            return annotationKey;
        }
    }

    private static byte[] ByteArrayToPrimitive(Byte[] nonNullArray) {
        byte[] primArray = new byte[nonNullArray.length];

        for (int i = 0; i < nonNullArray.length; i++) {
            if (nonNullArray[i] == null)
                throw new ReviewedStingException("nonNullArray[i] == null");
            primArray[i] = nonNullArray[i];
        }

        return primArray;
    }

    private void clearOldLocusFeatures(GenomeLoc curLoc) {
        Iterator<Map.Entry<GenomeLoc, PositionalRefSeqFeatures>> locusFeaturesIt = locusToRefSeqFeatures.entrySet().iterator();
        while (locusFeaturesIt.hasNext()) {
            Map.Entry<GenomeLoc, PositionalRefSeqFeatures> locusFeaturesEntry = locusFeaturesIt.next();
            if (curLoc.isPast(locusFeaturesEntry.getValue().getFurthestLocusUsingFeatures()))
                locusFeaturesIt.remove();
        }
    }

    public Integer reduce(Integer count, Integer total) {
        if (count != null)
            total = total + count;

        return total;
    }

    /**
     * @param result the number of MNPs processed.
     */
    public void onTraversalDone(Integer result) {
        System.out.println("Number of MNPs observed: " + result);
        writer.close();
    }

    private void writeVCF(VariantContext vc) {
        byte refBase;
        if (!vc.isIndel()) {
            Allele varAllele = vc.getReference();
            refBase = SNPallelePair.getSingleBase(varAllele);
        }
        else {
            refBase = vc.getReferenceBaseForIndel();
        }

        writer.add(vc, refBase);
    }

    /*
     Inner classes:
     */

    // Maps: RefSeq entry name -> features for ALL positions of a particular VariantContext MNP:

    private class RefSeqNameToFeatures {
        private Map<String, RefSeqFeatureList> nameToFeatures;

        public RefSeqNameToFeatures(VariantContext vc) {
            this.nameToFeatures = new HashMap<String, RefSeqFeatureList>();

            int MNPstart = vc.getStart();
            int MNPstop = vc.getEnd();
            int MNPlength = MNPstop - MNPstart + 1;

            for (int i = 0; i < MNPlength; i++) {
                int genomicPosition = MNPstart + i;
                GenomeLoc posLoc = locParser.createGenomeLoc(vc.getChr(), genomicPosition);

                PositionalRefSeqFeatures locFeatures = locusToRefSeqFeatures.getLocusFeatures(posLoc);
                if (locFeatures == null) // no features for posLoc
                    continue;

                for (Map.Entry<String, PositionalRefSeqFeature> nameToFeatureEntry : locFeatures.entrySet()) {
                    String name = nameToFeatureEntry.getKey();
                    PositionalRefSeqFeature posFeature = nameToFeatureEntry.getValue();

                    RefSeqFeatureList featureList = nameToFeatures.get(name);
                    if (featureList == null) {
                        featureList = new RefSeqFeatureList(MNPlength);
                        nameToFeatures.put(name, featureList);
                    }
                    featureList.updateFeatureAtPosition(i, posFeature);
                }
            }
        }

        public Set<Map.Entry<String, RefSeqFeatureList>> entrySet() {
            return nameToFeatures.entrySet();
        }
    }

    // For a particular RefSeq entry, contains the features for ALL positions of a particular VariantContext MNP

    private static class RefSeqFeatureList {
        private final static String CODON_FRAME_START = "(";
        private final static String CODON_FRAME_END = ")";
        private final static String CODON_DELIM = "|";

        private CodingRefSeqFeature[] refSeqFeatures;
        private String name2;
        private Boolean positiveStrand;

        private Map<Integer, List<Integer>> codonToIndices; // Map of: codon index -> MNP indices that refer to codon

        public RefSeqFeatureList(int MNPlength) {
            this.refSeqFeatures = new CodingRefSeqFeature[MNPlength];
            for (int i = 0; i < MNPlength; i++)
                this.refSeqFeatures[i] = null;

            this.name2 = null;
            this.positiveStrand = null;
            this.codonToIndices = new TreeMap<Integer, List<Integer>>();
        }

        public void updateFeatureAtPosition(int index, PositionalRefSeqFeature feature) {
            if (name2 == null) {
                name2 = feature.name2;
                positiveStrand = feature.positiveStrand;
            }
            else if (!name2.equals(feature.name2) || positiveStrand != feature.positiveStrand) {
                throw new UserException("Inconsistency between previous RefSeq entry and: " + feature);
            }

            CodingRefSeqFeature crsf = new CodingRefSeqFeature(feature);
            refSeqFeatures[index] = crsf;

            List<Integer> indicesWithCodon = codonToIndices.get(crsf.codonCoord);
            if (indicesWithCodon == null) {
                indicesWithCodon = new LinkedList<Integer>();
                codonToIndices.put(crsf.codonCoord, indicesWithCodon);
            }
            indicesWithCodon.add(index);
        }

        public Set<Map.Entry<Integer, List<Integer>>> codonIndicesEntrySet() {
            return codonToIndices.entrySet();
        }

        public String getCodonCoordString() {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < refSeqFeatures.length; i++) {
                CodingRefSeqFeature crsf = refSeqFeatures[i];
                if (crsf != null)
                    sb.append(crsf.codonCoord).append(CODON_FRAME_START).append(crsf.codingFrame).append(CODON_FRAME_END);
                if (i < refSeqFeatures.length - 1)
                    sb.append(CODON_DELIM);
            }

            return sb.toString();
        }
    }

    private static class CodingRefSeqFeature {
        protected int codonCoord;
        protected int codingFrame;
        protected String referenceCodon;
        protected String referenceAA;

        public CodingRefSeqFeature(PositionalRefSeqFeature feature) {
            this.codonCoord = feature.codonCoord;
            this.codingFrame = feature.codingFrame;
            this.referenceCodon = feature.referenceCodon.toUpperCase();
            this.referenceAA = feature.referenceAA;
        }
    }

    private static class CodonAnnotationsForAltAlleles {
        protected final static int MIN_CODON_INDEX = 0;
        protected final static int NUM_CODON_INDICES = 3;
        private final static String CODON_ANNOTATION_DELIM = ",";

        private List<SingleCodonAnnotationsForAlleles> alleleAnnotations;

        public CodonAnnotationsForAltAlleles(VariantContext vc, RefSeqFeatureList feature) {
            this.alleleAnnotations = new LinkedList<SingleCodonAnnotationsForAlleles>();

            int MNPstart = vc.getStart();
            int MNPstop = vc.getEnd();
            int MNPlength = MNPstop - MNPstart + 1;

            for (Map.Entry<Integer, List<Integer>> codonToIndicesEntry : feature.codonIndicesEntrySet()) {
                int codonIndex = codonToIndicesEntry.getKey();
                List<Integer> indices = codonToIndicesEntry.getValue();
                if (indices.isEmpty())
                    throw new ReviewedStingException("indices should not exist if it's empty!");

                for (int index : indices) {
                    int frame = feature.refSeqFeatures[index].codingFrame;
                    if (feature.refSeqFeatures[index].codonCoord != codonIndex)
                        throw new ReviewedStingException("LOGICAL ERROR: feature.refSeqFeatures[index].codonCoord != codonIndex");
                    if (frame < MIN_CODON_INDEX || frame >= NUM_CODON_INDICES)
                        throw new UserException("RefSeq codon frame not one of {0,1,2}");
                }
                CodingRefSeqFeature firstFeatureForCodon = feature.refSeqFeatures[indices.get(0)];
                String refCodon = firstFeatureForCodon.referenceCodon;

                SingleCodonAnnotationsForAlleles codonAnnotation = new SingleCodonAnnotationsForAlleles(codonIndex, vc.getAlternateAlleles(), MNPlength, refCodon, firstFeatureForCodon, indices, feature);
                alleleAnnotations.add(codonAnnotation);
            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            int index = 0;
            for (SingleCodonAnnotationsForAlleles codonToAlleles : alleleAnnotations) {
                sb.append(codonToAlleles);
                if (index < alleleAnnotations.size() - 1)
                    sb.append(CODON_ANNOTATION_DELIM);
                index++;
            }

            return sb.toString();
        }
    }

    private static class SingleCodonAnnotationsForAlleles {
        private final static String CODON_MAP_SYMBOL = "->";
        private final static String CODON_ANNOTATION_START = "[";
        private final static String CODON_ANNOTATION_END = "]";
        private final static String REF_CODON_INFO_DELIM = "|";
        private final static String ALLELE_ANNOTATION_DELIM = ",";
        private final static String ASSIGNMENT = ":";

        private int codonIndex;
        private String refCodon;
        private String refAA;

        private List<SingleCodonAnnotationsForAllele> annotationsForAlleles;

        public SingleCodonAnnotationsForAlleles(int codonIndex, Collection<Allele> altAlleles, int MNPlength, String refCodon, CodingRefSeqFeature firstFeatureForCodon, List<Integer> indices, RefSeqFeatureList feature) {
            if (refCodon.length() != CodonAnnotationsForAltAlleles.NUM_CODON_INDICES)
                throw new UserException("RefSeq reference codon " + refCodon + " is not of length " + CodonAnnotationsForAltAlleles.NUM_CODON_INDICES);

            AminoAcid refAA = AminoAcidTable.getEukaryoticAA(refCodon);
            if (!refAA.getCode().equals(firstFeatureForCodon.referenceAA))
                throw new UserException("RefSeq: translated reference codon= " + refAA + " != " + firstFeatureForCodon.referenceAA + " = reference AA");

            this.codonIndex = codonIndex;
            this.refCodon = refCodon;
            this.refAA = refAA.getCode();
            this.annotationsForAlleles = new LinkedList<SingleCodonAnnotationsForAllele>();

            for (Allele altAllele : altAlleles) {
                if (altAllele.length() != MNPlength)
                    throw new ReviewedStingException("length(altAllele) != length(MNP)");
                byte[] altBases = altAllele.getBases();

                Byte[] variantCodonArr = new Byte[CodonAnnotationsForAltAlleles.NUM_CODON_INDICES];
                for (int i = CodonAnnotationsForAltAlleles.MIN_CODON_INDEX; i < CodonAnnotationsForAltAlleles.NUM_CODON_INDICES; i++)
                    variantCodonArr[i] = null;

                for (int index : indices) {
                    int frame = feature.refSeqFeatures[index].codingFrame;
                    if (variantCodonArr[frame] != null)
                        throw new UserException("RefSeq assigns codon " + codonIndex + " twice at same frame: " + frame);

                    byte base = altBases[index];
                    if (!feature.positiveStrand) // negative strand codon
                        base = BaseUtils.simpleComplement(base);

                    variantCodonArr[frame] = base;
                }

                /* For missing frames, there MUST exist AT LEAST one index that refers to this codon,
                  so use it to derive the missing bases [ALREADY complemented if on the negative strand]:
                */
                for (int frame = CodonAnnotationsForAltAlleles.MIN_CODON_INDEX; frame < CodonAnnotationsForAltAlleles.NUM_CODON_INDICES; frame++) {
                    if (variantCodonArr[frame] == null)
                        variantCodonArr[frame] = (byte) refCodon.charAt(frame);
                }
                String variantCodon = new String(ByteArrayToPrimitive(variantCodonArr)).toUpperCase();

                SingleCodonAnnotationsForAllele alleleAnnotation = new SingleCodonAnnotationsForAllele(variantCodon, refCodon, refAA, codonIndex);
                annotationsForAlleles.add(alleleAnnotation);
            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append(codonIndex).append(CODON_MAP_SYMBOL).append(CODON_ANNOTATION_START);
            sb.append(REFSEQ_REF_CODON).append(ASSIGNMENT).append(refCodon).append(REF_CODON_INFO_DELIM);
            sb.append(REFSEQ_REF_AA).append(ASSIGNMENT).append(refAA).append(REF_CODON_INFO_DELIM);

            int index = 0;
            for (SingleCodonAnnotationsForAllele annotation : annotationsForAlleles) {
                sb.append(annotation);
                if (index < annotationsForAlleles.size() - 1)
                    sb.append(ALLELE_ANNOTATION_DELIM);
                index++;
            }
            sb.append(CODON_ANNOTATION_END);

            return sb.toString();
        }
    }

    private static class SingleCodonAnnotationsForAllele {
        private final static String ALLELE_START = "{";
        private final static String ALLELE_END = "}";
        private final static String CODON_INFO_DELIM = "|";
        private final static String ASSIGNMENT = ":";

        private String variantCodon;
        private String variantAA;
        private boolean changesAA;
        private String functionalClass;
        private String proteinCoordStr;

        public SingleCodonAnnotationsForAllele(String variantCodon, String refCodon, AminoAcid refAA, int codonIndex) {
            this.variantCodon = variantCodon;
            AminoAcid variantAA = AminoAcidTable.getEukaryoticAA(this.variantCodon);
            this.variantAA = variantAA.getCode();

            this.changesAA = !refAA.equals(variantAA);

            if (!this.variantCodon.equals(refCodon)) {
                if (changesAA) {
                    if (variantAA.isStop()) {
                        functionalClass = "nonsense";
                    }
                    else if (refAA.isStop()) {
                        functionalClass = "readthrough";
                    }
                    else {
                        functionalClass = "missense";
                    }
                }
                else { // the same aa:
                    functionalClass = "silent";
                }
            }
            else { // the same codon:
                functionalClass = "no_change";
            }

            this.proteinCoordStr = "p." + refAA.getLetter() + codonIndex + variantAA.getLetter();
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append(ALLELE_START);
            sb.append(REFSEQ_VARIANT_CODON).append(ASSIGNMENT).append(variantCodon).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_VARIANT_AA).append(ASSIGNMENT).append(variantAA).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_CHANGES_AA).append(ASSIGNMENT).append(changesAA).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_FUNCTIONAL_CLASS).append(ASSIGNMENT).append(functionalClass).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_PROTEIN_COORD_DESCRIPTION).append(ASSIGNMENT).append(proteinCoordStr);
            sb.append(ALLELE_END);

            return sb.toString();
        }
    }
}


// External classes:

class LocusToFeatures {
    private Map<GenomeLoc, PositionalRefSeqFeatures> locusToFeatures;

    public LocusToFeatures() {
        this.locusToFeatures = new TreeMap<GenomeLoc, PositionalRefSeqFeatures>();
    }

    public PositionalRefSeqFeatures getLocusFeatures(GenomeLoc loc) {
        return locusToFeatures.get(loc);
    }

    public void putLocusFeatures(GenomeLoc loc, AnnotatorInputTableFeature refSeqAnnotation, GenomeLoc locusUsingThis) {
        PositionalRefSeqFeatures locFeatures = locusToFeatures.get(loc);
        if (locFeatures == null) {
            locFeatures = new PositionalRefSeqFeatures(locusUsingThis);
            locusToFeatures.put(loc, locFeatures);
        }
        locFeatures.putFeature(refSeqAnnotation, locusUsingThis);
    }

    public Set<Map.Entry<GenomeLoc, PositionalRefSeqFeatures>> entrySet() {
        return locusToFeatures.entrySet();
    }

    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        for (Map.Entry<GenomeLoc, PositionalRefSeqFeatures> locFeatures : entrySet()) {
            GenomeLoc loc = locFeatures.getKey();
            PositionalRefSeqFeatures features = locFeatures.getValue();
            sb.append("Locus: ").append(loc).append("\n").append(features);
        }

        return sb.toString();
    }
}

class PositionalRefSeqFeatures {
    private final static String[] REQUIRE_COLUMNS =
            {AnnotateMNPsWalker.REFSEQ_NAME, AnnotateMNPsWalker.REFSEQ_POSITION_TYPE};

    private Map<String, PositionalRefSeqFeature> nameToFeature;
    private GenomeLoc furthestLocusUsingFeatures;

    public PositionalRefSeqFeatures(GenomeLoc locusUsingThis) {
        this.nameToFeature = new HashMap<String, PositionalRefSeqFeature>();
        this.furthestLocusUsingFeatures = locusUsingThis;
    }

    public void putFeature(AnnotatorInputTableFeature refSeqAnnotation, GenomeLoc locusUsingThis) {
        for (String column : REQUIRE_COLUMNS) {
            if (!refSeqAnnotation.containsColumnName(column))
                throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + column);
        }

        if (locusUsingThis.isPast(furthestLocusUsingFeatures))
            furthestLocusUsingFeatures = locusUsingThis;

        String posType = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_POSITION_TYPE);
        if (!posType.equals(AnnotateMNPsWalker.REFSEQ_CDS)) // only interested in coding sequence annotations
            return;

        PositionalRefSeqFeature newLocusFeature = new PositionalRefSeqFeature(refSeqAnnotation);

        String refSeqName = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_NAME);
        PositionalRefSeqFeature locusFeature = nameToFeature.get(refSeqName);
        if (locusFeature == null) {
            locusFeature = newLocusFeature;
            nameToFeature.put(refSeqName, locusFeature);
        }
        else if (!locusFeature.equals(newLocusFeature)) {
            throw new UserException("Inconsistency between previous RefSeq entry and: " + refSeqAnnotation);
        }

        locusFeature.updateFeature(refSeqAnnotation);
    }

    public GenomeLoc getFurthestLocusUsingFeatures() {
        return furthestLocusUsingFeatures;
    }

    public Set<Map.Entry<String, PositionalRefSeqFeature>> entrySet() {
        return nameToFeature.entrySet();
    }

    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        for (Map.Entry<String, PositionalRefSeqFeature> nameFeatureEntry : entrySet()) {
            String name = nameFeatureEntry.getKey();
            PositionalRefSeqFeature feature = nameFeatureEntry.getValue();
            sb.append(name).append(" -> [").append(feature).append("]\n");
        }

        return sb.toString();
    }
}

class PositionalRefSeqFeature {
    private final static String[] REQUIRE_COLUMNS =
            {AnnotateMNPsWalker.REFSEQ_NAME2, AnnotateMNPsWalker.REFSEQ_STRAND,
                    AnnotateMNPsWalker.REFSEQ_CODON_COORD, AnnotateMNPsWalker.REFSEQ_CODING_FRAME,
                    AnnotateMNPsWalker.REFSEQ_REF_CODON, AnnotateMNPsWalker.REFSEQ_REF_AA};

    protected String name2;
    protected boolean positiveStrand;
    protected int codonCoord;
    protected int codingFrame;
    protected String referenceCodon;
    protected String referenceAA;

    private Map<String, BaseAnnotations> baseToAnnotations;

    public PositionalRefSeqFeature(AnnotatorInputTableFeature refSeqAnnotation) {
        for (String column : REQUIRE_COLUMNS) {
            if (!refSeqAnnotation.containsColumnName(column))
                throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + column);
        }
        this.name2 = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_NAME2);
        this.positiveStrand = (refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_STRAND).equals(AnnotateMNPsWalker.REFSEQ_POS_STRAND));
        this.codonCoord = Integer.parseInt(refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_CODON_COORD));
        this.codingFrame = Integer.parseInt(refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_CODING_FRAME));
        this.referenceCodon = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_REF_CODON);
        this.referenceAA = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_REF_AA);

        this.baseToAnnotations = new HashMap<String, BaseAnnotations>();
    }

    public boolean equals(PositionalRefSeqFeature that) {
        return this.name2.equals(that.name2) && this.positiveStrand == that.positiveStrand && this.codonCoord == that.codonCoord && this.codingFrame == that.codingFrame
                && this.referenceCodon.equals(that.referenceCodon) && this.referenceAA.equals(that.referenceAA);
    }

    public void updateFeature(AnnotatorInputTableFeature refSeqAnnotation) {
        if (!refSeqAnnotation.containsColumnName(AnnotateMNPsWalker.REFSEQ_ALT_BASE))
            throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + AnnotateMNPsWalker.REFSEQ_ALT_BASE);
        String base = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_ALT_BASE);

        baseToAnnotations.put(base, new BaseAnnotations(refSeqAnnotation));
    }

    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        sb.append("name2= ").append(name2);
        sb.append(", positiveStrand= ").append(positiveStrand);
        sb.append(", codonCoord= ").append(codonCoord);
        sb.append(", codingFrame= ").append(codingFrame);
        sb.append(", referenceCodon= ").append(referenceCodon);
        sb.append(", referenceAA= ").append(referenceAA);

        sb.append(", baseAnnotations= {");
        for (Map.Entry<String, BaseAnnotations> baseToAnnotationsEntry : baseToAnnotations.entrySet()) {
            String base = baseToAnnotationsEntry.getKey();
            BaseAnnotations annotations = baseToAnnotationsEntry.getValue();
            sb.append(" ").append(base).append(" -> {").append(annotations).append("}");
        }
        sb.append(" }");

        return sb.toString();
    }
}

class BaseAnnotations {
    private final static String[] REQUIRE_COLUMNS =
            {AnnotateMNPsWalker.REFSEQ_VARIANT_CODON, AnnotateMNPsWalker.REFSEQ_VARIANT_AA,
                    AnnotateMNPsWalker.REFSEQ_CHANGES_AA, AnnotateMNPsWalker.REFSEQ_FUNCTIONAL_CLASS,
                    AnnotateMNPsWalker.REFSEQ_PROTEIN_COORD_DESCRIPTION};

    protected String variantCodon;
    protected String variantAA;
    protected boolean changesAA;
    protected String functionalClass;
    protected String proteinCoordStr;

    public BaseAnnotations(AnnotatorInputTableFeature refSeqAnnotation) {
        for (String column : REQUIRE_COLUMNS) {
            if (!refSeqAnnotation.containsColumnName(column))
                throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + column);
        }
        this.variantCodon = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_VARIANT_CODON);
        this.variantAA = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_VARIANT_AA);
        this.changesAA = Boolean.parseBoolean(refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_CHANGES_AA));
        this.functionalClass = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_FUNCTIONAL_CLASS);
        this.proteinCoordStr = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_PROTEIN_COORD_DESCRIPTION);
    }


    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        sb.append("variantCodon= ").append(variantCodon);
        sb.append(", variantAA= ").append(variantAA);
        sb.append(", changesAA= ").append(changesAA);
        sb.append(", functionalClass= ").append(functionalClass);
        sb.append(", proteinCoordStr= ").append(proteinCoordStr);

        return sb.toString();
    }
}