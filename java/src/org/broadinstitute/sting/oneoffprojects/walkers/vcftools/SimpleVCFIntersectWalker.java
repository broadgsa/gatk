package org.broadinstitute.sting.oneoffprojects.walkers.vcftools;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.ListUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.lang.Long;
import java.util.*;

/**
 * This walker intersects the samples for two VCF files and returns a set of calls determined
 * by command-line arguments. Intersect indicates that calls should be returned only where
 * the same locus appears in the two VCF files, while Union indicates that all calls should be
 * returned. Priority indicates which file to prefer when setting genotypes; if the "GQ" flag
 * is omitted, genotypes are merged according to highest genotype quality.
 *
 * @Author chartl
 * @Date Jan 29, 2010
 */
public class SimpleVCFIntersectWalker extends RodWalker<VCFRecordMerger,Long>{
    @Argument(fullName="locusUnion", shortName="union", doc="Will union loci rather than intersecting", required=false)
    boolean union = false;
    @Argument(fullName="mergeByQual", shortName="GQ", doc="Will merge genotypes according to quality rather than priority", required=false)
    boolean useQuality = false;

    private VCFRecordMerger merger; // no cost of repeated instantiation
    private VCFWriter writer;
    private boolean headerNotWritten;

    public void initialize() {
        merger = new VCFRecordMerger();
        writer = new VCFWriter(out);
        headerNotWritten = true;
    }

    public Long reduceInit() {
        return 0l;
    }

    public VCFRecordMerger map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return null;
        }

        RodVCF priorityCall = ( RodVCF ) tracker.lookup("priority",null);
        RodVCF otherCall = ( RodVCF ) tracker.lookup("other",null);

        if ( priorityCall == null && otherCall == null ) {
            return null;
        }

        if ( priorityCall == null ) {

            if ( ! union ) {
                return null;
            } else if ( merger.hasBeenInstantiated() ) {
                merger.update(null,otherCall.getRecord());
            } else {
                merger.hold(otherCall.getRecord());
            }

        } else if ( otherCall == null ) {

            if ( ! union ) {
                return null;
            } else if ( merger.hasBeenInstantiated() ) {
                merger.update(priorityCall.getRecord(),null);
            } else {
                merger.hold(priorityCall.getRecord());
            }

        } else {

            if ( merger.hasBeenInstantiated() ) {
                merger.update(priorityCall.getRecord(), otherCall.getRecord());
            } else {
                merger = instantiateMerger(priorityCall.getRecord(),otherCall.getRecord());
            }

        }

        return merger;
    }

    public Long reduce(VCFRecordMerger records, Long prevReduce) {
        if ( records == null ) {
            return prevReduce;
        } else {
            if ( records.hasBeenInstantiated() ) {
                if ( headerNotWritten ) {
                    writeHeader(records);
                }
                while ( records.hasHeldRecords() ) {
                    writer.addRecord(records.getHeldRecord());
                    prevReduce++;
                }

                writer.addRecord(records.getMergedRecord());
                prevReduce++;
            }
        }

        return prevReduce;
    }


    public void onTraversalDone(Long variantsAdded) {
        logger.info(variantsAdded);
    }

    private VCFRecordMerger instantiateMerger(VCFRecord priority, VCFRecord other) {
        List<String> samples = Arrays.asList(priority.getSampleNames());
        samples.retainAll(Arrays.asList(other.getSampleNames()));

        return new VCFRecordMerger(priority,other,useQuality,true,samples,merger.getAllHeldRecords());
    }

    private void writeHeader(VCFRecordMerger records) {
        VCFHeader header = new VCFHeader(new HashSet<VCFHeaderLine>(),records.getSamples());
        writer.writeHeader(header);
        headerNotWritten = false;
    }

}

class VCFRecordMerger {
    private List<VCFRecord> holder;
    private VCFGenotypeMerger genotypeMerger;
    private VCFRecord record1;
    private VCFRecord record2;
    private boolean isMergeByQual;
    private boolean record1Priority;
    private Collection<String> samples;

    public VCFRecordMerger(VCFRecord firstRecord, VCFRecord secondRecord, boolean useQual, boolean record1Priority, Collection<String> intersectedSamples, List<VCFRecord> heldRecords) {
        this.record1 = firstRecord;
        this.record2 = secondRecord;
        this.isMergeByQual = useQual;
        this.record1Priority = record1Priority;
        this.samples = intersectedSamples;
        instantiateGenotypeMerger();
        holder = heldRecords;
    }

    public Set<String> getSamples() {
        Set<String> sam = new HashSet<String>();
        for ( String s : samples ) {
            sam.add(s);
        }
        return sam;
    }

    public VCFRecordMerger() {
        holder = new ArrayList<VCFRecord>();
    }

    public List<VCFRecord> getAllHeldRecords() {
        return holder;
    }

    public void update(VCFRecord firstRecord, VCFRecord secondRecord) {
        record1 = firstRecord;
        record2 = secondRecord;
    }

    public void hold(VCFRecord record) {
        holder.add(record);
    }

    public boolean hasHeldRecords() {
        return holder.size() > 0;
    }

    public VCFRecord getHeldRecord() {
        if ( samples == null ) {
            throw new IllegalStateException("Held records are being accessed before intersection sample set have been created");
        }

        return downsample(holder.remove(0));
    }

    private void instantiateGenotypeMerger() {
        if ( record1 != null && record2 != null ) {
            ArrayList<String> keyset1 =  new ArrayList<String>(Arrays.asList(record1.getGenotypeFormatString().split(":")));
            ArrayList<String> keyset2 =  new ArrayList<String>(Arrays.asList(record2.getGenotypeFormatString().split(":")));
            keyset2.retainAll(keyset1);
            genotypeMerger = new VCFGenotypeMerger(keyset2);
        }
    }

    public VCFRecord getMergedRecord() {
        if ( record2 == null ) {
            return downsample(record1);
        } else if ( record1 == null ) {
            return downsample(record2);
        } else {
            HashMap<VCFHeader.HEADER_FIELDS,String> newFields = getFields(record1,record2);
            List<VCFGenotypeRecord> newGenotypes = mergeGenotypes(record1,record2);
            return new VCFRecord(newFields,genotypeMerger.getFormatString(),newGenotypes);
        }
    }

    private VCFRecord downsample(VCFRecord record) {
        HashMap<VCFHeader.HEADER_FIELDS,String> newFields = getFields(record,null);
        List<VCFGenotypeRecord> newGenotypes = new ArrayList<VCFGenotypeRecord>();
        for ( String s : samples ) {
            newGenotypes.add(record.getGenotype(s) );
        }

        return new VCFRecord(newFields,record.getGenotypeFormatString(),newGenotypes);
    }

    private List<VCFGenotypeRecord> mergeGenotypes(VCFRecord first, VCFRecord second) {
        List<VCFGenotypeRecord> genotypes = new ArrayList<VCFGenotypeRecord>();
        for ( String s : samples ) {
            if ( isMergeByQual ) {
                genotypes.add(genotypeMerger.mergeByQual(first.getGenotype(s), second.getGenotype(s)));
            } else if ( record1Priority ) {
                genotypes.add(genotypeMerger.resetKeys(first.getGenotype(s)));
            } else {
                genotypes.add(genotypeMerger.resetKeys(second.getGenotype(s)));
            }
        }

        return genotypes;
    }

    private HashMap<VCFHeader.HEADER_FIELDS,String> getFields(VCFRecord record1, VCFRecord record2) {
        HashMap<VCFHeader.HEADER_FIELDS,String> fields = new HashMap<VCFHeader.HEADER_FIELDS,String>();
        fields.put(VCFHeader.HEADER_FIELDS.CHROM, record1.getLocation().getContig());
        fields.put(VCFHeader.HEADER_FIELDS.POS, String.format("%d",record1.getLocation().getStart()));
        fields.put(VCFHeader.HEADER_FIELDS.ID, record1.getID());
        fields.put(VCFHeader.HEADER_FIELDS.ALT, listAsString(record1.getAlternateAlleleList()));
        fields.put(VCFHeader.HEADER_FIELDS.REF, record1.getReference());
        fields.put(VCFHeader.HEADER_FIELDS.QUAL, String.format("%f", record1.getQual()));
        if ( record2 != null && record1.isFiltered() && record2.isFiltered() ) {
            fields.put(VCFHeader.HEADER_FIELDS.FILTER,"BOTH_FILTERED");
        } else if ( record1.isFiltered() || ( record2 != null && record2.isFiltered() ) ) {
            fields.put(VCFHeader.HEADER_FIELDS.FILTER,"ONE_FILTERED");
        }
        return fields;
    }

    private String listAsString(List<String> lString) {
        StringBuilder b = new StringBuilder();
        b.append(lString.remove(0));
        for ( String s : lString ) {
            b.append(",");
            b.append(s);
        }

        return b.toString();
    }

    public boolean hasBeenInstantiated() {
        return genotypeMerger != null;
    }

}

class VCFGenotypeMerger {
    private List<String> intersectKeys;
    private String intersectFormatString;

    public VCFGenotypeMerger(List<String> keys) {
        intersectKeys = keys;
        intersectFormatString = keyJoin(keys);
    }

    public String getFormatString() {
        return intersectFormatString;
    }

    public VCFGenotypeRecord mergeByQual(VCFGenotypeRecord record1, VCFGenotypeRecord record2) {
        VCFGenotypeRecord toReturn;
        if ( record2 == null || record2.isNoCall() || record2.isFiltered() || record1.getNegLog10PError() > record2.getNegLog10PError() ) {
            toReturn = record1;
        } else {
            toReturn = record2;
        }

        resetKeys(toReturn);

        return toReturn;
    }

    public VCFGenotypeRecord resetKeys(VCFGenotypeRecord r) {

        HashMap<String,String> newFields = new HashMap<String,String>();

        for ( String key : intersectKeys ) {
            newFields.put(key,r.getFields().get(key));
        }

        r.replaceFields(newFields);

        return r;
    }


    private String keyJoin(List<String> sList) {
        StringBuffer format = new StringBuffer(sList.get(0));
        for ( int i = 1; i < sList.size(); i ++ ) {
            format.append(":");
            format.append(sList.get(i));
        }

        return format.toString();
    }


}