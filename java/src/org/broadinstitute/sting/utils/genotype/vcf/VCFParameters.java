package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;
import java.util.ArrayList;


/**
 * a helper class, which performs a lot of the safety checks on the parameters
 * we feed to the VCF (like ensuring the same position for each genotype in a call).
 */
class VCFParameters {
    private char referenceBase = '0';
    private int position = 0;
    private String contig = null;
    private boolean initialized = false;
    private List<VCFGenotypeRecord> genotypesRecord = new ArrayList<VCFGenotypeRecord>();
    private List<String> formatList = new ArrayList<String>();
    private List<VCFGenotypeEncoding> alternateBases = new ArrayList<VCFGenotypeEncoding>();

    public void setLocations(GenomeLoc location, char refBase) {
        // if we haven't set it up, we initialize the object
        if (!initialized) {
            initialized = true;
            this.contig = location.getContig();
            this.position = (int) location.getStart();
            if (location.getStart() != location.getStop()) {
                throw new IllegalArgumentException("The start and stop locations must be the same");
            }
            this.referenceBase = refBase;
        } else {
            if (!contig.equals(this.contig))
                throw new IllegalArgumentException("The contig name has to be the same at a single locus");
            if (location.getStart() != this.position)
                throw new IllegalArgumentException("The position has to be the same at a single locus");
            if (refBase != this.referenceBase)
                throw new IllegalArgumentException("The reference base name has to be the same at a single locus");
        }
    }

    /** @return get the position */
    public int getPosition() {
        return position;
    }

    /** @return get the contig name */
    public String getContig() {
        return contig;
    }

    /** @return get the reference base */
    public char getReferenceBase() {
        return referenceBase;
    }

    public void addGenotypeRecord(VCFGenotypeRecord record) {
        this.genotypesRecord.add(record);
    }

    public void addFormatItem(String item) {
        if (!formatList.contains(item))
            formatList.add(item);
    }

    public void addAlternateBase(VCFGenotypeEncoding base) {
        if ( !alternateBases.contains(base) &&
             !base.toString().equals(String.valueOf(getReferenceBase())) &&
             !base.toString().equals(VCFGenotypeRecord.EMPTY_ALLELE) )
            alternateBases.add(base);
    }

    public List<VCFGenotypeEncoding> getAlternateBases() {
        return alternateBases;
    }

    public String getFormatString() {
        return Utils.join(VCFRecord.FORMAT_FIELD_SEPERATOR, formatList);
    }

    public List<VCFGenotypeRecord> getGenotypesRecords() {
        return genotypesRecord;
    }
}
