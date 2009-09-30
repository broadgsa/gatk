package org.broadinstitute.sting.gatk.refdata;

import net.sf.samtools.util.SequenceUtil;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.BasicGenotype;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Example format:
 * 585 chr1 433 433 rs56289060  0  +  - - -/C  genomic  insertion unknown 0  0  unknown  between  1
 * 585 chr1 491 492 rs55998931  0  +  C C C/T  genomic  single   unknown 0 0 unknown exact 1
 * <p/>
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class rodDbSNP extends BasicReferenceOrderedDatum implements VariationRod, VariantBackedByGenotype, AllelicVariant {
    public GenomeLoc loc;       // genome location of SNP
                                // Reference sequence chromosome or scaffold
                                // Start and stop positions in chrom

    public String name;        // Reference SNP identifier or Affy SNP name
    public String strand;      // Which DNA strand contains the observed alleles

    public String refBases;        // the reference base according to NCBI, in the dbSNP file
    public String observed;    // The sequences of the observed alleles from rs-fasta files

    public String molType;     // Sample type from exemplar ss
    public String varType;     // The class of variant (simple, insertion, deletion, range, etc.)
                                // Can be 'unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion'
    public String validationStatus;    // The validation status of the SNP
                                        // one of set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap')

    public double avHet;        // The average heterozygosity from all observations
    public double avHetSE;      // The Standard Error for the average heterozygosity

    public String func;         // The functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)
                                // set('unknown','coding-synon','intron','cds-reference','near-gene-3','near-gene-5',
                                // 'nonsense','missense','frameshift','untranslated-3','untranslated-5','splice-3','splice-5')
    public String locType;      // How the variant affects the reference sequence
                                // enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion')

    public int weight;          // The quality of the alignment

    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public rodDbSNP(final String name) {
        super(name);
    }

    // ----------------------------------------------------------------------
    //
    // manipulating the SNP information
    //
    // ----------------------------------------------------------------------
    public GenomeLoc getLocation() {
        return loc;
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    @Override
    public String getReference() {
        return getRefBasesFWD();
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        return 4; // -log10(0.0001)
    }

    public boolean onFwdStrand() {
        return strand.equals("+");
    }

    /**
     * Returns bases in the reference allele as a String. String can be empty (as in insertion into
     * the reference), can contain a single character (as in SNP or one-base deletion), or multiple characters
     * (for longer indels).
     *
     * @return reference allele, forward strand
     */
    public String getRefBasesFWD() {
        // fix - at least this way we ensure that we'll get the other base compared to getAltBasesFWD()
        return (getAllelesFWD().get(0).equals(refBases)) ? getAllelesFWD().get(0) : getAllelesFWD().get(1);
        //if ( onFwdStrand() )
        //    return refBases;
        //else
        //    return SequenceUtil.reverseComplement(refBases);
    }

    /**
     * Returns reference (major) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * @return reference base on the forward strand
     */
    public char getRefSnpFWD() throws IllegalStateException {
        //System.out.printf("refbases is %s but %s%n", refBases, toString());
        if (isIndel()) throw new IllegalStateException("Variant is not a SNP");
        // fix - at least this way we ensure that we'll get the other base compared to getAltBasesFWD()
        List<String> alleles = getAllelesFWD();
        String val = (alleles.get(0).equals(refBases) ? alleles.get(0) : alleles.get(1));
        return val.charAt(0);
        // if ( onFwdStrand() ) return refBases.charAt(0);
        // else return SequenceUtil.reverseComplement(refBases).charAt(0);
    }

    public List<String> getAllelesFWD() {
        List<String> alleles = null;
        if (onFwdStrand())
            alleles = Arrays.asList(observed.split("/"));
        else
            alleles = Arrays.asList(SequenceUtil.reverseComplement(observed).split("/"));

        //System.out.printf("getAlleles %s on %s %b => %s %n", observed, strand, onFwdStrand(), Utils.join("/", alleles));
        return alleles;
    }

    public String getAllelesFWDString() {
        return Utils.join("", getAllelesFWD());
    }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    @Override
    public double getNonRefAlleleFrequency() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /** @return the VARIANT_TYPE of the current variant */
    @Override
    public VARIANT_TYPE getType() {
        return VARIANT_TYPE.SNP;
    }// ----------------------------------------------------------------------

    //
    // What kind of variant are we?
    //
    // ----------------------------------------------------------------------
    public boolean isSNP() {
        return varType.contains("single") && locType.contains("exact");
    }

    public boolean isInsertion() {
        return varType.contains("insertion");
    }

    public boolean isDeletion() {
        return varType.contains("deletion");
    }

    /**
     * get the base representation of this Variant
     *
     * @return a string, of ploidy
     */
    @Override
    public String getAlternateBases() {
        return getAllelesFWDString();
    }

    /**
     * gets the alternate bases.  Use this method if teh allele count is greater then 2
     *
     * @return
     */
    @Override
    public List<String> getAlternateBaseList() {
        List<String> list = new ArrayList<String>();
        list.add(this.getAlternateBases());
        return list;
    }

    public boolean isIndel() {
        return isInsertion() || isDeletion() || varType.contains("in-del");
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException in the case
     * of
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getAlternativeBaseForSNP() {
        return getAltSnpFWD(); /*
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        if (getAlternateBases().charAt(0) == this.getReference())
            return getAlternateBases().charAt(1);
        return getAlternateBases().charAt(0); */
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isReference() {
        return false;
    } // snp locations are never "reference", there's always a variant

    public boolean isHapmap() {
        return validationStatus.contains("by-hapmap");
    }

    public boolean is2Hit2Allele() {
        return validationStatus.contains("by-2hit-2allele");
    }

    // ----------------------------------------------------------------------
    //
    // formatting
    //
    // ----------------------------------------------------------------------
    public String toString() {
        return String.format("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%s\t%s\t%d",
                             getLocation().getContig(), getLocation().getStart(), getLocation().getStop() + 1,
                             name, strand, refBases, observed, molType,
                             varType, validationStatus, avHet, avHetSE, func, locType, weight);
    }

    public String toSimpleString() {
        return String.format("%s:%s:%s", name, observed, strand);
    }

    public String toMediumString() {
        String s = String.format("%s:%s:%s", getLocation().toString(), name, getAllelesFWDString());
        if (isSNP()) s += ":SNP";
        if (isIndel()) s += ":Indel";
        if (isHapmap()) s += ":Hapmap";
        if (is2Hit2Allele()) s += ":2Hit";
        return s;
    }

    public String repl() {
        return String.format("%d\t%s\t%d\t%d\t%s\t0\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%s\t%s\t%d",
                             585, getLocation().getContig(), getLocation().getStart() - 1, getLocation().getStop(),
                             name, strand, refBases, refBases, observed, molType,
                             varType, validationStatus, avHet, avHetSE, func, locType, weight);
    }

    public boolean parseLine(final Object header, final String[] parts) {
        try {
            String contig = parts[1];
            long start = Long.parseLong(parts[2]) + 1; // The final is 0 based
            long stop = Long.parseLong(parts[3]) + 1;  // The final is 0 based
            loc = GenomeLocParser.parseGenomeLoc(contig, start, Math.max(start, stop - 1));

            name = parts[4];
            strand = parts[6];
            refBases = parts[7];
            if (strand == "-")
                refBases = BaseUtils.simpleReverseComplement(refBases);
            observed = parts[9];
            molType = parts[10];
            varType = parts[11];
            validationStatus = parts[12];
            avHet = Double.parseDouble(parts[13]);
            avHetSE = Double.parseDouble(parts[14]);
            func = parts[15];
            locType = parts[16];
            weight = Integer.parseInt(parts[17]);
            //System.out.printf("Parsed %s%n", toString());
            return true;
        } catch (MalformedGenomeLocException ex) {
            // Just rethrow malformed genome locs; the ROD system itself will deal with these.
            throw ex;
        } catch (ArrayIndexOutOfBoundsException ex) {
            // Just rethrow malformed genome locs; the ROD system itself will deal with these.
            throw new RuntimeException("Badly formed dbSNP line: " + ex);
        } catch (RuntimeException e) {
            System.out.printf("  Exception caught during parsing DBSNP line %s%n", Utils.join(" <=> ", parts));
            throw e;
        }
    }

    public String getAltBasesFWD() {
        List<String> alleles = getAllelesFWD();
        return (alleles.get(0).equals(refBases) ? alleles.get(1) : alleles.get(0));
    }

    public char getAltSnpFWD() throws IllegalStateException {
        if (!isSNP())
            throw new IllegalStateException("I'm not a SNP");
        return getAltBasesFWD().charAt(0);
    }

    public double getConsensusConfidence() {
        // TODO Auto-generated method stub
        return Double.MAX_VALUE;
    }

    public List<String> getGenotype() throws IllegalStateException {
        return Arrays.asList(Utils.join("", getAllelesFWD()));
    }

    public double getMAF() {
        // Fixme: update to actually get MAF
        //return avHet;
        return -1;
    }

    public double getHeterozygosity() {
        return avHet;
    }

    public int getPloidy() throws IllegalStateException {
        // TODO Auto-generated method stub
        return 0;
    }

    public double getVariationConfidence() {
        // TODO Auto-generated method stub
        return Double.MAX_VALUE;
    }

    public boolean isGenotype() {
        // TODO Auto-generated method stub
        return false;
    }

    public boolean isBiallelic() {
        // TODO Auto-generated method stub
        return observed.indexOf('/') == observed.lastIndexOf('/');
    }

    public int length() {
        return (int) (loc.getStop() - loc.getStart() + 1);
    }

    /**
     * get the genotype
     *
     * @return a map in lexigraphical order of the genotypes
     */
    @Override
    public org.broadinstitute.sting.utils.genotype.Genotype getCalledGenotype() {
        return new BasicGenotype(this.getLocation(), this.getAltBasesFWD(), this.getRefSnpFWD(), this.getConsensusConfidence());
    }

    /**
     * get the likelihoods
     *
     * @return an array in lexigraphical order of the likelihoods
     */
    @Override
    public List<org.broadinstitute.sting.utils.genotype.Genotype> getGenotypes() {
        List<org.broadinstitute.sting.utils.genotype.Genotype> list = new ArrayList<org.broadinstitute.sting.utils.genotype.Genotype>();
        list.add(new BasicGenotype(this.getLocation(), this.getAltBasesFWD(), this.getRefSnpFWD(), this.getConsensusConfidence()));
        return list;
    }

    /**
     * do we have the specified genotype?  not all backedByGenotypes
     * have all the genotype data.
     *
     * @param x the genotype
     *
     * @return true if available, false otherwise
     */
    @Override
    public boolean hasGenotype(DiploidGenotype x) {
        return (!x.toString().equals(this.getAltBasesFWD())) ? false : true;
    }

    public static rodDbSNP getFirstRealSNP(RODRecordList<ReferenceOrderedDatum> dbsnpList) {
	if ( dbsnpList == null )
	    return null;

	rodDbSNP dbsnp = null;
	for ( ReferenceOrderedDatum d : dbsnpList ) {
	    if ( ((rodDbSNP)d).isSNP() ) {
		dbsnp = (rodDbSNP)d;
		break;
	    }
	}

	return dbsnp;
    }
}
