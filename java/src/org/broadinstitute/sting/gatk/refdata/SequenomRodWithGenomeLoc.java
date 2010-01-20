package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 19, 2010
 * Time: 10:24:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class SequenomRodWithGenomeLoc extends BasicReferenceOrderedDatum implements ReferenceOrderedDatum {
    private final String[] SEQUENOM_HEADER_FIELDS = { "#Family ID", "Individual ID", "Paternal ID", "Maternal ID", "Sex", "Phenotype" } ;

    private ArrayList<SequenomVariantInfo> variants;
    private SequenomVariantInfo currentVariant;
    private ListIterator<SequenomVariantInfo> variantIterator;
    private HashSet headerEntries;

    // // // CONSTRUCTOR // // //

    public SequenomRodWithGenomeLoc(String name) {
        super(name);
    }

    @Override
    public Object initialize(final File seqFile) throws FileNotFoundException {
        if ( ! seqFile.exists() ) {
            throw new FileNotFoundException("File "+seqFile.getAbsolutePath()+" does not exist.");
        }

        headerEntries = new HashSet<String>(Arrays.asList(SEQUENOM_HEADER_FIELDS));

        variants = parseSequenomFile(seqFile);
        if ( variants != null ) {
            variantIterator = variants.listIterator();
            currentVariant = variantIterator.next();
        }

        assertNotNull();

        return null;
    }

    private void assertNotNull() {
        if ( currentVariant == null ) {
            throw new UnsupportedOperationException ( "Current sequenom variant information was set to null" );
        }
    }

    @Override
    public boolean parseLine(Object obj, String[] args) {
        if ( variantIterator.hasNext() ) {
            currentVariant = variantIterator.next();
            return true;
        } else {
            return false;
        }
    }

    @Override
    public GenomeLoc getLocation() {
        return currentVariant.getLocation();
    }

    @Override
    public String toString() {
        return currentVariant == null ? "" : currentVariant.toString();
    }

    public String getVariantName() {
        return currentVariant.getName();
    }

    public ArrayList<String> getVariantSampleNames() {
        return currentVariant.getSampleNames();
    }

    public ArrayList<Genotype> getGenotypes() {
        return currentVariant.getGenotypes();
    }

    private ArrayList<SequenomVariantInfo> parseSequenomFile(File file) {
        try {
            BufferedReader reader = new BufferedReader( new FileReader ( file ) );
            String header = reader.readLine();
            ArrayList<SequenomVariantInfo> seqVars = instantiateVariantListFromHeader(header);
            ArrayList<Integer> snpOffsets = getSNPOffsetsFromHeader(header);

            String line = null;
            do {
                line = reader.readLine();
                incorporateInfo(seqVars,snpOffsets,line);
            } while ( line != null );

            java.util.Collections.sort(seqVars); // because the comparable uses the GenomeLoc comparable; these
                                                 // are sorted in standard reference order

            return seqVars;

        } catch ( FileNotFoundException e ) {
            throw new StingException("File "+file.getAbsolutePath()+" could not be found. This was checked earlier. Should never happen.",e);
        } catch ( IOException e ) {
            throw new StingException("Error reading file "+file.getAbsolutePath()+".",e);
        }
    }

    private void incorporateInfo(List<SequenomVariantInfo> vars, List<Integer> offsets, String seqLine) {
        String[] genotypes = seqLine.split("\t");
        String individualName = genotypes[1];

        int snpNumber = 0;
        for ( int i : offsets ) {
            vars.get(snpNumber).addGenotypeEntry(genotypes[i], individualName);
            snpNumber++;
        }
    }

    private ArrayList<SequenomVariantInfo> instantiateVariantListFromHeader(String header) {
        ArrayList<SequenomVariantInfo> seqVars = new ArrayList<SequenomVariantInfo>();
        String[] headerFields = header.split("\t");

        for ( String field : headerFields ) {
            if ( ! headerEntries.contains(field) ) {
                // not a standard header, so a variant
                seqVars.add(new SequenomVariantInfo(field));
            }
        }

        return seqVars;
    }

    private ArrayList<Integer> getSNPOffsetsFromHeader(String header) {
        ArrayList<Integer> offsets = new ArrayList<Integer>();
        String[] headerFields = header.split("\t");

        int offset = 0;
        for ( String field : headerFields ) {
            if ( ! headerEntries.contains(field) ) {
                offsets.add(offset);
            }
            offset++;
        }

        return offsets;
    }
}

class SequenomVariantInfo implements Comparable {
    private String variantName;
    private GenomeLoc loc;
    private ArrayList<Genotype> genotypes;
    private ArrayList<String> sampleNames;

    private ArrayList<Allele> deletionHolder = new ArrayList<Allele>();
    private ArrayList<String> sampleHolder = new ArrayList<String>();
    private int siteDeletionLength = -1;

    public GenomeLoc getLocation() {
        return loc;
    }

    public String getName() {
        return variantName;
    }

    public ArrayList<String> getSampleNames() {
        return sampleNames;
    }

    public ArrayList<Genotype> getGenotypes() {
        return genotypes;
    }

    // CONSTRUCTOR

    public SequenomVariantInfo(String variantName) {
        this.variantName = variantName;
        this.parseNameToLoc();
    }

    private void parseNameToLoc() {
        String chrom = this.variantName.split("_c")[1].split("_")[0];
        String pos = this.variantName.split("_p")[1].split("_")[0];
        this.loc = GenomeLocParser.parseGenomeLoc(chrom+":"+pos);
    }

    public void addGenotypeEntry(String genotypeString, String sampleName) {
        String[] alleleStrs = genotypeString.split(" ");
        // identify if we're dealing with a deletion
        if ( genotypeString.contains("-") ) {
            this.addDeletion(alleleStrs, sampleName);
        } else {
            // simple SNP or indel (easier to handle)
            this.addIndelOrSNP(alleleStrs,sampleName);
        }
    }

    private void addIndelOrSNP(String[] alleleStrings, String sampleName) {
        ArrayList<Allele> alleles = new ArrayList<Allele>(2);

        if ( alleleStrings[0].length() > 1 || alleleStrings[1].length() > 1 ) {
            // insertion
            for ( String alStr : alleleStrings ) {
                if ( alStr.length() > 1 ) {
                    alleles.add(new Allele(Allele.AlleleType.INSERTION,alStr));
                } else {
                    alleles.add(new Allele(Allele.AlleleType.REFERENCE, alStr));
                }
            }
        } else {
            // SNP
            for ( String alStr : alleleStrings ) {
                alleles.add(new Allele(Allele.AlleleType.UNKNOWN_POINT_ALLELE,alStr));
            }
        }
    }

    private void addDeletion(String[] alleleStrings, String sampleName) {
        String alleleStr1 = alleleStrings[0];
        String alleleStr2 = alleleStrings[1];
        Allele allele1 = null;
        Allele allele2 = null;

        if ( alleleStr1.contains("-") && alleleStr2.contains("-") ) {
            // homozygous deletion
            this.addHomDeletion(allele1,allele2, sampleName);
        } else {
            // heterozygous deletion
            if ( alleleStr1.contains("-") ) {
                this.addHetDeletion(allele1,allele2, alleleStr1, alleleStr2, sampleName);
            } else {
                this.addHetDeletion(allele2,allele1, alleleStr2, alleleStr1, sampleName); // note the order change
            }
        }
    }

    private void addHetDeletion(Allele del, Allele ref, String delStr, String refStr, String sampleName) {
        del = new Allele(Allele.AlleleType.DELETION,"");
        ref = new Allele(Allele.AlleleType.REFERENCE,refStr.substring(0,1));
        this.setDeletionLength(del,refStr.length());
        if ( ! deletionHolder.isEmpty() ) {
            siteDeletionLength = refStr.length();
            this.addHeldDeletions();
        }
        Genotype indel = new Genotype(Arrays.asList(ref,del), sampleName, 20.0);
        this.setIndelGenotypeLength(indel,siteDeletionLength);
        this.genotypes.add(indel);
        this.sampleNames.add(sampleName);
    }

    private void addHomDeletion(Allele allele1, Allele allele2, String sampleName) {
        allele1 = new Allele(Allele.AlleleType.DELETION,"");
        allele2 = new Allele(Allele.AlleleType.DELETION,"");
        if ( siteDeletionLength != -1 ) {
            this.setDeletionLength(allele1,siteDeletionLength);
            this.setDeletionLength(allele2,siteDeletionLength);
            Genotype indel = new Genotype(Arrays.asList(allele1,allele2), sampleName, 20.0);
            this.setIndelGenotypeLength(indel, siteDeletionLength);
            this.genotypes.add(indel);
            this.sampleNames.add(sampleName);
        } else {
            deletionHolder.add(allele1);
            deletionHolder.add(allele2);
            sampleHolder.add(sampleName);
        }
    }

    private void setIndelGenotypeLength(Genotype g, int length) {
        g.setAttribute(Genotype.StandardAttributes.DELETION_LENGTH,length);
    }

    private void addHeldDeletions() {
        Allele del1;
        Allele del2;
        int startingSize = deletionHolder.size();
        for ( int i = 0; i < startingSize ; i+=2 ) {
            del1 = deletionHolder.get(i);
            del2 = deletionHolder.get(i+1);
            this.addHomDeletion(del1,del2,sampleHolder.get(i/2));
            if ( deletionHolder.size() != startingSize ) {
                throw new StingException("Halting algorithm -- possible infinite loop");
            }
        }
        deletionHolder.clear();
        sampleHolder.clear();
    }

    public int compareTo(Object obj) {
        if ( ! ( obj instanceof SequenomVariantInfo ) ) {
            return 1;
        }

        return loc.compareTo(((SequenomVariantInfo) obj).getLocation());
    }

    private void setDeletionLength(Allele al, int length) {
        // Todo -- once alleles support deletion lengths add that information
        // Todo -- into the object; for now this can just return
        return;
    }
}