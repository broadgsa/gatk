package org.broadinstitute.sting.gatk.refdata;
/*
import org.broadinstitute.sting.oneoffprojects.variantcontext.Allele;
import org.broadinstitute.sting.oneoffprojects.variantcontext.Genotype;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.CommandLineGATK;

import java.io.*;
import java.util.*;

import net.sf.samtools.SAMFileHeader;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 19, 2010
 * Time: 10:24:18 AM
 * To change this template use File | Settings | File Templates.
 *
public class PlinkRod extends BasicReferenceOrderedDatum implements ReferenceOrderedDatum {
    private final Set<String> headerEntries = new HashSet<String>(Arrays.asList("#Family ID","Individual ID","Sex",
                "Paternal ID","Maternal ID","Phenotype", "FID","IID","PAT","MAT","SEX","PHENOTYPE"));
    private final byte SNP_MAJOR_MODE = 1;

    private ArrayList<PlinkVariantInfo> variants;
    private PlinkVariantInfo currentVariant;
    private ListIterator<PlinkVariantInfo> variantIterator;

    private PlinkFileType plinkFileType;

    public enum PlinkFileType {
        STANDARD_PED,RAW_PED,BINARY_PED
    }

//       CONSTRUCTOR

    public PlinkRod(String name) {
        super(name);
    }

    @Override
    public Object initialize(final File plinkFile) throws FileNotFoundException {
        if ( ! plinkFile.exists() ) {
            throw new FileNotFoundException("File "+plinkFile.getAbsolutePath()+" does not exist.");
        }

        variants = parsePlinkFile(plinkFile);
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

    public boolean variantIsSNP() {
        return currentVariant.isSNP();
    }


//     AM I PARSING A TEXT OR A BINARY FILE ??

    private ArrayList<PlinkVariantInfo> parsePlinkFile(File file) {
        String[] splitFileName = file.getName().split("\\.");
        String extension = splitFileName[splitFileName.length-1];
        if ( extension.equals("ped") || extension.equals("raw") ) {
            return parseTextFormattedPlinkFile(file);
        } else if ( extension.equals("bed") || extension.equals("bim") || extension.equals("fam") ) {
            plinkFileType = PlinkFileType.BINARY_PED;
            return parseBinaryFormattedPlinkFile(file);
        } else {
            System.out.println("Warning -- Plink file does not have a standard extension (ped/raw for text, bed/bim/fam for binary) -- assuming ped format");
            return parseTextFormattedPlinkFile(file);
        }

    }

    /* *** *** *** *** *** ** *** ** *** ** *** ** *** ** ***
     * *    PARSING    STANDARD   TEXT   PED    FILES    * **
     * *** *** *** *** *** ** *** ** *** ** *** ** *** ** ***/
/*
    private ArrayList<PlinkVariantInfo> parseTextFormattedPlinkFile( File file ) {
        try {
            BufferedReader reader = new BufferedReader( new FileReader ( file ) );
            String header = reader.readLine();
            ArrayList<PlinkVariantInfo> seqVars = instantiateVariantListFromHeader(header);
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

    private void incorporateInfo(List<PlinkVariantInfo> vars, List<Integer> offsets, String plinkLine) {
        if ( plinkLine == null ) {
            return;
        }

        String[] plinkInfo;
        if ( plinkFileType == PlinkFileType.STANDARD_PED) {
            plinkInfo = plinkLine.split("\t");
        } else {
            throw new StingException("Plink file is likely of .raw or recoded format. Please use an uncoded .ped file.");
        }

        String individualName = plinkInfo[1];

        int snpNumber = 0;

        if ( plinkFileType == PlinkFileType.STANDARD_PED ) {
            for ( int i : offsets ) {
                vars.get(snpNumber).addGenotypeEntry(plinkInfo[i], individualName);
                snpNumber++;
            }
        }
    }

    private ArrayList<PlinkVariantInfo> instantiateVariantListFromHeader(String header) {
        if ( header.startsWith("#") ) {// first line is a comment; what we're used to seeing
            plinkFileType = PlinkFileType.STANDARD_PED;
        } else {// first line is the raw header; comes from de-binary-ing a .bed file
            plinkFileType = PlinkFileType.RAW_PED;
            throw new StingException("Plink file is likely of .raw or recoded format. Please use an uncoded .ped file.");
        }

        ArrayList<PlinkVariantInfo> seqVars = new ArrayList<PlinkVariantInfo>();
        String[] headerFields;

        if ( plinkFileType == PlinkFileType.STANDARD_PED ) {
            headerFields = header.split("\t");
        } else {
            throw new StingException("Plink file is likely of .raw or recoded format. Please use an uncoded .ped file.");
        }

        for ( String field : headerFields ) {
            if ( ! headerEntries.contains(field) ) {
                 // not a standard header, so a variant
                seqVars.add(new PlinkVariantInfo(field));
            }
        }

        return seqVars;
    }

    private ArrayList<Integer> getSNPOffsetsFromHeader(String header) {
        ArrayList<Integer> offsets = new ArrayList<Integer>();
        String[] headerFields;
        if ( plinkFileType == PlinkFileType.STANDARD_PED ) {
            headerFields = header.split("\t+");
        } else {
            headerFields = header.split("\\s+");
        }

        int offset = 0;

        for ( String field : headerFields ) {
            if ( ! headerEntries.contains(field) ) {
                offsets.add(offset);
            }
            offset++;
        }

        return offsets;
    }

    /* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
     * *    PARSING     BINARY    PLINK   BED/BIM/FAM   FILES      *  *
     * *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***/
 /*
    private ArrayList<PlinkVariantInfo> parseBinaryFormattedPlinkFile(File file) {
        PlinkBinaryTrifecta binaryFiles = getPlinkBinaryTriplet(file);
        ArrayList<PlinkVariantInfo> parsedVariants = instantiateVariantsFromBimFile(binaryFiles.bimFile);
        ArrayList<String> sampleNames = getSampleNameOrderingFromFamFile(binaryFiles.famFile);
        ArrayList<PlinkVariantInfo> updatedVariants = getGenotypesFromBedFile(parsedVariants,sampleNames,binaryFiles.bedFile);

        java.util.Collections.sort(updatedVariants);

        return updatedVariants;
    }

    private PlinkBinaryTrifecta getPlinkBinaryTriplet(File file) {
         // just gonna parse the name
        PlinkBinaryTrifecta trifecta = new PlinkBinaryTrifecta();
        String absolute_path = file.getAbsolutePath();
        String[] directory_tree = absolute_path.split("/");
        String file_name = directory_tree[directory_tree.length-1].split("\\.")[0];
        StringBuilder pathBuilder = new StringBuilder();
        for ( int i = 0; i < directory_tree.length - 1; i ++ ) {
            pathBuilder.append(String.format("%s/",directory_tree[i]));
        }
        String path = pathBuilder.toString();
        trifecta.bedFile = new File(path+file_name+".bed");
        trifecta.bimFile = new File(path+file_name+".bim");
        trifecta.famFile = new File(path+file_name+".fam");

        return trifecta;

    }

    private ArrayList<PlinkVariantInfo> instantiateVariantsFromBimFile(File bimFile) {
        BufferedReader reader;
        try {
            reader = new BufferedReader( new FileReader( bimFile ));
        } catch ( FileNotFoundException e) {
            throw new StingException("The SNP information file accompanying the binary ped file was not found (the .bim file). "+
                                     "Please ensure that it is in the same directory as the .bed and .fam files. The file we "+
                                     "Were looking for was "+bimFile.getAbsolutePath(),e);
        }

        ArrayList<PlinkVariantInfo> variants = new ArrayList<PlinkVariantInfo>();

        try {
            String line = null;
            do {
                line = reader.readLine();
                if ( line != null ) {
                    String[] snpInfo = line.split("\\s+");
                    PlinkVariantInfo variant = new PlinkVariantInfo(snpInfo[1],true);
                    variant.setGenomeLoc(GenomeLocParser.parseGenomeLoc(snpInfo[0],Long.valueOf(snpInfo[3]), Long.valueOf(snpInfo[3])));
                    variant.setAlleles(snpInfo[4],snpInfo[5]);
                    variants.add(variant);
                }
            } while ( line != null );
        } catch ( IOException e ) {
            throw new StingException("There was an error reading the .bim file "+bimFile.getAbsolutePath(), e);
        }

        return variants;
    }

    private ArrayList<String> getSampleNameOrderingFromFamFile(File famFile) {
        BufferedReader reader;
        try {
            reader = new BufferedReader( new FileReader( famFile ));
        } catch ( FileNotFoundException e) {
            throw new StingException("The Family information file accompanying the binary ped file was not found (the .fam file). "+
                                     "Please ensure that it is in the same directory as the .bed and .bim files. The file we "+
                                     "Were looking for was "+famFile.getAbsolutePath(),e);
        }

        ArrayList<String> sampleNames = new ArrayList<String>();

        try {
            String line = null;
            do {
                line = reader.readLine();
                if ( line != null ) {
                    sampleNames.add(line.split("\\s+")[1]);
                }
            } while ( line != null );
        } catch (IOException e) {
            throw new StingException("There was an error reading the .fam file "+famFile.getAbsolutePath(),e);
        }

        return sampleNames;
    }

    private ArrayList<PlinkVariantInfo> getGenotypesFromBedFile(ArrayList<PlinkVariantInfo> variants, ArrayList<String> samples, File bedFile) {
        FileInputStream inStream;
        try {
            inStream = new FileInputStream(bedFile);
        } catch (FileNotFoundException e) {
            throw new StingException("The Binary pedigree file file accompanying the family file was not found (the .bed file). "+
                                     "Please ensure that it is in the same directory as the .bim and .fam files. The file we "+
                                     "Were looking for was "+bedFile.getAbsolutePath(),e);
        }

        try {
            byte genotype = -1;
            long bytesRead = 0;
            int snpOffset = 0;
            int sampleOffset = 0;
            boolean snpMajorMode = true;
            do {
                genotype = (byte) inStream.read();
                bytesRead++;
                if ( genotype != -1 ) {
                    if ( bytesRead > 3 ) {
                        addGenotypeByte(genotype,variants,samples,snpOffset,sampleOffset, snpMajorMode);

                        if ( snpMajorMode ) {
                            sampleOffset = sampleOffset + 4;
                            if ( sampleOffset > samples.size() -1 ) {
                                snpOffset ++;
                                sampleOffset = 0;
                            }
                        } else {
                            snpOffset = snpOffset + 4;
                            if ( snpOffset > variants.size() -1 ) {
                                sampleOffset ++;
                                snpOffset = 0;
                            }
                        }

                    } else {
                        if ( bytesRead == 3) {
                            snpMajorMode = genotype == SNP_MAJOR_MODE;
                        }
                    }
                }
            } while ( genotype != -1 );
        } catch ( IOException e) {
            throw new StingException("Error reading binary ped file "+bedFile.getAbsolutePath(), e);
        }

        return variants;
    }

    private void addGenotypeByte(byte genotype, ArrayList<PlinkVariantInfo> variants, ArrayList<String> sampleNames, int snpOffset, int sampleOffset, boolean major) {
        // four genotypes encoded in this byte
        int[] genotypes = parseGenotypes(genotype);
        for ( int g : genotypes ) {
            variants.get(snpOffset).addBinaryGenotypeEntry(g,sampleNames.get(sampleOffset));

            if ( major ) {
                sampleOffset++;
                if ( sampleOffset > sampleNames.size()-1 ) {//  using offsets for comparison; size 5 == offset 4
                    return;
                }
            } else {
                snpOffset++;
                if( snpOffset > variants.size()-1 ) {
                    return;
                }
            }
        }
    }

    private int[] parseGenotypes(byte genotype) {
        int[] genotypes = new int[4];
        genotypes[0] = ( genotype & 3 );
        genotypes[1] = ( ( genotype & 12 ) >>> 2 );
        genotypes[2] = ( ( genotype & 48 ) >>> 4 );
        genotypes[3] = ( ( genotype & 192 ) >>> 6 );
        return genotypes;
    }
}

class PlinkVariantInfo implements Comparable {

    private enum AlleleType {
        INSERTION,DELETION
    }

    private String variantName;
    private GenomeLoc loc;
    private ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
    private ArrayList<String> sampleNames = new ArrayList<String>();

    private ArrayList<Allele> indelHolder;
    private ArrayList<String> sampleHolder;
    private int siteIndelLength;
    private AlleleType indelType;

    // for binary parsing

    private String locAllele1;
    private String locAllele2;

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

    public boolean isSNP() {
        return this.indelType == null;
    }

    public void setGenomeLoc(GenomeLoc loc) {
        this.loc = loc;
    }

    public void setAlleles(String al1, String al2) {
        if ( al1.equals("0") ) {
            // encoding for a site at which no variants were detected
            locAllele1 = al2;
        } else {
            locAllele1 = al1;
        }
        locAllele2 = al2;
        if ( ! isSNP() ) {
            siteIndelLength = Math.max(locAllele1.length(),locAllele2.length());
        }

    }

    // CONSTRUCTOR

    public PlinkVariantInfo(String variantName) {
        this.variantName = variantName;
        parseName();
    }

    public PlinkVariantInfo(String variantName, boolean onlyLookForIndelInfo ) {
        this.variantName = variantName;
        if ( onlyLookForIndelInfo ) {
            parseNameForIndels();
        } else {
            parseName();
        }
    }

    private void parseName() {
        String chrom = this.variantName.split("\\|c")[1].split("_")[0];
        String pos = this.variantName.split("_p")[1].split("_")[0];
        this.loc = GenomeLocParser.parseGenomeLoc(chrom+":"+pos);
        this.parseNameForIndels();
    }

    private void parseNameForIndels() {
        if ( this.variantName.contains("_gI") || this.variantName.contains("_gD") ) {
            this.instantiateIndel();
        }
    }

    private void instantiateIndel() {
        this.indelHolder = new ArrayList<Allele>();
        this.sampleHolder = new ArrayList<String>();
        this.siteIndelLength = -1;
        this.indelType = this.variantName.contains("_gI") ? AlleleType.INSERTION : AlleleType.DELETION;
    }

    public void addGenotypeEntry(String genotypeString, String sampleName) {
        // identify if we're dealing with a deletion
        if ( this.isSNP() ) {
            this.addSNP(genotypeString.split("\\s+"),sampleName);
        } else {
            this.addIndel(genotypeString.split("\\s+"),sampleName);
        }
    }

    public void addBinaryGenotypeEntry( int genoTYPE, String sampleName ) {
        String[] alleleStr = new String[2];
        if ( genoTYPE == 0 ) {
            alleleStr[0] = locAllele1;
            alleleStr[1] = locAllele1;
        } else if (genoTYPE == 2) {
            alleleStr[0] = locAllele1;
            alleleStr[1] = locAllele2;
        } else if (genoTYPE == 3 ) {
            alleleStr[0] = locAllele2;
            alleleStr[1] = locAllele2;
        } else {
            alleleStr[0] = "0";
            alleleStr[1] = "0";
        }

        if ( this.isSNP() ) {
            this.addSNP(alleleStr,sampleName);
        } else {
            this.addIndel(alleleStr,sampleName);
        }
    }

    private void addSNP(String[] alleleStrings, String sampleName) {
        ArrayList<Allele> alleles = new ArrayList<Allele>(2);

        for ( String alStr : alleleStrings ) {
            alleles.add(new Allele(alStr.getBytes()));
        }

        genotypes.add(new Genotype(alleles,sampleName,20.0) );
        sampleNames.add(sampleName);
    }

    private void addIndel(String[] alleleStrings, String sampleName) {
        String alleleStr1 = alleleStrings[0];
        String alleleStr2 = alleleStrings[1];
        if ( alleleStr1.contains("-") ^ alleleStr2.contains("-") ) {
            // heterozygous indel
            if ( alleleStr1.contains("-") ) {
                this.addHetIndel(alleleStr2,sampleName) ;
            } else {
                this.addHetIndel(alleleStr1,sampleName);
            }
        } else {
            this.addHomIndel(alleleStr1, alleleStr2, sampleName);
        }
    }

    private void addHetIndel(String baseStr, String sampleName) {
        Allele ref;
        Allele alt;

        if ( indelType == AlleleType.INSERTION ) {
            ref = new Allele("-",true);
            alt = new Allele(baseStr.getBytes(),false);
        } else {
            alt = new Allele("-",false);
            ref = new Allele(baseStr.getBytes(),true);
        }

        // this.setIndelLength(alt,baseStr.length());

        if ( ! indelHolder.isEmpty() ) {
            siteIndelLength = baseStr.length();
            this.addHeldIndels();
        }

        Genotype indel = new Genotype(Arrays.asList(ref,alt), sampleName, 20.0);
        // this.setIndelGenotypeLength(indel,siteIndelLength);
        this.genotypes.add(indel);
        this.sampleNames.add(sampleName);
    }

    private void addHomIndel(String strand1, String strand2, String sampleName) {
        Allele allele1;
        Allele allele2;
        boolean reference;
        if ( indelType == AlleleType.DELETION ) {
            if ( strand1.contains("-") ) {
                // homozygous deletion
                allele1 = new Allele("-",false);
                allele2 = new Allele("-",false);
                reference = false;
            } else { // homozygous reference at a deletion variant site
                allele1 = new Allele(strand1.getBytes(),true);
                allele2 = new Allele(strand2.getBytes(),true);
                reference = true;
            }
        } else {
            if ( strand1.contains("-") ) {
                // homozygous reference
                allele1 = new Allele("-",true);
                allele2 = new Allele("-",true);
                reference = true;
            } else {
                allele1 = new Allele(strand1.getBytes(),false);
                allele2 = new Allele(strand2.getBytes(),false);
                reference = false;
            }
        }

        if ( reference ) {
            if ( ! indelHolder.isEmpty() ) {
                siteIndelLength = strand1.length();
                this.addHeldIndels();
            }
        }

        if ( reference || siteIndelLength != -1 ) { // if we're ref or know the insertion/deletion length of the site
            Genotype gen = new Genotype(Arrays.asList(allele1,allele2), sampleName, 20.0);
            // setIndelGenotypeLength(gen,siteIndelLength);
            this.genotypes.add(gen);
            this.sampleNames.add(sampleName);
        } else { // hold on the variants until we *do* know the in/del length at this site
            this.indelHolder.add(allele1);
            this.indelHolder.add(allele2);
            this.sampleHolder.add(sampleName);
        }

    }

/*    private void setIndelGenotypeLength(Genotype g, int length) {
        g.setAttribute(Genotype.StandardAttributes.DELETION_LENGTH,length);
    }

    private void addHeldIndels() {
        Allele del1;
        Allele del2;
        int startingSize = indelHolder.size();
        for ( int i = 0; i < startingSize ; i+=2 ) {
            del1 = indelHolder.get(i);
            del2 = indelHolder.get(i+1);
            this.addIndelFromCache(del1,del2,sampleHolder.get(i/2));
            if ( indelHolder.size() != startingSize ) {
                throw new StingException("Halting algorithm -- possible infinite loop");
            }
        }
        indelHolder.clear();
        sampleHolder.clear();
    }

    private void addIndelFromCache ( Allele indel1, Allele indel2, String sampleName ) {
        this.setIndelLength(indel1,siteIndelLength);
        this.setIndelLength(indel2,siteIndelLength);
        Genotype indel = new Genotype(Arrays.asList(indel1,indel2),sampleName, 20.0);
        // this.setIndelGenotypeLength(indel,siteIndelLength);
        this.genotypes.add(indel);
        this.sampleNames.add(sampleName);
    }

    public int compareTo(Object obj) {
        if ( ! ( obj instanceof PlinkVariantInfo) ) {
            return 1;
        }

        return loc.compareTo(((PlinkVariantInfo) obj).getLocation());
    }

    private void setIndelLength(Allele al, int length) {
        // Todo -- once alleles support deletion lengths add that information
        // Todo -- into the object; for now this can just return
        return;
    }
}

class PlinkBinaryTrifecta {

    public PlinkBinaryTrifecta() {

    }

    public File bedFile;
    public File bimFile;
    public File famFile;

} */