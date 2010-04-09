package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl, ebanks
 *
 */
public class PlinkRod extends BasicReferenceOrderedDatum implements Iterator<PlinkRod> {

    public static final String SEQUENOM_NO_CALL = "0";
    public static final String SEQUENOM_NO_BASE = "-";

    private final Set<String> headerEntries = new HashSet<String>(Arrays.asList("#Family ID","Individual ID","Sex",
                "Paternal ID","Maternal ID","Phenotype", "FID","IID","PAT","MAT","SEX","PHENOTYPE"));
    private final byte SNP_MAJOR_MODE = 1;

    private PlinkVariantInfo currentVariant;
    private ListIterator<PlinkVariantInfo> variantIterator;

    private PlinkFileType plinkFileType;

    private List<String> sampleNames;

    public enum PlinkFileType {
        STANDARD_PED, RAW_PED, BINARY_PED
    }

    public PlinkRod(String name) {
        super(name);
    }

    public PlinkRod(String name, PlinkVariantInfo record, ListIterator<PlinkVariantInfo> iter, List<String> sampleNames) {
        super(name);
        currentVariant = record;
        variantIterator = iter;
        this.sampleNames = sampleNames;
    }

    @Override
    public Object initialize(final File plinkFile) throws FileNotFoundException {
        if ( ! plinkFile.exists() ) {
            throw new FileNotFoundException("File "+plinkFile.getAbsolutePath()+" does not exist.");
        }

        ArrayList<PlinkVariantInfo> variants = parsePlinkFile(plinkFile);
        if ( variants == null )
            throw new IllegalStateException("Error parsing Plink file");

        variantIterator = variants.listIterator();
        return null;
    }

    public static PlinkRod createIterator(String name, File file) {
        PlinkRod plink = new PlinkRod(name);
        try {
            plink.initialize(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to find file " + file);
        }
        return plink;
    }

    private void assertNotNull() {
        if ( currentVariant == null ) {
            throw new IllegalStateException("Current variant information is null");
        }
    }

    public boolean hasNext() {
        return variantIterator.hasNext();
    }

    public PlinkRod next() {
        if ( !this.hasNext() )
            throw new NoSuchElementException("PlinkRod next called on iterator with no more elements");

        // get the next record
        currentVariant = variantIterator.next();
        return new PlinkRod(name, currentVariant, variantIterator, sampleNames);
    }

    @Override
    public boolean parseLine(Object obj, String[] args) {
        throw new UnsupportedOperationException("PlinkRod does not support the parseLine method");
    }

    public void remove() {
        throw new UnsupportedOperationException("The remove operation is not supported for a PlinkRod");
    }

    @Override
    public GenomeLoc getLocation() {
        assertNotNull();
        return currentVariant.getLocation();
    }

    @Override
    public String toString() {
        assertNotNull();
        return currentVariant.toString();
    }

    public String getVariantName() {
        assertNotNull();
        return currentVariant.getName();

    }

    public Map<String, List<byte[]>> getGenotypes() {
        return currentVariant.getGenotypes();
    }

    public List<String> getSampleNames() {
        return sampleNames;
    }

    public boolean isIndel() {
        assertNotNull();
        return currentVariant.isIndel();
    }

    public boolean isInsertion() {
        assertNotNull();
        return currentVariant.isInsertion();
    }

    public int getLength() {
        assertNotNull();
        return currentVariant.getLength();
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

    private ArrayList<PlinkVariantInfo> parseTextFormattedPlinkFile( File file ) {
        try {
            BufferedReader reader = new BufferedReader( new FileReader ( file ) );
            String header = reader.readLine();
            ArrayList<PlinkVariantInfo> seqVars = instantiateVariantListFromHeader(header);
            ArrayList<Integer> snpOffsets = getSNPOffsetsFromHeader(header);

            sampleNames = new ArrayList<String>();

            String line;
            long counter = 0;
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
        if ( plinkFileType != PlinkFileType.STANDARD_PED )
            throw new StingException("Plink file is likely of .raw or recoded format. Please use an uncoded .ped file.");

        plinkInfo = plinkLine.split("\t");
        String individualName = plinkInfo[1];
        sampleNames.add(individualName);

        int snpNumber = 0;
        for ( int i : offsets ) {
            vars.get(snpNumber).addGenotypeEntry(plinkInfo[i].split("\\s+"), individualName);
            snpNumber++;
        }
    }

    private ArrayList<PlinkVariantInfo> instantiateVariantListFromHeader(String header) {
        // if the first line is not a comment (what we're used to seeing),
        // then it's the raw header (comes from de-binary-ing a .bed file)
        if ( !header.startsWith("#") )
            throw new StingException("Plink file is likely of .raw or recoded format. Please use an uncoded .ped file.");

        plinkFileType = PlinkFileType.STANDARD_PED;

        ArrayList<PlinkVariantInfo> seqVars = new ArrayList<PlinkVariantInfo>();
        String[] headerFields = header.split("\t");

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

    private ArrayList<PlinkVariantInfo> parseBinaryFormattedPlinkFile(File file) {
        PlinkBinaryTrifecta binaryFiles = getPlinkBinaryTriplet(file);
        ArrayList<PlinkVariantInfo> parsedVariants = instantiateVariantsFromBimFile(binaryFiles.bimFile);
        sampleNames = getSampleNameOrderingFromFamFile(binaryFiles.famFile);
        ArrayList<PlinkVariantInfo> updatedVariants = getGenotypesFromBedFile(parsedVariants, binaryFiles.bedFile);

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
            String line = reader.readLine();
            while ( line != null ) {
                String[] snpInfo = line.split("\\s+");
                PlinkVariantInfo variant = new PlinkVariantInfo(snpInfo[1]);
                variant.setGenomeLoc(GenomeLocParser.parseGenomeLoc(snpInfo[0],Long.valueOf(snpInfo[3]), Long.valueOf(snpInfo[3])));
                variant.setAlleles(snpInfo[4],snpInfo[5]);
                variants.add(variant);

                line = reader.readLine();
            }
        } catch ( IOException e ) {
            throw new StingException("There was an error reading the .bim file "+bimFile.getAbsolutePath(), e);
        }

        return variants;
    }

    private static ArrayList<String> getSampleNameOrderingFromFamFile(File famFile) {
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
            String line;
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

    private ArrayList<PlinkVariantInfo> getGenotypesFromBedFile(ArrayList<PlinkVariantInfo> variants, File bedFile) {
        FileInputStream inStream;
        try {
            inStream = new FileInputStream(bedFile);
        } catch (FileNotFoundException e) {
            throw new StingException("The Binary pedigree file file accompanying the family file was not found (the .bed file). "+
                                     "Please ensure that it is in the same directory as the .bim and .fam files. The file we "+
                                     "Were looking for was "+bedFile.getAbsolutePath(),e);
        }

        try {
            byte genotype;
            long bytesRead = 0;
            int snpOffset = 0;
            int sampleOffset = 0;
            boolean snpMajorMode = true;
            do {
                genotype = (byte) inStream.read();
                bytesRead++;
                if ( genotype != -1 ) {
                    if ( bytesRead > 3 ) {
                        addGenotypeByte(genotype,variants,snpOffset,sampleOffset, snpMajorMode);

                        if ( snpMajorMode ) {
                            sampleOffset = sampleOffset + 4;
                            if ( sampleOffset > sampleNames.size() -1 ) {
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

    private void addGenotypeByte(byte genotype, ArrayList<PlinkVariantInfo> variants, int snpOffset, int sampleOffset, boolean major) {
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

    private String variantName;
    private GenomeLoc loc;
    private Map<String, List<byte[]>> genotypes = new HashMap<String, List<byte[]>>();

    // for indels
    private boolean isIndel = false;
    private boolean isInsertion = false;
    private int length = 1;

    // for binary parsing
    private String locAllele1;
    private String locAllele2;


    public PlinkVariantInfo(String variantName) {
        this.variantName = variantName;
        parseName();
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public String getName() {
        return variantName;
    }

    public Map<String, List<byte[]>> getGenotypes() {
        return genotypes;
    }

    public boolean isIndel() {
        return isIndel;
    }

    public boolean isInsertion() {
        return isInsertion;
    }

    public int getLength() {
        return length;
    }

    public void setGenomeLoc(GenomeLoc loc) {
        this.loc = loc;
    }

    private void parseName() {
        int chromIdx = variantName.indexOf("|c");
        if ( chromIdx == -1 )
            throw new IllegalArgumentException("Variant name " + variantName + " does not adhere to required convention (...|c...)");
        String[] pieces = variantName.substring(chromIdx+2).split("_");
        if ( pieces.length < 2 )
            throw new IllegalArgumentException("Variant name " + variantName + " does not adhere to required convention (...|c..._p...)");

        String chrom = pieces[0].trim();
        if ( pieces[1].charAt(0) != 'p' )
            throw new IllegalArgumentException("Variant name " + variantName + " does not adhere to required convention (...|c..._p...)");

        String pos = pieces[1].substring(1).trim();
        loc = GenomeLocParser.parseGenomeLoc(chrom+":"+pos);

        if ( pieces.length > 2 && (pieces[2].startsWith("gI") || pieces[2].startsWith("gD")) ) {
            // it's an indel
            isIndel = true;
            isInsertion = pieces[2].startsWith("gI");
            try {
                // length of insertion on reference is still 1
                if ( !isInsertion )
                    length = Integer.parseInt(pieces[2].substring(2));
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Variant name " + variantName + " does not adhere to required convention (...|c..._p..._g[I/D][length])");
            }
        }
    }

    public void setAlleles(String al1, String al2) {
        if ( al1.equals(PlinkRod.SEQUENOM_NO_CALL) ) {
            // encoding for a site at which no variants were detected
            locAllele1 = al2;
        } else {
            locAllele1 = al1;
        }
        locAllele2 = al2;
    }

    public void addGenotypeEntry(String[] alleleStrings, String sampleName) {

        ArrayList<byte[]> alleles = new ArrayList<byte[]>(2);

        for ( String alleleString : alleleStrings ) {
            if ( alleleString.equals(PlinkRod.SEQUENOM_NO_CALL) )
                alleles.add(Allele.NO_CALL_STRING.getBytes());
            else
                alleles.add(alleleString.getBytes());
        }

        genotypes.put(sampleName, alleles);
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

        addGenotypeEntry(alleleStr, sampleName);
    }

    public int compareTo(Object obj) {
        if ( ! ( obj instanceof PlinkVariantInfo) ) {
            return 1;
        }

        return loc.compareTo(((PlinkVariantInfo) obj).getLocation());
    }
}

class PlinkBinaryTrifecta {

    public PlinkBinaryTrifecta() {}

    public File bedFile;
    public File bimFile;
    public File famFile;

}