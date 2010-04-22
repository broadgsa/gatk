package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.util.HashMap;

/**
 * A simple {codon -> amino acid name} lookup table.
 * Handles differences between mitochondrial and nuclear genes.
 */
public class AminoAcidTable {

    protected static HashMap<String, String[]> aminoAcidTable = new HashMap<String, String[]>();
    protected static HashMap<String, String[]> mtAminoAcidTable = new HashMap<String, String[]>(); //mitochondrial

    protected static final String[] ISOLEUCINE = new String[] {"I" , "Isoleucine", "Ile"};
    protected static final String[] LEUCINE = new String[] {"L" , "Leucine", "Leu"};

    protected static final String[] VALINE =    new String[] {"V" , "Valine", "Val"};
    protected static final String[] PHENYLALANINE =    new String[] {"F" , "Phenylalanine", "Phe"};
    protected static final String[] METHIONINE = new String[] {"M" , "Methionine", "Met"};
    protected static final String[] CYSTEINE = new String[] {"C" , "Cysteine", "Cys"};
    protected static final String[] ALANINE = new String[] {"A" , "Alanine", "Ala"};
    protected static final String[] STOP_CODON = new String[] {"*" , "Stop Codon", "Amb"};
    protected static final String[] GLYCINE = new String[] {"G" , "Glycine", "Gly"};
    protected static final String[] PROLINE = new String[] {"P" , "Proline", "Pro"};
    protected static final String[] THEONINE = new String[] {"T" , "Threonine", "Thr"};
    protected static final String[] SERINE = new String[] {"S" , "Serine", "Ser"};
    protected static final String[] TYROSINE = new String[] {"Y" , "Tyrosine", "Tyr"};
    protected static final String[] TRYPTOPHAN = new String[] {"W" , "Tryptophan", "Trp"};
    protected static final String[] GLUTAMINE = new String[] {"Q" , "Glutamine", "Gln"};
    protected static final String[] ASPARAGINE = new String[] {"N" , "Asparagine", "Asn"};
    protected static final String[] HISTIDINE = new String[] {"H" , "Histidine", "His"};
    protected static final String[] GLUTAMIC_ACID = new String[] {"E" , "Glutamic acid", "Glu"};
    protected static final String[] ASPARTIC_ACID = new String[] {"D" , "Aspartic acid", "Asp"};
    protected static final String[] LYSINE = new String[] {"K" , "Lysine", "Lys"};
    protected static final String[] ARGININE = new String[] {"R" , "Arginine", "Arg"};

    static {

        //code gotten from: http://algoart.com/aatable.htm
        aminoAcidTable.put("ATT", ISOLEUCINE);
        aminoAcidTable.put("ATC", ISOLEUCINE);
        aminoAcidTable.put("ATA", ISOLEUCINE);


        aminoAcidTable.put("CTT", LEUCINE);
        aminoAcidTable.put("CTC", LEUCINE);
        aminoAcidTable.put("CTA", LEUCINE);
        aminoAcidTable.put("CTG", LEUCINE);
        aminoAcidTable.put("TTA", LEUCINE);
        aminoAcidTable.put("TTG", LEUCINE);


        aminoAcidTable.put("GTT", VALINE);
        aminoAcidTable.put("GTC", VALINE);
        aminoAcidTable.put("GTA", VALINE);
        aminoAcidTable.put("GTG", VALINE);


        aminoAcidTable.put("TTT", PHENYLALANINE);
        aminoAcidTable.put("TTC", PHENYLALANINE);


        aminoAcidTable.put("ATG", METHIONINE);

        aminoAcidTable.put("TGT", CYSTEINE);
        aminoAcidTable.put("TGC", CYSTEINE);

        aminoAcidTable.put("GCT", ALANINE);
        aminoAcidTable.put("GCC", ALANINE);
        aminoAcidTable.put("GCA", ALANINE);
        aminoAcidTable.put("GCG", ALANINE);


        aminoAcidTable.put("GGT", GLYCINE);
        aminoAcidTable.put("GGC", GLYCINE);
        aminoAcidTable.put("GGA", GLYCINE);
        aminoAcidTable.put("GGG", GLYCINE);


        aminoAcidTable.put("CCT", PROLINE);
        aminoAcidTable.put("CCC", PROLINE);
        aminoAcidTable.put("CCA", PROLINE);
        aminoAcidTable.put("CCG", PROLINE);




        aminoAcidTable.put("ACT", THEONINE);
        aminoAcidTable.put("ACC", THEONINE);
        aminoAcidTable.put("ACA", THEONINE);
        aminoAcidTable.put("ACG", THEONINE);



        aminoAcidTable.put("TCT", SERINE);
        aminoAcidTable.put("TCC", SERINE);
        aminoAcidTable.put("TCA", SERINE);
        aminoAcidTable.put("TCG", SERINE);
        aminoAcidTable.put("AGT", SERINE);
        aminoAcidTable.put("AGC", SERINE);

        aminoAcidTable.put("TAT", TYROSINE);
        aminoAcidTable.put("TAC", TYROSINE);



        aminoAcidTable.put("TGG", TRYPTOPHAN);


        aminoAcidTable.put("CAA", GLUTAMINE);
        aminoAcidTable.put("CAG", GLUTAMINE);


        aminoAcidTable.put("AAT", ASPARAGINE);
        aminoAcidTable.put("AAC", ASPARAGINE);


        aminoAcidTable.put("CAT", HISTIDINE);
        aminoAcidTable.put("CAC", HISTIDINE);


        aminoAcidTable.put("GAA", GLUTAMIC_ACID);
        aminoAcidTable.put("GAG", GLUTAMIC_ACID);



        aminoAcidTable.put("GAT", ASPARTIC_ACID);
        aminoAcidTable.put("GAC", ASPARTIC_ACID);


        aminoAcidTable.put("AAA", LYSINE);
        aminoAcidTable.put("AAG", LYSINE);


        aminoAcidTable.put("CGT", ARGININE);
        aminoAcidTable.put("CGC", ARGININE);
        aminoAcidTable.put("CGA", ARGININE);
        aminoAcidTable.put("CGG", ARGININE);
        aminoAcidTable.put("AGA", ARGININE);
        aminoAcidTable.put("AGG", ARGININE);


        aminoAcidTable.put("TAA", STOP_CODON );
        aminoAcidTable.put("TAG", STOP_CODON);
        aminoAcidTable.put("TGA", STOP_CODON);


        //populate the mitochondrial ammino-acid table:
        mtAminoAcidTable.putAll(aminoAcidTable);
        mtAminoAcidTable.put("AGA", STOP_CODON);
        mtAminoAcidTable.put("AGG", STOP_CODON);
        mtAminoAcidTable.put("ATA", METHIONINE);
        mtAminoAcidTable.put("TGA", TRYPTOPHAN);
    }



    /**
     * Returns the 1-letter code for the matching amino acid.
     *
     * @param codon The 3-letter mRNA nucleotide codon 5' to 3'. Expects T's instead of U's. Not case sensitive.
     * @param mitochondrial Whether this is from a mitochondrial gene (mitochondria have a slightly different codon table).
     *
     * @return The 1 letter code of the amino acid.
     */
    public static String getAA1LetterCode(String codon, boolean mitochondrial) {
        codon = codon.toUpperCase();
        if(!aminoAcidTable.containsKey(codon)) {
            throw new IllegalArgumentException("Invalid codon: " + codon);
        }

        return aminoAcidTable.get(codon)[0];
    }

    /**
     * Returns the 3-letter code for the matching amino acid.
     *
     * @param codon The 3-letter mRNA nucleotide codon 5' to 3'. Expects T's instead of U's. Not case sensitive.
     * @param mitochondrial Whether this is from a mitochondrial gene (mitochondria have a slightly different codon table).
     *
     * @return The 3 letter code of the amino acid.
     */
    public static String getAA3LetterCode(String codon, boolean mitochondrial) {
        codon = codon.toUpperCase();
        if(!aminoAcidTable.containsKey(codon)) {
            throw new IllegalArgumentException("Invalid codon: " + codon);
        }

        return aminoAcidTable.get(codon)[2];
    }


    /**
     * Returns the full name of the matching amino acid.
     *
     * @param codon The 3-letter mRNA nucleotide codon 5' to 3'. Expects T's instead of U's. Not case sensitive.
     * @param mitochondrial Whether this is from a mitochondrial gene (mitochondria have a slightly different codon table).
     *
     * @return The full name of the amino acid.
     */
    public static String getAAName(String codon, boolean mitochondrial) {
        codon = codon.toUpperCase();
        if(!aminoAcidTable.containsKey(codon)) {
            throw new IllegalArgumentException("Invalid codon: " + codon);
        }

        return aminoAcidTable.get(codon)[1];
    }
}
