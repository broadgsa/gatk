package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.util.HashMap;

/**
 * A simple {codon -> amino acid name} lookup table.
 * Handles differences between mitochondrial and nuclear genes.
 */
public class AminoAcidTable {


    protected static final AminoAcid ISOLEUCINE = new AminoAcid("I" , "Isoleucine", "Ile");
    protected static final AminoAcid LEUCINE = new AminoAcid("L" , "Leucine", "Leu");
    protected static final AminoAcid VALINE =    new AminoAcid("V" , "Valine", "Val");
    protected static final AminoAcid PHENYLALANINE =    new AminoAcid("F" , "Phenylalanine", "Phe");
    protected static final AminoAcid METHIONINE = new AminoAcid("M" , "Methionine", "Met");
    protected static final AminoAcid CYSTEINE = new AminoAcid("C" , "Cysteine", "Cys");
    protected static final AminoAcid ALANINE = new AminoAcid("A" , "Alanine", "Ala");
    protected static final AminoAcid STOP_CODON = new AminoAcid("*" , "Stop Codon", "Stop");
    protected static final AminoAcid GLYCINE = new AminoAcid("G" , "Glycine", "Gly");
    protected static final AminoAcid PROLINE = new AminoAcid("P" , "Proline", "Pro");
    protected static final AminoAcid THEONINE = new AminoAcid("T" , "Threonine", "Thr");
    protected static final AminoAcid SERINE = new AminoAcid("S" , "Serine", "Ser");
    protected static final AminoAcid TYROSINE = new AminoAcid("Y" , "Tyrosine", "Tyr");
    protected static final AminoAcid TRYPTOPHAN = new AminoAcid("W" , "Tryptophan", "Trp");
    protected static final AminoAcid GLUTAMINE = new AminoAcid("Q" , "Glutamine", "Gln");
    protected static final AminoAcid ASPARAGINE = new AminoAcid("N" , "Asparagine", "Asn");
    protected static final AminoAcid HISTIDINE = new AminoAcid("H" , "Histidine", "His");
    protected static final AminoAcid GLUTAMIC_ACID = new AminoAcid("E" , "Glutamic acid", "Glu");
    protected static final AminoAcid ASPARTIC_ACID = new AminoAcid("D" , "Aspartic acid", "Asp");
    protected static final AminoAcid LYSINE = new AminoAcid("K" , "Lysine", "Lys");
    protected static final AminoAcid ARGININE = new AminoAcid("R" , "Arginine", "Arg");

    protected static HashMap<String, AminoAcid> aminoAcidTable = new HashMap<String, AminoAcid>();
    protected static HashMap<String, AminoAcid> mitochondrialAminoAcidTable = new HashMap<String, AminoAcid>();

    static {
        //populate the tables
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


        //populate the mitochondrial AA table
        mitochondrialAminoAcidTable.putAll(aminoAcidTable);
        mitochondrialAminoAcidTable.put("AGA", STOP_CODON);
        mitochondrialAminoAcidTable.put("AGG", STOP_CODON);
        mitochondrialAminoAcidTable.put("ATA", METHIONINE);
        mitochondrialAminoAcidTable.put("TGA", TRYPTOPHAN);
    }




    /**
     * Returns the the matching amino acid.
     *
     * @param codon The 3-letter mRNA nucleotide codon 5' to 3'. Expects T's instead of U's. Not case sensitive.
     * @param mitochondrial Whether this is from a mitochondrial gene (mitochondria have a slightly different codon table).
     *
     * @return The amino acid matching the given codon.
     */
    public static AminoAcid getAA(String codon, boolean mitochondrial) {
        codon = codon.toUpperCase();
        if(!aminoAcidTable.containsKey(codon)) {
            throw new IllegalArgumentException("Invalid codon: " + codon);
        }

        return aminoAcidTable.get(codon);
    }
}
