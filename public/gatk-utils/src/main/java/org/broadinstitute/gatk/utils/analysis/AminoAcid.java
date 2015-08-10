/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.analysis;

/*
 * Copyright (c) 2010 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author chartl
 * @since June 28, 2010
 */

public enum AminoAcid {
    
    Alanine("Alanine","Ala","A",new String[]{"GCA","GCC","GCG","GCT"}),
    Arganine("Arganine","Arg","R",new String[]{"AGA","AGG","CGA","CGC","CGG","CGT"}),
    Asparagine("Asparagine","Asn","N",new String[]{"AAC","AAT"}),
    Aspartic_acid("Aspartic acid","Asp","D",new String[]{"GAT","GAC"}),
    Cysteine("Cysteine","Cys","C",new String[]{"TGC","TGC"}),
    Glutamic_acid("Glutamic acid","Glu","E",new String[]{"GAA","GAG"}),
    Glutamine("Glutamine","Gln","Q",new String[]{"CAA","CAG"}),
    Glycine("Glycine","Gly","G",new String[]{"GGA","GGC","GGG","GGT"}),
    Histidine("Histidine","His","H",new String[]{"CAC","CAT"}),
    Isoleucine("Isoleucine","Ile","I",new String[]{"ATA","ATC","ATT"}),
    Leucine("Leucine","Leu","L",new String[]{"CTA","CTC","CTG","CTT","TTA","TTG"}),
    Lysine("Lysine","Lys","K", new String[]{"AAA","AAG"}),
    Methionine("Methionine","Met","M",new String[]{"ATG"}),
    Phenylalanine("Phenylalanine","Phe","F",new String[]{"TTC","TTT"}),
    Proline("Proline","Pro","P",new String[]{"CCA","CCC","CCG","CCT"}),
    Serine("Serine","Ser","S",new String[]{"AGC","AGT","TCA","TCC","TCG","TCT"}),
    Stop_codon("Stop codon","Stop","*",new String[]{"TAA","TAG","TGA"}),
    Threonine("Threonine","Thr","T",new String[]{"ACA","ACC","ACG","ACT"}),
    Tryptophan("Tryptophan","Trp","W",new String[]{"TGG"}),
    Tyrosine("Tyrosine","Tyr","Y",new String[]{"TAC","TAT"}),
    Valine("Valine","Val","V",new String[]{"GTA","GTC","GTG","GTT"});

    String[] codons;
    String fullName;
    String code;
    String letter;

    AminoAcid(String name, String shortName, String abbrev, String[] myCodons) {
        codons = myCodons;
        fullName = name;
        code = shortName;
        letter = abbrev;
    }

    public String getName() {
        return fullName;
    }

    public String getLetter() {
        return letter;
    }

    public String getCode() {
        return code;
    }

    public boolean isStop() {
        return this == Stop_codon;
    }

    public String toString() {
        return getName();
    }
    
}
