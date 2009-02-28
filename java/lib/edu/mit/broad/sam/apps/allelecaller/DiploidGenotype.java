package edu.mit.broad.sam.apps.allelecaller;

public enum DiploidGenotype {
    AA('A','A'),
    AC('A','C'),
    AG('A','G'),
    AT('A','T'),
    CC('C','C'),
    CG('C','G'),
    CT('C','T'),
    GG('G','G'),
    GT('G','T'),
    TT('T','T');

    private final char allele1;
    private final char allele2;

    private DiploidGenotype(final char allele1, final char allele2) {
        this.allele1 = allele1;
        this.allele2 = allele2;
    }

    public char getAllele1() { return allele1; }
    public char getAllele2() { return allele2; }
    public boolean isHet() { return this.allele1 != this.allele2; }
    public boolean isHom() { return this.allele1 == this.allele2; }
}
