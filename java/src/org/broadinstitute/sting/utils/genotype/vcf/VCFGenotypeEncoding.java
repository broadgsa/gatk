package org.broadinstitute.sting.utils.genotype.vcf;


/**
 * @author aaron
 *         <p/>
 *         Class VCFGenotypeEncoding
 *         <p/>
 *         basic encoding class for genotype fields in VCF
 */
public class VCFGenotypeEncoding {
    public enum TYPE {
        SINGLE_BASE,
        INSERTION,
        DELETION,
        UNCALLED
    }

    // our length (0 for SINGLE_BASE), our bases, and our type
    private final int mLength;
    private final String mBases;
    private final TYPE mType;

    // public constructor, that parses out the base string
    public VCFGenotypeEncoding(String baseString) {
        if ((baseString.length() == 1)) {
            // are we an empty (no-call) genotype?
            if (baseString.equals(VCFGenotypeRecord.EMPTY_ALLELE)) {
                mBases = VCFGenotypeRecord.EMPTY_ALLELE;
                mLength = 0;
                mType = TYPE.UNCALLED;
            } else if (!validBases(baseString)) {
                throw new IllegalArgumentException("Alleles of length 1 must be one of A,C,G,T, " + baseString + " was passed in");
            } else { // we're a valid base
                mBases = baseString.toUpperCase();
                mLength = 0;
                mType = TYPE.SINGLE_BASE;
            }
        } else { // deletion or insertion
            if (baseString.length() < 1 || (baseString.toUpperCase().charAt(0) != 'D' && baseString.toUpperCase().charAt(0) != 'I')) {
                throw new IllegalArgumentException("Genotype encoding of " + baseString + " was passed in, but is not a valid deletion, insertion, base, or no call (.)");
            }
            if (baseString.toUpperCase().charAt(0) == 'D') {
                mLength = Integer.valueOf(baseString.substring(1, baseString.length()));
                mBases = "";
                mType = TYPE.DELETION;
            } else { // we're an I
                mBases = baseString.substring(1, baseString.length()).toUpperCase();
                if (!validBases(mBases))
                    throw new IllegalArgumentException("The insertion base string contained invalid bases -> " + baseString);
                mLength = mBases.length();
                mType = TYPE.INSERTION;
            }
        }
    }

    public int getLength() {
        return mLength;
    }

    public String getBases() {
        return mBases;
    }

    public TYPE getType() {
        return mType;
    }

    public boolean equals(Object obj) {
        if ( obj == null )
            return false;
        if ( obj instanceof VCFGenotypeEncoding ) {
            VCFGenotypeEncoding d = (VCFGenotypeEncoding) obj;
            return (mType == d.mType) && (mBases.equals(d.mBases)) && (mLength == d.mLength);
        }
        if ( mType == TYPE.UNCALLED && obj.toString().equals(VCFGenotypeRecord.EMPTY_ALLELE) )
            return true;
        return false;
    }

    public int hashCode() {
        // our underlying data is immutable, so this is safe (we won't strand a value in a hashtable somewhere
        // when the data changes underneath, altering this value).
        String str = this.mBases + String.valueOf(this.mLength) + this.mType.toString();
        return str.hashCode();
    }

    /**
     * dump the string representation of this genotype encoding
     *
     * @return string representation
     */
    public String toString() {
        StringBuilder builder = new StringBuilder();
        switch (mType) {
            case SINGLE_BASE:
            case UNCALLED:
                builder.append(mBases);
                break;
            case INSERTION:
                builder.append("I");
                builder.append(mBases);
                break;
            case DELETION:
                builder.append("D");
                builder.append(mLength);
                break;
        }
        return builder.toString();
    }

    /**
     * ensure that string contains valid bases
     *
     * @param bases the bases to check
     *
     * @return true if they're all either A,C,G,T; false otherwise
     */
    private static boolean validBases(String bases) {
        for (char c : bases.toUpperCase().toCharArray()) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                return false;
        }
        return true;
    }
}