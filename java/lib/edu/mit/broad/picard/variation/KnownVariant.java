package edu.mit.broad.picard.variation;

/**
 * Utility class to hold data about a population or somatic variant.
 *
 * IMPORTANT!  Regardless of the coordinate system of the data from which it is drawn, the data
 * in this class should be 1-based.  Start and end coordinates should be as follows:
 *  For SNPs, start and end should be the same base.
 *  For insertions and deletions, the base on either side of the affected reference sequence
 *    will be the start and end.  For insertions, this means they will always be 1 base apart.
 */
public class KnownVariant implements Comparable<KnownVariant>
{
    private final String name;
    private final int sequenceIndex;
    private final int startPos;
    private final int endPos;
    private final VariantType type;
    private final boolean validated;
    private transient String referenceSequence;

    /**
     * Constructor
     * 
     * @param name
     * @param sequenceIndex
     * @param startPos
     * @param endPos
     * @param type
     * @param validated
     */
    public KnownVariant(String name, int sequenceIndex, int startPos, int endPos,
                        VariantType type, boolean validated)
    {
        this.name = name;
        this.sequenceIndex = sequenceIndex;
        this.startPos = startPos;
        this.endPos = endPos;
        this.type = type;
        this.validated = validated;
    }

    /**
     * Compares this object with the specified object for order. Returns a negative integer, zero, or a positive
     * integer as this object is less than, equal to, or greater than the specified object.
     *
     * @param that  The KnownVariant to compare
     * @return  a negative integer, zero, or a positive integer as this object is less than, equal to,
     *          or greater than the specified object
     */
    public int compareTo(KnownVariant that)
    {
        if (this.getSequenceIndex() != that.getSequenceIndex())
        {
            return (this.getSequenceIndex() > that.getSequenceIndex()) ? 1 : -1;
        }
        else if (this.getStartPos() != that.getStartPos())
        {
            return (this.getStartPos() > that.getStartPos()) ? 1 : -1;
        }
        else if (this.getEndPos() != that.getEndPos())
        {
            return (this.getEndPos() > that.getEndPos()) ? 1 : -1;
        }
        else if (!this.getName().equals(that.getName()))
        {
            return this.getName().compareTo(that.getName());
        }
        else if (this.getType() != that.getType())
        {
            return this.getType().compareTo(that.getType());
        }
        else if (this.isValidated() != that.isValidated())
        {
            return this.isValidated() ? 1 : -1;
        }
        return 0;
    }

    public boolean equals(Object o)
    {
        if (!(o instanceof KnownVariant)) {
            return false;
        }
        KnownVariant that = (KnownVariant)o;
        return (this.name.equals(that.name) &&
                this.sequenceIndex == that.sequenceIndex &&
                this.startPos == that.startPos &&
                this.endPos == that.endPos &&
                this.type == that.type &&
                this.validated == that.validated);
    }

    public int hasCode()
    {
        int result = 17;
        result = 37*result + name.hashCode();
        result = 37*result + sequenceIndex;
        result = 37*result + startPos;
        result = 37*result + endPos;
        result = 37*result + type.hashCode();
        result = 37*result + (validated ? 1 : 0);
        return result;
    }

    public String getName() { return name; }
    public int getSequenceIndex() { return sequenceIndex; }
    public String getRefrenceSequence() { return referenceSequence; }
    public void setRefrenceSequence(String referenceSequence) { this.referenceSequence = referenceSequence; }
    public int getStartPos() { return startPos; }
    public int getEndPos() { return endPos; }
    public VariantType getType() { return type; }
    public boolean isValidated() { return validated; }

}
