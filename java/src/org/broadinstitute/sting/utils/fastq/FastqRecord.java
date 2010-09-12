package org.broadinstitute.sting.utils.fastq;

public class FastqRecord {
    private String seqHeader;
    private String seqLine;
    private String qualHeader;
    private String qualLine;

    public FastqRecord(String seqHeader, String seqLine, String qualHeader, String qualLine) {
        setReadHeader(seqHeader);
        setReadString(seqLine);
        setBaseQualityHeader(qualHeader);
        setBaseQualityString(qualLine);
    }

    public void setReadHeader(String seqHeader) {
        this.seqHeader = seqHeader.replaceFirst("@", "");
    }

    public void setReadString(String seqLine) { this.seqLine = seqLine; }

    public void setBaseQualityHeader(String qualHeader) {
        this.qualHeader = qualHeader.replaceFirst("\\+", "");
    }

    public void setBaseQualityString(String qualLine) { this.qualLine = qualLine; }

    public String getReadHeader() { return seqHeader; }
    public String getReadString() { return seqLine; }
    public String getBaseQualityHeader() { return qualHeader; }
    public String getBaseQualityString() { return qualLine; }

    public String format() {
        return String.format("@%s\n%s\n+%s\n%s", seqHeader, seqLine, qualHeader, qualLine);
    }

    public String toString() {
        return String.format("%s : %s %s", seqHeader, seqLine, qualLine);
    }
}
