/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.snpEff;

import org.broad.tribble.Feature;

import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.EffectType;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.EffectImpact;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.ChangeType;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.Zygosity;

public class SnpEffFeature implements Feature {

    private String contig;
    private long position;
    private String reference;
    private String change;
    private ChangeType changeType;
    private Zygosity zygosity;
    private Double quality;
    private Long coverage;
    private String warnings;
    private String geneID;
    private String geneName;
    private String bioType;
    private String transcriptID;
    private String exonID;
    private Integer exonRank;
    private boolean isNonCodingGene;
    private EffectType effect;
    private String effectExtraInformation;
    private String oldAndNewAA;
    private String oldAndNewCodon;
    private Integer codonNum;
    private Integer cdsSize;
    private String codonsAround;
    private String aasAround;
    private String customIntervalID;

    public SnpEffFeature ( String contig,
                           long position,
                           String reference,
                           String change,
                           ChangeType changeType,
                           Zygosity zygosity,
                           Double quality,
                           Long coverage,
                           String warnings,
                           String geneID,
                           String geneName,
                           String bioType,
                           String transcriptID,
                           String exonID,
                           Integer exonRank,
                           boolean isNonCodingGene,
                           EffectType effect,
                           String effectExtraInformation,
                           String oldAndNewAA,
                           String oldAndNewCodon,
                           Integer codonNum,
                           Integer cdsSize,
                           String codonsAround,
                           String aasAround,
                           String customIntervalID ) {

        this.contig = contig;
        this.position = position;
        this.reference = reference;
        this.change = change;
        this.changeType = changeType;
        this.zygosity = zygosity;
        this.quality = quality;
        this.coverage = coverage;
        this.warnings = warnings;
        this.geneID = geneID;
        this.geneName = geneName;
        this.bioType = bioType;
        this.transcriptID = transcriptID;
        this.exonID = exonID;
        this.exonRank = exonRank;
        this.isNonCodingGene = isNonCodingGene;
        this.effect = effect;
        this.effectExtraInformation = effectExtraInformation;
        this.oldAndNewAA = oldAndNewAA;
        this.oldAndNewCodon = oldAndNewCodon;
        this.codonNum = codonNum;
        this.cdsSize = cdsSize;
        this.codonsAround = codonsAround;
        this.aasAround = aasAround;
        this.customIntervalID = customIntervalID;
    }

    public boolean isHigherImpactThan ( SnpEffFeature other ) {
        if ( ! isNonCodingGene() && other.isNonCodingGene() ) {
            return true;
        }
        else if ( isNonCodingGene() && ! other.isNonCodingGene() ) {
            return false;
        }

        return getEffectImpact().isHigherImpactThan(other.getEffectImpact());
    }

    public String getChr() {
        return contig;
    }

    public int getStart() {
        return (int)position;
    }

    public int getEnd() {
        return (int)position;
    }

    public boolean hasReference() {
        return reference != null;
    }

    public String getReference() {
        return reference;
    }

    public boolean hasChange() {
        return change != null;
    }

    public String getChange() {
        return change;
    }

    public boolean hasChangeType() {
        return changeType != null;
    }

    public ChangeType getChangeType() {
        return changeType;
    }

    public boolean hasZygosity() {
        return zygosity != null;
    }

    public Zygosity getZygosity() {
        return zygosity;
    }

    public boolean hasQuality() {
        return quality != null;
    }

    public Double getQuality() {
        return quality;
    }

    public boolean hasCoverage() {
        return coverage != null;
    }

    public Long getCoverage() {
        return coverage;
    }

    public boolean hasWarnings() {
        return warnings != null;
    }

    public String getWarnings() {
        return warnings;
    }

    public boolean hasGeneID() {
        return geneID != null;
    }

    public String getGeneID() {
        return geneID;
    }

    public boolean hasGeneName() {
        return geneName != null;
    }

    public String getGeneName() {
        return geneName;
    }

    public boolean hasBioType() {
        return bioType != null;
    }

    public String getBioType() {
        return bioType;
    }

    public boolean hasTranscriptID() {
        return transcriptID != null;
    }

    public String getTranscriptID() {
        return transcriptID;
    }

    public boolean hasExonID() {
        return exonID != null;
    }

    public String getExonID() {
        return exonID;
    }

    public boolean hasExonRank() {
        return exonRank != null;
    }

    public Integer getExonRank() {
        return exonRank;
    }

    public boolean isNonCodingGene() {
        return isNonCodingGene;
    }

    public EffectType getEffect() {
        return effect;
    }

    public EffectImpact getEffectImpact() {
        return effect.getImpact();
    }

    public boolean hasEffectExtraInformation() {
        return effectExtraInformation != null;
    }

    public String getEffectExtraInformation() {
        return effectExtraInformation;
    }

    public boolean hasOldAndNewAA() {
        return oldAndNewAA != null;
    }

    public String getOldAndNewAA() {
        return oldAndNewAA;
    }

    public boolean hasOldAndNewCodon() {
        return oldAndNewCodon != null;
    }

    public String getOldAndNewCodon() {
        return oldAndNewCodon;
    }

    public boolean hasCodonNum() {
        return codonNum != null;
    }

    public Integer getCodonNum() {
        return codonNum;
    }

    public boolean hasCdsSize() {
        return cdsSize != null;
    }

    public Integer getCdsSize() {
        return cdsSize;
    }

    public boolean hasCodonsAround() {
        return codonsAround != null;
    }

    public String getCodonsAround() {
        return codonsAround;
    }

    public boolean hadAasAround() {
        return aasAround != null;
    }

    public String getAasAround() {
        return aasAround;
    }

    public boolean hasCustomIntervalID() {
        return customIntervalID != null;
    }

    public String getCustomIntervalID() {
        return customIntervalID;
    }

    public boolean equals ( Object o ) {
        if ( o == null || ! (o instanceof SnpEffFeature) ) {
            return false;
        }

        SnpEffFeature other = (SnpEffFeature)o;

        return contig.equals(other.contig) &&
               position == other.position &&
               (reference == null ? other.reference == null : reference.equals(other.reference)) &&
               (change == null ? other.change == null : change.equals(other.change)) &&
               changeType == other.changeType &&
               zygosity == other.zygosity &&
               (quality == null ? other.quality == null : quality.equals(other.quality)) &&
               (coverage == null ? other.coverage == null : coverage.equals(other.coverage)) &&
               (warnings == null ? other.warnings == null : warnings.equals(other.warnings)) &&
               (geneID == null ? other.geneID == null : geneID.equals(other.geneID)) &&
               (geneName == null ? other.geneName == null : geneName.equals(other.geneName)) &&
               (bioType == null ? other.bioType == null : bioType.equals(other.bioType)) &&
               (transcriptID == null ? other.transcriptID == null : transcriptID.equals(other.transcriptID)) &&
               (exonID == null ? other.exonID == null : exonID.equals(other.exonID)) &&
               (exonRank == null ? other.exonRank == null : exonRank.equals(other.exonRank)) &&
               isNonCodingGene == other.isNonCodingGene &&
               effect == other.effect &&
               (effectExtraInformation == null ? other.effectExtraInformation == null : effectExtraInformation.equals(other.effectExtraInformation)) &&
               (oldAndNewAA == null ? other.oldAndNewAA == null : oldAndNewAA.equals(other.oldAndNewAA)) &&
               (oldAndNewCodon == null ? other.oldAndNewCodon == null : oldAndNewCodon.equals(other.oldAndNewCodon)) &&
               (codonNum == null ? other.codonNum == null : codonNum.equals(other.codonNum)) &&
               (cdsSize == null ? other.cdsSize == null : cdsSize.equals(other.cdsSize)) &&
               (codonsAround == null ? other.codonsAround == null : codonsAround.equals(other.codonsAround)) &&
               (aasAround == null ? other.aasAround == null : aasAround.equals(other.aasAround)) &&
               (customIntervalID == null ? other.customIntervalID == null : customIntervalID.equals(other.customIntervalID));
    }

    public String toString() {
        return "[Contig: " + contig +
               " Position: " + position +
               " Reference: " + reference +
               " Change: " + change +
               " Change Type: " + changeType +
               " Zygosity: " + zygosity +
               " Quality: " + quality +
               " Coverage: " + coverage +
               " Warnings: " + warnings +
               " Gene ID: " + geneID +
               " Gene Name: " + geneName +
               " Bio Type: " + bioType +
               " Transcript ID: " + transcriptID +
               " Exon ID: " + exonID +
               " Exon Rank: " + exonRank +
               " Non-Coding Gene: " + isNonCodingGene +
               " Effect: " + effect +
               " Effect Extra Information: " + effectExtraInformation +
               " Old/New AA: " + oldAndNewAA +
               " Old/New Codon: " + oldAndNewCodon +
               " Codon Num: " + codonNum +
               " CDS Size: " + cdsSize +
               " Codons Around: " + codonsAround +
               " AAs Around: " + aasAround +
               " Custom Interval ID: " + customIntervalID +
               "]";
    }
}
