package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Apr 1, 2009
 * Time: 11:02:03 AM
 * To change this template use File | Settings | File Templates.
 */
public interface AllelicVariant extends ReferenceOrderedDatum {
    // ----------------------------------------------------------------------
    //
    // manipulating the SNP information
    //
    // ----------------------------------------------------------------------
    /** Location of this variant on the reference (on the forward strand).
     *
     * @return
     */
    GenomeLoc getLocation();

    /** Returns bases in the reference allele as a String. String can be empty (as in insertion into
     * the reference), can contain a single character (as in SNP or one-base deletion), or multiple characters
     * (for longer indels).
     *
     * @return reference allele, forward strand
     */
    String getRefBasesFWD();

    /** Returns reference (major) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * @return reference base on the forward strand
     */
    char getRefSnpFWD() throws IllegalStateException;

    /** Returns bases in the alternative allele as a String. String can be empty (as in deletion from
     * the reference), can contain a single character (as in SNP or one-base insertion), or multiple characters
     * (for longer indels).
     *
     * @return alternative allele, forward strand
     */
    String getAltBasesFWD();

    /** Returns alternative (minor) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * @return alternative allele base on the forward starnd
     */
    char getAltSnpFWD() throws IllegalStateException;
    
    /** Returns true if all observed alleles are reference alleles. All is<Variant> methods (where Variant=SNP,Insertion, etc) should
     * return false at such site to ensure consistency. This method is included for use with genotyping calls (isGenotype()==true), it makes
     * no sense for, e.g. dbSNP and should return false for the latter.
     * @return
     */
    boolean isReference();

    /** Is this variant a SNP?
     *
     * @return true or false
     */
    boolean isSNP();

    /** Is this variant an insertion? The contract requires isIndel() to return true
     * if this method returns true.
     *
     * @return true or false
     */
    boolean isInsertion();

    /** Is this variant a deletion? The contract requires isIndel() to return true
     * if isDeletion() returns true.
     *
     * @return true or false
     */
    boolean isDeletion();

    /** Is this variant an insertion or a deletion? The contract requires
     * this to be true if either isInsertion() or isDeletion() returns true. However,
     * this method is currently allowed to return true even if neither of isInsertion()
     * and isDeletion() does.
     * @return
     */
    boolean isIndel();

    /** Returns minor allele frequency.
     *
     * @return
     */
    double getMAF() ;

    /** Returns heterozygosity, a more accessible general feature of a variant.
     *
     * @return
     */
    double getHeterozygosity() ;
    
    /** Is this variant an actual genotype (such as individual call from sequencing, HapMap chip etc), or
     * population allelic variant (call from pooled sequencing, dbSNP site etc). Only if variant is a genotype, there
     * is a meaningful question of, e.g., whether it is a het, or homogeneous non-ref.
     *
      * @return true if this variant is an actual genotype.
     */
    boolean isGenotype();

    /** Returns phred-mapped confidence in variation event (e.g. MAQ's SNP confidence, or AlleleCaller's best vs. ref).
     *
     * @return
     */
    double getVariationConfidence();

    /** Returns phred-mapped confidence in called genotype (e.g. MAQ's consensus confidence, or AlleleCaller's
     * best vs next-best.
     * @return
     */
    double getConsensusConfidence();

    /** Returns actual observed genotype. Allowed to return more than two alleles (@see #getPloidy()). If this variant
     * is not a genotype, should throw an IllegalStateException.
     * @return
     */
    List<String> getGenotype() throws IllegalStateException;

    /** Return actual number of observed alleles (chromosomes) in the genotype. If this variant is not a genotype,
     * should throw IllegalStateException.
     * @return
     */
    int getPloidy() throws IllegalStateException;
    
    /** Returns true if the site has at most two known or observed alleles (ref and non-ref), or false if there are > 2 allelic variants known or observed. When
     * the implementing class is a genotype, alleles should be always counted including the reference allele whether it was observed in the particular
     * individual or not: i.e. if the reference is 'C', then both 'CA' and 'AA' genotypes must be reported as bi-allelic, while 'AT' is <i>not</i> bi-allelic (since there are
     * two allelic variants, 'A' and 'T' <i>in addition</i> to the (known) reference variant 'C'). 
     * @return
     */
    boolean isBiallelic();

    /** returns the length of the variant.  For SNPs this is just 1.
     */
    int length();
}
