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
public interface Genotype extends Comparable<ReferenceOrderedDatum> {
    // ----------------------------------------------------------------------
    //
    // manipulating the SNP information
    //
    // ----------------------------------------------------------------------
    /** Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     *  is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    GenomeLoc getLocation();

    /** Returns actual observed alleles on the fwd strand. Allowed to return more than two alleles (@see #getPloidy()). 
     * @return list of alleles
     */
    List<String> getFWDAlleles() ;

    
    /** Returns bases on the fwd strand in the true reference allele (regardless of whether the reference allele 
     * was observed at this site in this particular genotype!) 
     * as a String. String can be empty (if alt allele is insertion into
     * the reference), can contain a single character (as in SNP or one-base deletion), or multiple characters
     * (if alt is a longer deletion). For the actually observed alleles, see #getFWDAlleles()
     *
     * @return reference allele, forward strand
     */
    String getFWDRefBases();
    
    /** Returns reference base (regardless of whether the reference allele was actually observed. This method 
     * is similar to getFWDRefBases() except it returns a single char. For non-single letter (non-SNP) variants, behavior
     * of this method is not specified and can be defined by implementing classes.
     * @return
     */
    char getRef() ;

    /** Is this variant a SNP?
    *
    * @return true or false
    */
   boolean isSNP();

   /** Returns true if all observed alleles are reference ones. All is<Variant> methods (where Variant=SNP,Insertion, etc) should
    * return false at a site where isReference() is true to ensure consistency. 
    * @return
    */
   boolean isReference();

   
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
  
   
   /** Returns phred-scaled confidence for the variation event (e.g. MAQ's SNP confidence, or AlleleCaller's best vs. ref).
   *
   * @return
   */
  double getVariantConfidence();

  /** Returns phred-scaled confidence in called genotype (e.g. MAQ's consensus confidence, or AlleleCaller's
   * best vs next-best.
   * @return
   */
  double getConsensusConfidence();


  /** Returns true if the site has at most two known or observed alleles (ref and non-ref), or false if there are > 2 allelic variants known or observed. When
   * the implementing class is a genotype, alleles should be always counted including the reference allele whether it was observed in the particular
   * individual or not: i.e. if the reference is 'C', then both 'CA' and 'AA' genotypes must be reported as bi-allelic, while 'AT' is <i>not</i> bi-allelic (since there are
   * two allelic variants, 'A' and 'T' <i>in addition</i> to the (known) reference variant 'C'). 
   * @return
   */
  boolean isBiallelic() ;
  
  /** Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
   *  
   * @return
   */
  boolean isHom() ;

  /** Returns true if observed alleles differ (regardless of whether they are ref or alt)
   *  
   * @return
   */
  boolean isHet() ;

  /** Returns reference (major) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * @return reference base on the forward strand
     */
    //char getRefSnpFWD() throws IllegalStateException;

    /** Returns bases in the alternative allele as a String. String can be empty (as in deletion from
     * the reference), can contain a single character (as in SNP or one-base insertion), or multiple characters
     * (for longer indels).
     *
     * @return alternative allele, forward strand
     */
    //String getAltBasesFWD();

    /** Returns alternative (minor) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * @return alternative allele base on the forward starnd
     */
  //  char getAltSnpFWD() throws IllegalStateException;
    



    /** Return actual number of observed alleles (chromosomes) in the genotype.
     * @return
     */
//    int getPloidy() throws IllegalStateException;
    
 }
