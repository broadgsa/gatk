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
     * For point genotypes (those that report themselves as isPointGenotype()), the alleles are observed bases. For indels
     * (i.e. if isIndelGenotype() is true), a no-event allele is represented by an empty string, 
     * while the event allele is represented by '+' (insertion) or '-' (deletion) followed by actual inserted or deleted bases. 
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
     * is similar to getFWDRefBases() except it returns a single char. For non-point variants (i.e. indels), behavior
     * of this method is not specified and can be defined by implementing classes.
     * @return
     */
    char getRef() ;

    /** Returns true if the object represents a point genotype (e.g. reference base vs. a SNP),
     * regardless of what the actual observation is (see isReference() and isSNP()).
     * @return
     */
    boolean isPointGenotype() ;
    
    /** Returns true if the object represents an extended (indel) genotype (e.g. insertion/deletion event vs. no event),
     * regardless of what the actual observation is (see isReference(), isIndel(), isInsertion(), isDeletion()).
     * @return
     */
    boolean isIndelGenotype() ;

    
    /** Is this variant a SNP? Returns true only if this is point genotype AND a non-reference base was observed, false otherwise.
    *
    * @return true or false
    */
   boolean isSNP();

   /** Returns true if all observed alleles are reference ones (only reference base was observed if this is a point genotype,
    * or no indel events were observed if this is an indel genotype). All is"Variant" methods (where Variant=SNP,Insertion, etc) should
    * return false at a site where isReference() is true to ensure consistency. 
    * @return
    */
   boolean isReference();

   
   /** Is this variant an insertion? Returns true if this is indel genotype AND only insertion event(s) was observed,
    * false otherwise. The contract requires isIndel() to return true
    * whenever isInsertion() method returns true.
    *
    * @return true or false
    */
   boolean isInsertion();

   
   /** Is this variant a deletion? Returns true if this is indel genotype AND only deletion event(s) was observed,
    * false otherwise. The contract requires isIndel() to return true
    * whenever isDeletion() method returns true.
    *
    * @return true or false
    */
   boolean isDeletion();

  
   /** Is this variant an insertion or a deletion? The contract requires
    * this to be true if either isInsertion() or isDeletion() returns true. However,
    * this method is also allowed to return true even if neither isInsertion()
    * nor isDeletion() do (e.g. a het genotype with one insertion and one deletion). Always returns
    * false if the object does not represent indel genotype. 
    * @return
    */
   boolean isIndel();
  
   
   /** Returns phred-scaled confidence for the variation event (e.g. MAQ's SNP confidence, or AlleleCaller's best vs. ref).
   *
   * @return
   */
  double getVariantConfidence();

  /** Returns phred-mapped confidence in called genotype (e.g. MAQ's consensus confidence, or AlleleCaller's
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
    
  int length();


    /** Return actual number of observed alleles (chromosomes) in the genotype.
     * @return
     */
//    int getPloidy() throws IllegalStateException;
    
 }
