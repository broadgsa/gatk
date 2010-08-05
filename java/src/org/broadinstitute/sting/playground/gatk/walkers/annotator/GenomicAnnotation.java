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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.util.*;
import java.util.Map.Entry;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;

/**
 * This plugin for {@link VariantAnnotatorEngine} serves as the core
 * of the {@link GenomicAnnotator}. It finds all records in the -B input files
 * that match the given variant's position and, optionally, the variant's reference and alternate alleles.
 *
 * For details, see:  http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
public class GenomicAnnotation implements InfoFieldAnnotation {

    public static final String CHRPOS_COLUMN = "chrpos";
    public static final String HAPLOTYPE_REFERENCE_COLUMN = "haplotypeReference";
    public static final String HAPLOTYPE_ALTERNATE_COLUMN = "haplotypeAlternate";
    public static final String HAPLOTYPE_STRAND_COLUMN = "haplotypeStrand";

    public static final String NUM_MATCHES_SPECIAL_INFO_FIELD = "numMatchingRecords";

    /** Characters that aren't allowed within VCF info field key-value pairs */
    public static final char[] ILLEGAL_INFO_FIELD_VALUES            = {  ' ', '=', ';' };
    /** Replacement for each character in ILLEGAL_INFO_FIELD_VALUES */
    public static final char[] ILLEGAL_INFO_FIELD_VALUE_SUBSTITUTES = {  '_', '-',  '!'  };

    /**
     * For each -B input file, for each record which overlaps the current locus, generates a
     * set of annotations of the form:
     *
     * bindingName.columnName1=columnValue, bindingName.columnName2=columnValue2, etc.
     *
     * For example: dbSNP.avHet=0.7, dbSNP.ref_allele=A, etc.
     *
     * @return The following is an explanation of this method's return value:
     *
     * The annotations from a matching in a particular file are stored in a Map<String, String>
     * where the key is bindingName.columnName and the value is the columnValue.
     * Since a single input file can have multiple records that overlap the current
     * locus (eg. dbSNP can have multiple entries for the same genomic position), a different
     * Map<String, String> is created for each matching record in a particular file.
     * The set of matching records for each file is then represented as a List<Map<String, String>>
     *
     * The return value of this method is a Map<String, Object> of the form:
     *     rodName1 -> List<Map<String, String>>
     *     rodName2 -> List<Map<String, String>>
     *     rodName3 -> List<Map<String, String>>
     *     ...
     * Where the rodNames are the -B binding names for each file that were specified on the command line (eg. -B bindingName,AnnotatorInputTable,/path/to/file).
     *
     * NOTE: The lists (List<Map<String, String>>) are guaranteed to have size > 0
     * because a  rodName -> List<Map<String, String>>  entry will only
     * be created in Map<String, Object> if the List has at least one element.
     */
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
            final ReferenceContext ref,
            final Map<String, StratifiedAlignmentContext> stratifiedContexts,
            final VariantContext vc) {

        //iterate over each record that overlaps the current locus, and, if it passes certain filters,
        //add its values to the list of annotations for this locus.
        final Map<String, Object> annotations = new HashMap<String, Object>();
        for(final GATKFeature gatkFeature : tracker.getAllRods())
        {
            final String name = gatkFeature.getName();
            if( name.equals("variant") || name.equals("interval") ) {
                continue;
            }

            if( ! (gatkFeature.getUnderlyingObject() instanceof AnnotatorInputTableFeature) ) {
                continue; //GenericAnnotation only works with TabularRODs because it needs to be able to select individual columns.
            }

            final Map<String, String> annotationsForRecord = convertRecordToAnnotations( gatkFeature.getName(), ((AnnotatorInputTableFeature) gatkFeature.getUnderlyingObject()).getColumnValues());

            //If this record contains the HAPLOTYPE_REFERENCE_COLUMN and/or HAPLOTYPE_ALTERNATE_COLUMN, check whether the
            //alleles specified match the the variant's reference allele and alternate allele.
            //If they don't match, this record will be skipped, and its values will not be used for annotations.
            //
            //If one of these columns doesn't exist in the current rod, or if its value is * (star), then this is treated as an automatic match.
            //Otherwise, the HAPLOTYPE_REFERENCE_COLUMN is only considered to be matching the variant's reference if the string values of the two
            //are exactly equal (case-insensitive).

            //The HAPLOTYPE_REFERENCE_COLUMN is matches the variant's reference allele based on a case-insensitive string comparison.
            //The HAPLOTYPE_ALTERNATE_COLUMN is can optionally list more than allele separated by one of these chars: ,\/:|
            //The matches if any of the
            String hapAltValue = annotationsForRecord.get( generateInfoFieldKey(name, HAPLOTYPE_ALTERNATE_COLUMN) );
            if(hapAltValue != null)
            {
                if(!hapAltValue.equals("*"))
                {
                    Set<Allele> alternateAlleles = vc.getAlternateAlleles();
                    //if(alternateAlleles.isEmpty()) {
                        //handle a site that has been called monomorphic reference
                        //alternateAlleles.add(vc.getReference());
                        //continue;            //TODO If this site is monomorphic in the VC, and the current record specifies a particular alternate allele, skip this record. Right?
                    //} else
                    if(alternateAlleles.size() > 1) {
                        throw new StingException("Record [" + vc + "] contains " + alternateAlleles.size() + " alternate alleles. GenomicAnnotion currently only supports annotating 1 alternate allele.");
                    }

                    boolean positiveStrand = true; //if HAPLOTYPE_STRAND_COLUMN isn't specified, assume positive strand.
                    String hapStrandValue = annotationsForRecord.get( generateInfoFieldKey(name, HAPLOTYPE_STRAND_COLUMN) );
                    if(hapStrandValue != null ) {
                        hapStrandValue = hapStrandValue.trim().toLowerCase();
                        if(hapStrandValue.equals("-") || hapStrandValue.equals("r")) {
                            positiveStrand = false;
                        } else if(!hapStrandValue.equals("+") && !hapStrandValue.equals("f")) {
                            throw new StingException("Record (" + gatkFeature.getUnderlyingObject() + ") in " + name + " has an invalid value for " + HAPLOTYPE_STRAND_COLUMN + ". This value is: \"" + hapStrandValue + "\"");
                        }
                    }


                    Allele vcAlt;
                    if(alternateAlleles.isEmpty()) {
                        vcAlt = vc.getReference();
                    } else {
                        vcAlt = alternateAlleles.iterator().next();
                    }

                    boolean matchFound = false;
                    for(String hapAlt : hapAltValue.split("[,\\\\/:|]")) {
                        if(!positiveStrand) {
                            hapAlt = BaseUtils.simpleReverseComplement(hapAlt);
                        }

                        if(!hapAlt.isEmpty() && vcAlt.basesMatch(hapAlt)) {
                            matchFound = true;
                            break;
                        }
                    }
                    if(!matchFound) {
                        continue; //skip record - none of its alternate alleles match the variant's alternate allele
                    }
                }
            }

            String hapRefValue = annotationsForRecord.get( generateInfoFieldKey(name, HAPLOTYPE_REFERENCE_COLUMN) );
            if(hapRefValue != null)
            {
                hapRefValue = hapRefValue.trim();
                if(!hapRefValue.equals("*"))
                {
                    //match against hapolotypeReference.
                    Allele vcRef = vc.getReference();
                    if(!vcRef.basesMatch(hapRefValue)) {
                        continue; //skip record
                    }
                }
            }

            //filters passed, so add this record.
            List<Map<String, String>> listOfMatchingRecords = (List<Map<String, String>>) annotations.get( name );
            if(listOfMatchingRecords == null) {
                listOfMatchingRecords = new LinkedList<Map<String,String>>();
                listOfMatchingRecords.add( annotationsForRecord );
                annotations.put(name, listOfMatchingRecords);
            } else {
                listOfMatchingRecords.add( annotationsForRecord );
            }
        }

        return annotations;
    }




    /**
     * Converts the given record to a set of key-value pairs of the form:
     *   bindingName.columnName1=column1Value, bindingName.columnName2=column2Value
     *   (eg. dbSNP.avHet=0.7, dbSNP.ref_allele=A)
     *
     * @param record AnnotatorInputTableFeature corresponding to one record in one -B input file.
     * @param bindingName The binding name of the given AnnotatorInputTableFeature.
     * @return The map of columnName -> columnValue pairs.
     */
    public static Map<String, String> convertRecordToAnnotations( String bindingName, Map<String, String> record) {
        final Map<String, String> result = new HashMap<String, String>();

        for(final Entry<String, String> entry : record.entrySet()) {
            final String value = entry.getValue();
            if(!value.trim().isEmpty()) {
                result.put( generateInfoFieldKey(bindingName, entry.getKey()), scrubInfoFieldValue(entry.getValue()));
            }
        }

        return result;
    }

    /**
     * Combines the 2 values into a full key.
     * @param rodBindingName
     * @param columnName
     * @return
     */
    public static String generateInfoFieldKey(String rodBindingName, String columnName ) {
        return rodBindingName + '.' + columnName;
    }



    /**
     * Replaces any characters that are not allowed in the info field of a VCF file.
     *
     * @param value
     * @return the value with any illegal characters replaced by legal ones.
     */
    private static String scrubInfoFieldValue(String value) {
        for(int i = 0; i < GenomicAnnotation.ILLEGAL_INFO_FIELD_VALUES.length; i++) {
            value = value.replace(GenomicAnnotation.ILLEGAL_INFO_FIELD_VALUES[i], GenomicAnnotation.ILLEGAL_INFO_FIELD_VALUE_SUBSTITUTES[i]);
        }

        return value;
    }



    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("GenericAnnotation", 1, VCFHeaderLineType.Integer, "For each variant in the 'variants' ROD, finds all entries in the other -B files that overlap the variant's position. "));
    }

    public List<String> getKeyNames() {
        return Arrays.asList("GenericAnnotation");
    }

}
