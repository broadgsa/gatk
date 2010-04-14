package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.AnnotatorROD;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.TabularROD;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;




/**
 * This plugin for {@link VariantAnnotatorEngine} serves as the core
 * of the {@link GenomicAnnotator}. It finds all records in the -B input files
 * that match the given variant's position and, optionally, it's reference and alternate alleles.
 * Whether or not matching is done by reference and alternate alleles for a particular input file
 * based solely on whether the given -B input has columns named "haplotypeReference" and
 * "haplotypeAlternate".
 */
public class GenomicAnnotation implements InfoFieldAnnotation {

    private static final String HAPLOTYPE_REFERENCE_COLUMN = AnnotatorROD.HAPLOTYPE_REFERENCE_COLUMN;
    private static final String HAPLOTYPE_ALTERNATE_COLUMN = AnnotatorROD.HAPLOTYPE_ALTERNATE_COLUMN;
    private static final String HAPLOTYPE_STRAND_COLUMN = AnnotatorROD.HAPLOTYPE_STRAND_COLUMN;


    /**
     * For each ROD (aka. record) which overlaps the current locus, generates a
     * set of annotations of the form:
     *
     * thisRodName.fieldName1=fieldValue, thisRodName.fieldName1=fieldValue
     *
     * Examples of generated annotations are: dbSNP.avHet=0.7, dbSNP.ref_allele=A, etc.
     *
     * @return The following is an explanation of this method's return value:
     * The annotations (aka. columnName=fieldValue pairs) from a matching record in a particular file are stored in a Map<String, String>.
     * Since a single input file can have multiple records that overlap the current
     * locus (eg. dbSNP can have multiple entries for the same location), a different
     * Map<String, String> is created for each matching record in a particular file.
     * The list of matching records for each file is then represented as a List<Map<String, String>>
     *
     * The return value of this method is a Map<String, Object> of the form:
     *     rodName1 -> List<Map<String, String>>
     *     rodName2 -> List<Map<String, String>>
     *     rodName3 -> List<Map<String, String>>
     *     ...
     * Where the rodNames are the -B binding names for each file that were specified on the command line.
     *
     * NOTE: The lists (List<Map<String, String>>) are guaranteed to have size > 0.
     * The reason is that a rodName -> List<Map<String, String>> entry will only
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
            final ReferenceOrderedDatum rod = (ReferenceOrderedDatum) gatkFeature.getUnderlyingObject();
            final String name = rod.getName();
            if( name.equals("variant") || name.equals("interval") ) {
                continue;
            }

            if( ! (rod instanceof TabularROD) ) {
                continue; //GenericAnnotation only works with TabularRODs because it needs to be able to select individual columns.
            }

            TabularROD tabularRod = (TabularROD) rod;
            final Map<String, String> annotationsForRecord = convertRecordToAnnotations( tabularRod );

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
                        throw new StingException("Record (" + rod + ") in " + name + " contains " + alternateAlleles.size() + " alternate alleles. GenomicAnnotion currently only supports annotating 1 alternate allele.");
                    }

                    boolean positiveStrand = true; //if HAPLOTYPE_STRAND_COLUMN isn't specified, assume positive strand.
                    String hapStrandValue = annotationsForRecord.get( generateInfoFieldKey(name, HAPLOTYPE_STRAND_COLUMN) );
                    if(hapStrandValue != null ) {
                        hapStrandValue = hapStrandValue.trim().toLowerCase();
                        if(hapStrandValue.equals("-") || hapStrandValue.equals("r")) {
                            positiveStrand = false;
                        } else if(!hapStrandValue.equals("+") && !hapStrandValue.equals("f")) {
                            throw new StingException("Record (" + rod + ") in " + name + " has an invalid value for " + HAPLOTYPE_STRAND_COLUMN + ". This value is: \"" + hapStrandValue + "\"");
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
                        //TODO check the intersection of vcRef and hapRefValue
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
     * Converts the ROD to a set of key-value pairs of the form:
     *   thisRodName.fieldName1=fieldValue, thisRodName.fieldName1=fieldValue
     *   (eg. dbSNP.avHet=0.7, dbSNP.ref_allele=A)
     *
     * @param rod A TabularROD corresponding to one record in one input file.
     *
     * @return The map of column-name -> value pairs.
     */
    private Map<String, String> convertRecordToAnnotations( final TabularROD rod ) {
        final String rodName = rod.getName(); //aka the rod binding

        final Map<String, String> result = new HashMap<String, String>();
        for(final Entry<String, String> entry : rod.entrySet()) {
            final String value = entry.getValue();
            if(!value.trim().isEmpty()) {
                result.put( generateInfoFieldKey(rodName, entry.getKey()), entry.getValue());
            }
        }

        return result;
    }

    public static String generateInfoFieldKey(String rodBindingName, String columnName ) {
        return rodBindingName + "." + columnName;
    }


    /**
     * Parses the columns arg and returns a Map of columns hashed by their binding name.
     * For example:
     *   The command line:
     *      -s dbSnp.valid,dbsnp.avHet -s refGene.txStart,refGene.txEnd
     *
     *   will be passed to this method as:
     *       ["dbSnp.valid,dbsnp.avHet", "refGene.txStart,refGene.txEnd"]
     *
     *   resulting in a return value of:
     *      {
     *       "dbSnp" -> "dbSnp.valid" ,
     *       "dbSnp" -> "dbsnp.avHet" ,
     *       "refGene" -> "refGene.txStart",
     *       "refGene" -> "refGene.txEnd"
     *      }
     *
     * @param columnsArg The -s command line arg value.
     *
     * @return Map representing a parsed version of this arg - see above.
     */
    public static Map<String, Set<String>> parseColumnsArg(String[] columnsArg) {
        Map<String, Set<String>> result = new HashMap<String, Set<String>>();

        for(String s : columnsArg) {
            for(String columnSpecifier : s.split(",") ) {
                String[] rodNameColumnName = columnSpecifier.split("\\.");
                if(rodNameColumnName.length != 2) {
                    throw new IllegalArgumentException("The following column specifier in the -s arg is invalid: [" + columnSpecifier + "]. It must be of the form 'bindingName.columnName'.");
                }
                String rodName = rodNameColumnName[0];
                //String columnName = rodNameColumnName[1];

                Set<String> requestedColumns = result.get(rodName);
                if(requestedColumns == null) {
                    requestedColumns = new HashSet<String>();
                    result.put(rodName, requestedColumns);
                }
                requestedColumns.add(columnSpecifier);
            }
        }

        return result;
    }

    public VCFInfoHeaderLine getDescription() {
        return new VCFInfoHeaderLine("GenericAnnotation", 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "For each variant in the 'variants' ROD, finds all entries in the other -B files that overlap the variant's position. ");
    }

    public String getKeyName() {
        return "GenericAnnotation";
    }

}
