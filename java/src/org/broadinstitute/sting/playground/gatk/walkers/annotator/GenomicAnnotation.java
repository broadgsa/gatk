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
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.TabularROD;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;




/**
 * TODO comment
 */
public class GenomicAnnotation implements InfoFieldAnnotation {

    private static final String HAP_REF = "hap_ref";
    private static final String HAP_ALT = "hap_alt";


    /**
     * For each ROD (aka. record) which overlaps the current locus, generates a
     * set of annotations of the form:
     *
     * thisRodName.fieldName1=fieldValue, thisRodName.fieldName1=fieldValue (eg. dbSNP.avHet=0.7, dbSNP.ref_allele=A),
     *
     * These annotations are stored in a Map<String, String>.
     *
     * Since a single input file can have multiple records that overlap the current
     * locus (eg. dbSNP can have multiple entries for the same location), a different
     * Map<String, String> is created for each of these, resulting in a List<Map<String, String>>
     * for each input file.
     *
     * The return value of this method is a Map of the form:
     *     rodName1 -> List<Map<String, String>>
     *     rodName2 -> List<Map<String, String>>
     *     rodName3 -> List<Map<String, String>>
     *     ...
     *
     * The List values are guaranteed to have size > 0, and in most cases will have size == 1.
     */
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
            final ReferenceContext ref,
            final Map<String, StratifiedAlignmentContext> stratifiedContexts,
            final VariantContext vc) {

        final Map<String, Object> annotations = new HashMap<String, Object>();
        for(final GATKFeature gatkFeature : tracker.getAllRods())
        {
            final ReferenceOrderedDatum rod = (ReferenceOrderedDatum) gatkFeature.getUnderlyingObject();
            final String name = rod.getName();
            if( name.equals("variant") || name.equals("interval") ) {
                continue;
            }

            if( ! (rod instanceof TabularROD) ) {
                continue; //GenericAnnotation only works with TabularRODs so that it can select individual columns
            }

            TabularROD tabularRod = (TabularROD) rod;
            final Map<String, String> annotationsForRecord = convertRecordToAnnotations( tabularRod );

            //filter by matching the vc's reference allele and alternate allele to the record's HAP_REF and HAP_ALT columns.
            final String hapAltValue = annotationsForRecord.get( generateInfoFieldKey(name, HAP_ALT) );
            if(hapAltValue != null && !hapAltValue.equals("*")) {
                Set<Allele> alternateAlleles = vc.getAlternateAlleles();
                if(alternateAlleles.isEmpty()) {
                    //TODO If this site is monomorphic in the VC, and the current record specifies a particular alternate allele, skip this record. Right?
                    continue;
                } else if(alternateAlleles.size() > 1) {
                    throw new StingException("Record (" + rod + ") in " + name + " contains " + alternateAlleles.size() + " alternate alleles. GenomicAnnotion currently only supports annotating 1 alternate allele.");
                }

                //TODO handle strand somehow?

                //decide whether to skip this record.
                Allele vcAlt = alternateAlleles.iterator().next();
                if(!vcAlt.basesMatch(hapAltValue)) {
                    //TODO check the intersection of vcAlt and hapAltValue
                    continue; //skip record
                }
            }

            final String hapRefValue = annotationsForRecord.get( generateInfoFieldKey(name, HAP_REF) );
            if(hapRefValue != null && !hapRefValue.equals("*")) {
                //match against hapRef
                Allele vcRef = vc.getReference();
                if(!vcRef.basesMatch(hapRefValue)) {
                    //TODO check the intersection of vcRef and hapRefValue
                    continue; //skip record
                }
            }

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
            result.put( generateInfoFieldKey(rodName, entry.getKey()), entry.getValue());
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
