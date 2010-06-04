package org.broadinstitute.sting.gatk.refdata.features.vcf4;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.util.LineReader;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFReaderUtils;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;

import java.io.IOException;
import java.util.*;

/**
 * a feature codec for the VCF 4 specification.  Our aim is to read in the records and convert to VariantContext as
 * quickly as possible, relying on VariantContext to do the validation of any contradictory (or malformed) record parameters.
 */
public class VCF4Codec implements FeatureCodec {

    // we have to store the list of strings that make up the header until they're needed
    private List<String> headerStrings = new ArrayList<String>();
    private VCFHeader header = null;

    public VCF4Codec() {
        throw new StingException("DON'T USE THIS");
    }

    /**
     * this a super hack like method to parse out what we need from a variant context
     * @param line the line to parse
     * @return
     */
    @Override
    public Feature decode(String line) {

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new IllegalStateException("VCF Header cannot be null");

        // split the line on whitespace (Jim's parser will be faster, but it's broken right now)
        String[] result = line.split("\\s+");

        // check to make sure the split resulted in the correct number of fields (8 + (1 + genotytpe counts if it has genotypes)
        if (result.length != header.getColumnCount()) throw new IllegalArgumentException("we expected " + header.getColumnCount() + " columns and we got " + result.length+ " for line " + line);

        // our genotype names
        Iterator<String> iter = header.getGenotypeSamples().iterator();

        // out genotype map, sample name to genotype
        Map<String, Genotype> genotypes = new LinkedHashMap<String,Genotype>();

        // our allele list, add the reference and the alts
        List<Allele> alleles = new ArrayList<Allele>();
        String[] alts = result[4].split(",");
        for (String alt : alts)
            alleles.add(new Allele(alt,false));
        alleles.add(new Allele(result[3],true));

        // parse out each of the genotypes
        for (int genotypeIndex = 9; genotypeIndex < header.getColumnCount(); genotypeIndex++) {
            if (!iter.hasNext()) throw new StingException("Wrong number of samples!");
            String sample = iter.next();
            genotypes.put(sample,createGenotypeFromString(sample,result[genotypeIndex],result[8].split(":"),alts,result[3]));
        }

        // make a new set of all the filters
        Set<String> filters = new TreeSet<String>();
        filters.addAll(Arrays.asList(result[5].split(",")));

        // create, validate, and return the record
        VCF4Record rec = new VCF4Record(result[2],
                GenomeLocParser.createGenomeLoc(result[0],Long.valueOf(result[1])),
                        Collections.unmodifiableCollection(alleles),
                        genotypes,
                        Double.valueOf(result[5]),
                        filters,
                        new HashMap<String,Object>());
        rec.validate();
        return rec;
    }


    /**
     * this method is a big hack, since I haven't gotten to updating the VCF header for the 4.0 updates
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    @Override
    public int readHeader(LineReader reader) {
        String line = "";
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("##")) {
                    headerStrings.add(line);
                }
                else if (line.startsWith("#")) {
                    headerStrings.add(line);
                    String[] genotypes = line.split("\\s+");
                    Set<String> genotypeSampleNames = new TreeSet<String>();
                    for (int x = 8; x < genotypes.length; x++)
                        genotypeSampleNames.add(genotypes[x]);
                    // this should be the next line -> header = VCFReaderUtils.createHeader(headerStrings);
                    header = new VCFHeader(new HashSet<VCFHeaderLine>(),genotypeSampleNames);
                    return headerStrings.size();
                }
                else {
                    throw new CodecLineParsingException("We never saw the required header line (starting with one #) for the input VCF file");
                }

            }
        } catch (IOException e) {
            throw new RuntimeException("IO Exception ", e);
        }
        throw new CodecLineParsingException("We never saw the required header line (starting with one #) for the input VCF file");
    }


    /**
     * create the genotype object from the VCF record
     * @param name the name of the sample
     * @param VCFEntry the entry; the text containing all the fields corresponding to the format fields
     * @param formatStrings the format string
     * @param altAlleleList our alt alleles
     * @param reference the reference base(s)
     * @return a Genotype object
     */
    private Genotype createGenotypeFromString(String name, String VCFEntry, String[] formatStrings, String[] altAlleleList, String reference) {
        // split the text entry into parts
        String genotypeSplit[] = VCFEntry.split(":");

        Set<Allele> aList = new TreeSet<Allele>();
        Map<String, Object> attributes = new LinkedHashMap<String,Object>();

        // for each entry in the vcf field (we drive by this so that dropped fields aren't processed
        for (int index = 0; index < genotypeSplit.length; index++) {
            if (formatStrings[index].toUpperCase().equals("GT")) {
                String[] genotypes = genotypeSplit[index].split("[\\\\|\\/]+");
                for (String g : genotypes) {
                    int altIndex = Integer.valueOf(g);
                    if (altIndex == 0)
                        aList.add(new Allele(reference,true));
                    else
                        aList.add(new Allele(altAlleleList[altIndex-1]));
                }
            } else {
                attributes.put(formatStrings[index],genotypeSplit[index]);
            }
        }
        return new Genotype(name,new ArrayList(aList),0.0,new HashSet<String>(),attributes,false);
    }

    /**
     *
     * @return the type of record
     */
    @Override
    public Class getFeatureType() {
        return VCF4Record.class;
    }
}
