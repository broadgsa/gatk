package edu.mit.broad.sting.utils;

import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.picard.util.TabbedTextFileParser;

import java.io.File;
import java.io.InputStream;
import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.util.Iterator;
import java.util.HashMap;

/**
 * Example format:
 *   585 chr1 433 433 rs56289060  0  +  - - -/C  genomic  insertion unknown 0  0  unknown  between  1
 *   585 chr1 491 492 rs55998931  0  +  C C C/T  genomic  single   unknown 0 0 unknown exact 1
 *
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class rodDbSNP extends ReferenceOrderedDatum {
    public  String contig;      // Reference sequence chromosome or scaffold
    public  long start, stop;   // Start and stop positions in chrom

    public  String name;        // Reference SNP identifier or Affy SNP name
    public  String strand;      // Which DNA strand contains the observed alleles
    public  String observed;    // The sequences of the observed alleles from rs-fasta files
    public  char[] observedBases;    // The sequences of the observed alleles from rs-fasta files
    public  String molType;     // Sample type from exemplar ss
    public  String varType;     // The class of variant (simple, insertion, deletion, range, etc.)
                                // Can be 'unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion'
    public String validationStatus;    // The validation status of the SNP
                                        // one of set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap')

    public double avHet;        // The average heterozygosity from all observations
    public double avHetSE;      // The Standard Error for the average heterozygosity

    public String func;         // The functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)
                                // set('unknown','coding-synon','intron','cds-reference','near-gene-3','near-gene-5',
                                // 'nonsense','missense','frameshift','untranslated-3','untranslated-5','splice-3','splice-5')
    public String locType;      // How the variant affects the reference sequence
                                // enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion')

    public int weight;          // The quality of the alignment

    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public rodDbSNP() {}

    public String getContig() { return this.contig; }
    public long getStart() { return start; }
    public long getStop() { return stop; }

    // ----------------------------------------------------------------------
    //
    // formatting
    //
    // ----------------------------------------------------------------------
    public String toString() {
        return String.format("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%s\t%s\t%d",
                contig, start, stop, name, strand, observed, molType,
                varType, validationStatus, avHet, avHetSE, func, locType, weight );
    }

    public String toSimpleString() {
        return String.format("%s:%s", name, observed);
    }

    public String repl() {
        return String.format("%d\t%s\t%d\t%d\t%s\t0\t%s\tX\tX\t%s\t%s\t%s\t%s\t%f\t%f\t%s\t%s\t%d",
                585, contig, start-1, stop-1, name, strand, observed, molType,
                varType, validationStatus, avHet, avHetSE, func, locType, weight );
    }

    public void parseLine(final String[] parts) {
        try {
            contig = parts[1];
            start = Long.parseLong(parts[2]) + 1; // The final is 0 based

            stop = Long.parseLong(parts[3]) + 1;  // The final is 0 based
            name = parts[4];
            strand = parts[6];
            observed = parts[9];
            molType = parts[10];
            varType = parts[11];
            validationStatus = parts[12];
            avHet = Double.parseDouble(parts[13]);
            avHetSE = Double.parseDouble(parts[14]);
            func = parts[15];
            locType = parts[16];
            weight = Integer.parseInt(parts[17]);

            // Cut up the observed bases string into an array of individual bases
            String[] bases = observed.split("/");
            observedBases = new char[bases.length];
            for ( String elt : bases ) {
                observedBases[0] = (char)elt.getBytes()[0];
                //System.out.printf("  Bases %s %d %c%n", elt, elt.getBytes()[0], (char)elt.getBytes()[0]);
            }
            //System.out.printf("  => Observed bases are %s%n", Utils.join(" B ", bases));
        } catch ( RuntimeException e ) {
            System.out.printf("  Exception caught during parsing GFFLine %s%n", Utils.join(" <=> ", parts));
            throw e;
        }
    }
}