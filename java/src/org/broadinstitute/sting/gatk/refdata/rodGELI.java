package org.broadinstitute.sting.gatk.refdata;

import java.util.*;
import java.io.IOException;
import java.io.File;

import edu.mit.broad.picard.genotype.geli.GeliFileReader;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import edu.mit.broad.picard.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import net.sf.samtools.util.CloseableIterator;


/**
 *  This class wraps Picard Geli CHiP data and presents it as a ROD.
 */

public class rodGELI extends BasicReferenceOrderedDatum {
	// ----------------------------------------------------------------------
	//
	// Constructors
	//
	// ----------------------------------------------------------------------

	private GenotypeLikelihoods gh = null;

	public rodGELI(final String name, GenotypeLikelihoods gh) {
		super(name);
        this.gh = gh;
    }

	@Override
	public GenomeLoc getLocation() {
		return GenomeLocParser.createGenomeLoc(gh.getSequenceIndex(), gh.getPosition());
	}

	/** Required by ReferenceOrderedDatum interface. This implementation provides its own iterator,
	 * so this method does nothing at all (always returns false).
	 *
	 */
	@Override
	public boolean parseLine(Object header, String[] parts) throws IOException {
		return false;
	}

	@Override
	public String toString() {
        final StringBuilder builder = new StringBuilder();
        builder.append(gh.getSequenceName() + "\t");
        builder.append(gh.getPosition() + "\t");
        builder.append(Character.toString((char)(gh.getReferenceBase() & 0xff)) + "\t");
        builder.append(gh.getNumReads() + "\t");
        builder.append(gh.getMaxMappingQuality() + "\t");
        builder.append(gh.getBestGenotype().name() + "\t");
        builder.append(gh.getBestToReferenceLod() + "\t");
        builder.append(gh.getBestToSecondBestLod() + "\t");
        builder.append("\t"); // no dbSNP info in GenotypeLikelihoods class
        for (final DiploidGenotype genotype : DiploidGenotype.values())
            builder.append(gh.getLikelihood(genotype) + "\t");
        builder.append("\n");
        return builder.toString();
	}

	private static class rodGELIIterator implements Iterator<rodGELI> {

		private String rodName = null;
		private GeliFileReader parser = null;
        private CloseableIterator<GenotypeLikelihoods> iterator = null;

		rodGELIIterator(String name, File f) {
            rodName = name;
			parser = new GeliFileReader(f);
			iterator = parser.iterator();
		}

		public boolean hasNext() {
			return iterator.hasNext();
		}

		public rodGELI next() {
            return new rodGELI(rodName, iterator.next());
		}

		public void remove() {
			throw new UnsupportedOperationException("'remove' operation is not supported for GELIs");
		}

	}

	public static Iterator<rodGELI> createIterator(String name, File file) {
		return new rodGELI.rodGELIIterator(name,file);
	}

    /****
     * We no longer want to use this class as truth data (sigh)
     *
    public List<String> getFWDAlleles() {
        return new ArrayList<String>();
    }

    public String getFWDRefBases() {
        return "";
    }

    public char getRef() {
        return (char)(gh.getReferenceBase() & 0xff);
    }

    public boolean isPointGenotype() { return true; }
    public boolean isIndelGenotype() { return false; }
    public boolean isSNP() { return true; }
    public boolean isReference() { return gh.isHomozygousReference(); }
    public boolean isInsertion() { return false; }
    public boolean isDeletion() { return false; }
    public boolean isIndel() { return false; }
    public boolean isBiallelic() { return false; }
    public boolean isHom() { return gh.isHomozyous(); }
    public boolean isHet() { return !isHom(); }

    public double getVariantConfidence() {
        return gh.getBestToReferenceLod();
    }

    public double getConsensusConfidence() {
        return gh.getBestToSecondBestLod();
    }
    *
    */

    public static void main(String argv[]) {
        String testFile = "NA12878.geli";

        Iterator<rodGELI> it = createIterator("test-geli", new File(testFile));

        net.sf.picard.reference.ReferenceSequenceFileWalker reference = new net.sf.picard.reference.ReferenceSequenceFileWalker(new File(  "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));

        if ( reference.getSequenceDictionary() == null ) {
            System.out.println("No reference sequence dictionary found. Abort.");
            System.exit(1);
        }

        GenomeLocParser.setupRefContigOrdering(reference.getSequenceDictionary());

        int counter = 0;

        while ( it.hasNext() && counter < 500 ) {
            rodGELI rg = it.next();
            System.out.println(rg.toString());
            counter++;
        }
    }
}