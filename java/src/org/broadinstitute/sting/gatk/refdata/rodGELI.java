package org.broadinstitute.sting.gatk.refdata;

import java.util.Iterator;
import java.io.IOException;
import java.io.File;

import edu.mit.broad.picard.genotype.geli.GeliFileReader;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
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
        return gh.toString();
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
}