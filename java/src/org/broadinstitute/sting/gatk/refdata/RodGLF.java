package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.genotype.glf.GLFReader;
import org.broadinstitute.sting.utils.genotype.glf.GLFRecord;
import org.broadinstitute.sting.utils.genotype.glf.SinglePointCall;
import org.broadinstitute.sting.utils.genotype.glf.VariableLengthCall;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;


/**
 * @author aaron
 *         <p/>
 *         Class RodGLF
 *         <p/>
 *         the rod class for GLF data.
 */
public class RodGLF extends BasicReferenceOrderedDatum {
    public GLFRecord mRecord;
    private GenomeLoc mLoc;
    private static GLFRODIterator mWrap;
    private String contig;
    private int contigLength;

    public RodGLF(String name) {
        super(name);
    }

    @Override
    public boolean parseLine(Object header, String[] parts) throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String toString() {
        final StringBuilder builder = new StringBuilder();
        builder.append(contig);
        builder.append("\t");
        builder.append(mRecord.getOffset());
        builder.append("\t");
        builder.append(mRecord.getRefBase());
        builder.append("\t");
        builder.append(mRecord.getReadDepth());
        builder.append("\t");
        builder.append(mRecord.getRmsMapQ());
        builder.append("\t");
        builder.append(mRecord.getMinimumLikelihood());

        if (mRecord.getRecordType() == GLFRecord.RECORD_TYPE.SINGLE) {
            for (double d : ((SinglePointCall) mRecord).getLikelihoods()) {
                builder.append("\t");
                builder.append(d);
            }
        } else if (mRecord.getRecordType() == GLFRecord.RECORD_TYPE.VARIABLE) {
            VariableLengthCall call = (VariableLengthCall) mRecord;
            builder.append("\t");
            builder.append(call.getLkHom1());
            builder.append("\t");
            builder.append(call.getLkHom2());
            builder.append("\t");
            builder.append(call.getLkHet());
            builder.append("\t");
            builder.append(call.getIndelLen1());
            builder.append("\t");
            builder.append(call.getIndelSeq1());
            builder.append("\t");
            builder.append(call.getIndelLen2());
            builder.append("\t");
            builder.append(call.getIndelSeq2());
        }
        return builder.toString();
    }

    @Override
    public GenomeLoc getLocation() {
        return mLoc;
    }

    public static Iterator<RodGLF> createIterator(String name, File file) {
        if (mWrap == null) {
            return new RodGLF.GLFWrapper(name, file);
        }
        return mWrap;
    }

    static void setGLFWrapper(GLFRODIterator iter) {
        mWrap = iter;
    }

    private static class GLFWrapper implements GLFRODIterator {
        private GLFReader mReader;
        private String mName;

        public GLFWrapper(String name, File file) {
            mName = name;
            mReader = new GLFReader(file);
        }

        @Override
        public boolean hasNext() {
            return (mReader.hasNext());
        }

        @Override
        public RodGLF next() {
            RodGLF glf = new RodGLF(mName);
            glf.mRecord = mReader.next();
            glf.mLoc = GenomeLocParser.createGenomeLoc(mReader.getReferenceName(), mReader.getCurrentLocation());
            glf.contig = mReader.getReferenceName();
            glf.contigLength = mReader.getReferenceLength();
            return glf;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("GLF Rods don't support the remove() function");
        }
    }
}

// for testing
interface GLFRODIterator extends Iterator<RodGLF> {
}
