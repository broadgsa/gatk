package org.broadinstitute.sting.gatk.iterators;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Baseclass used to describe a read transformer like BAQ and BQSR
 *
 * Read transformers are plugable infrastructure that modify read state
 * either on input, on output, or within walkers themselves.
 *
 * The function apply() is called on each read seen by the GATK (after passing
 * all ReadFilters) and it can do as it sees fit (without modifying the alignment)
 * to the read to change qualities, add tags, etc.
 *
 * Initialize is called once right before the GATK traversal begins providing
 * the ReadTransformer with the ability to collect and initialize data from the
 * engine.
 *
 * Note that all ReadTransformers within the classpath are created and initialized.  If one
 * shouldn't be run it should look at the command line options of the engine and override
 * the enabled.
 *
 * @since 8/31/12
 * @author depristo
 */
abstract public class ReadTransformer {
    /**
     * When should this read transform be applied?
     */
    private ApplicationTime applicationTime;

    /**
     * Keep track of whether we've been initialized already, and ensure it's not called more than once.
     */
    private boolean initialized = false;

    protected ReadTransformer() {}

    /**
     * Master initialization routine.  Called to setup a ReadTransform, using it's overloaded initializeSub routine.
     *
     * @param overrideTime if not null, we will run this ReadTransform at the time provided, regardless of the timing of this read transformer itself
     * @param engine the engine, for initializing values
     * @param walker the walker we intend to run
     */
    @Requires({"initialized == false", "engine != null", "walker != null"})
    @Ensures("initialized == true")
    public final void initialize(final ApplicationTime overrideTime, final GenomeAnalysisEngine engine, final Walker walker) {
        if ( engine == null ) throw new IllegalArgumentException("engine cannot be null");
        if ( walker == null ) throw new IllegalArgumentException("walker cannot be null");

        this.applicationTime = initializeSub(engine, walker);
        if ( overrideTime != null ) this.applicationTime = overrideTime;
        initialized = true;
    }

    /**
     * Subclasses must override this to initialize themselves
     *
     * @param engine the engine, for initializing values
     * @param walker the walker we intend to run
     * @return the point of time we'd like this read transform to be run
     */
    @Requires({"engine != null", "walker != null"})
    @Ensures("result != null")
    protected abstract ApplicationTime initializeSub(final GenomeAnalysisEngine engine, final Walker walker);

    /**
     * Should this ReadTransformer be activated?  Called after initialize, which allows this
     * read transformer to look at its arguments and decide if it should be active.  All
     * ReadTransformers must override this, as by default they are not enabled.
     *
     * @return true if this ReadTransformer should be used on the read stream
     */
    public boolean enabled() {
        return false;
    }

    /**
     * Has this transformer been initialized?
     *
     * @return true if it has
     */
    public final boolean isInitialized() {
        return initialized;
    }

    /**
     * When should we apply this read transformer?
     *
     * @return true if yes
     */
    public final ApplicationTime getApplicationTime() {
        return applicationTime;
    }

    /**
     * Primary interface function for a read transform to actually do some work
     *
     * The function apply() is called on each read seen by the GATK (after passing
     * all ReadFilters) and it can do as it sees fit (without modifying the alignment)
     * to the read to change qualities, add tags, etc.
     *
     * @param read the read to transform
     * @return the transformed read
     */
    @Requires("read != null")
    @Ensures("result != null")
    abstract public GATKSAMRecord apply(final GATKSAMRecord read);

    @Override
    public String toString() {
        return getClass().getSimpleName();
    }

    /**
     * When should a read transformer be applied?
     */
    public static enum ApplicationTime {
        /**
         * Walker does not tolerate this read transformer
         */
        FORBIDDEN,

        /**
         * apply the transformation to the incoming reads, the default
         */
        ON_INPUT,

        /**
         * apply the transformation to the outgoing read stream
         */
        ON_OUTPUT,

        /**
         * the walker will deal with the calculation itself
         */
        HANDLED_IN_WALKER
    }
}
