package org.broadinstitute.sting.queue.extensions.gatk

/**
 * Splits intervals by contig instead of evenly.
 */
class ContigScatterFunction extends IntervalScatterFunction {
  splitIntervalsScript = "splitIntervalsByContig.py"
}
