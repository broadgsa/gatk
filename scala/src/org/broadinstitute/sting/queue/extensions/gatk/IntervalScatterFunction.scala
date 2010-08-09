package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.queue.function.scattergather.ScatterFunction

/**
 * An interval scatter function that allows the script to be swapped out.
 * The syntax of the script must be:
 * <splitIntervalsScript> <intervals_file> <split_intervals_1> [.. <split_intervals_n>]
 */
class IntervalScatterFunction extends ScatterFunction {
  @Argument(doc="Interval split script")
  var splitIntervalsScript: String = "splitIntervals.sh"

  def commandLine = "%s %s%s".format(splitIntervalsScript, originalInput, repeat(" ", scatterParts))
}
