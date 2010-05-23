package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.graph.SimpleDirectedGraph
import org.broadinstitute.sting.queue.util.ProcessUtils
import org.broadinstitute.sting.queue.engine.Pipeline

/**
 * Runs jobs one at a time locally
 */
class SimpleJobScheduler(jobGraph: SimpleDirectedGraph[ResourceNode, ResourceEdge]) extends TopologicalJobScheduler(jobGraph) {
  protected def traversedExec(exec: ExecEdge) = {
    var commandString = exec.commandString
    logger.info(commandString)

    // TODO: Pre-print the commands?
    if (!Pipeline.dryRun)
      ProcessUtils.runCommandAndWait(commandString, Map.empty[String, String])
  }
}
