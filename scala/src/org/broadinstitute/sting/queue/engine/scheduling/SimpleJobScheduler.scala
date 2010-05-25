package org.broadinstitute.sting.queue.engine.scheduling

import org.broadinstitute.sting.queue.util.ProcessUtils
import org.broadinstitute.sting.queue.QArguments
import org.jgrapht.DirectedGraph

/**
 * Runs jobs one at a time locally
 */
class SimpleJobScheduler(jobGraph: DirectedGraph[ResourceNode, ResourceEdge], qArgs: QArguments)
        extends TopologicalJobScheduler(jobGraph, qArgs) {
  protected def traversedExec(exec: ExecEdge) = {
    var commandString = exec.commandString
    logger.info(commandString)

    if (!qArgs.dryRun)
      ProcessUtils.runCommandAndWait(commandString)
  }
}
