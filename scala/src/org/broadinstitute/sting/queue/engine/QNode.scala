package org.broadinstitute.sting.queue.engine

import java.io.File

/**
 * Represents a state between QFunctions the directed acyclic QGraph
 * @param files The set of files that represent this node state.
 */
case class QNode (val files: Set[File])
