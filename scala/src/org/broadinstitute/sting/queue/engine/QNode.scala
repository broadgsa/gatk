package org.broadinstitute.sting.queue.engine

/**
 * Represents a state between QFunctions the directed acyclic QGraph
 */
case class QNode (val items: Set[Any])
