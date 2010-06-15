package org.broadinstitute.sting.queue.engine

/**
 * Represents a state between QFunctions the directed acyclic QGraph
 */
case class QNode (private val items: Set[Any])
