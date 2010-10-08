package org.broadinstitute.sting.queue.pipeline

import org.broadinstitute.sting.queue.function.CommandLineFunction

trait QPipeline {
  var commands: List[CommandLineFunction] = Nil

  def generateCommands

  def track // todo -- implement me

  def removeIntermediate // todo -- implement me
}