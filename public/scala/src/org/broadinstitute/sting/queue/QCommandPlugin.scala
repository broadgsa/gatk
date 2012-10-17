package org.broadinstitute.sting.queue

import engine.QStatusMessenger
import util.RemoteFileConverter

trait QCommandPlugin {
  def statusMessenger: QStatusMessenger = null
  def remoteFileConverter: RemoteFileConverter = null
}
