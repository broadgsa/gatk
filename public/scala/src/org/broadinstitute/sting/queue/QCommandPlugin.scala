package org.broadinstitute.sting.queue

import engine.QStatusMessenger
import util.RemoteFileConverter

trait QCommandPlugin {
  def statusMessenger: QStatusMessenger = null
  def remoteFileConverter: RemoteFileConverter = null
  def qScriptClass: Class[_ <: QScript] = classOf[QScript]
  def initScript(script: QScript) {}
}
