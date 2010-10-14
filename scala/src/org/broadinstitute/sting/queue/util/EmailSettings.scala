package org.broadinstitute.sting.queue.util

import java.io.File
import org.broadinstitute.sting.commandline.Argument

/**
 * Email server settings.
 */
class EmailSettings {
  @Argument(doc="Email SMTP host. Defaults to localhost.", shortName="emailHost", fullName="emailSmtpHost", required=false)
  var host = "localhost"

  @Argument(doc="Email SMTP port. Defaults to 465 for ssl, otherwise 25.", shortName="emailPort", fullName="emailSmtpPort", required=false)
  var port = 25

  @Argument(doc="Email should use TLS. Defaults to false.", shortName="emailTLS", fullName="emailUseTLS", required=false)
  var tls = false

  @Argument(doc="Email should use SSL. Defaults to false.", shortName="emailSSL", fullName="emailUseSSL", required=false)
  var ssl = false

  @Argument(doc="Email SMTP username. Defaults to none.", shortName="emailUser", fullName="emailUsername", required=false)
  var username: String = _

  @Argument(doc="Email SMTP password. Defaults to none. Not secure! See emailPassFile.", shortName="emailPass", fullName="emailPassword", required=false)
  var password: String = _

  @Argument(doc="Email SMTP password file. Defaults to none.", shortName="emailPassFile", fullName="emailPasswordFile", required=false)
  var passwordFile: File = _
}
