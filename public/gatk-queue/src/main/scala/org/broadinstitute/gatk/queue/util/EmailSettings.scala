/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.util

import java.io.File
import org.broadinstitute.gatk.utils.commandline.Argument

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

  @Argument(doc="Email SMTP password file. Defaults to none.", shortName="emailPassFile", fullName="emailPasswordFile", required=false)
  var passwordFile: File = _

  @Argument(doc="Email SMTP password. Defaults to none. Not secure! See emailPassFile.", shortName="emailPass", fullName="emailPassword", required=false)
  var password: String = _
}
