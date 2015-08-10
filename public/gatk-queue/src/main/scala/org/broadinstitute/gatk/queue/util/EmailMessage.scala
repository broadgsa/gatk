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

import org.apache.commons.mail.{MultiPartEmail, EmailAttachment}
import java.io.{IOException, FileReader, File}
import javax.mail.internet.InternetAddress
import scala.collection.JavaConversions._

/**
 * Encapsulates a message to be sent over email.
 */
class EmailMessage extends Logging {
  var from: String = _
  var to: Seq[String] = Nil
  var cc: Seq[String] = Nil
  var bcc: Seq[String] = Nil
  var subject: String = _
  var body: String = _
  var attachments: Seq[File] = Nil

  /**
   * Sends the email and throws an exception if the email can't be sent.
   * @param settings The server settings for the email.
   */
  def send(settings: EmailSettings) = {
    val email = new MultiPartEmail

    email.setHostName(settings.host)
    email.setSmtpPort(settings.port)
    email.setTLS(settings.tls)
    if (settings.ssl) {
      email.setSSL(true)
      email.setSslSmtpPort(settings.port.toString)
    }

    if (settings.username != null && (settings.password != null || settings.passwordFile != null)) {
      val password = {
        if (settings.passwordFile != null) {
          val reader = new FileReader(settings.passwordFile)
          try {
            org.apache.commons.io.IOUtils.toString(reader).replaceAll("\\r|\\n", "")
          } finally {
            org.apache.commons.io.IOUtils.closeQuietly(reader)
          }
        } else {
          settings.password
        }
      }
      email.setAuthentication(settings.username, password)
    }

    email.setFrom(this.from)
    if (this.subject != null)
      email.setSubject(this.subject)
    if (this.body != null)
      email.setMsg(this.body)
    if (this.to.size > 0)
      email.setTo(convert(this.to))
    if (this.cc.size > 0)
      email.setCc(convert(this.cc))
    if (this.bcc.size > 0)
      email.setBcc(convert(this.bcc))

    for (file <- this.attachments) {
      val attachment = new EmailAttachment
      attachment.setDisposition(EmailAttachment.ATTACHMENT)
      attachment.setPath(file.getAbsolutePath)
      attachment.setDescription(file.getAbsolutePath)
      attachment.setName(file.getName)
      email.attach(attachment)
    }

    email.send
  }

  /**
   * Tries twice 30 seconds apart to send the email.  Then logs the message if it can't be sent.
   * @param settings The server settings for the email.
   */
  def trySend(settings: EmailSettings) = {
    try {
      Retry.attempt(() => send(settings), .5)
    } catch {
      case e: RetryException=> logger.error("Error sending message: %n%s".format(this.toString), e)
    }
  }

  /**
   * Converts the email addresses to a collection of InternetAddress which can bypass client side validation,
   * specifically that the domain is specified.
   * @param addresses Seq of email addresses.
   * @return java.util.List of InternetAddress'es
   */
  private def convert(addresses: Seq[String]): java.util.List[InternetAddress] = {
    addresses.map(address => new InternetAddress(address, false))
  }

  override def toString = {
    """
    |From: %s
    |To: %s
    |Cc: %s
    |Bcc: %s
    |Subject: %s
    |
    |%s
    |
    |Attachments:
    |%s
    |""".stripMargin.format(
      this.from, this.to.mkString(", "),
      this.cc.mkString(", "), this.bcc.mkString(", "),
      this.subject, this.body,
      this.attachments.map(_.getAbsolutePath).mkString("%n".format()))
  }
}
