package org.broadinstitute.sting.queue.util

import collection.JavaConversions._
import org.reflections.util.ManifestAwareClasspathHelper
import java.io.File
import javax.print.URIException

/**
 * Builds the correct class path by examining the manifests
 */
object ClasspathUtils {
  def manifestAwareClassPath = {
    var urls = ManifestAwareClasspathHelper.getUrlsForManifestCurrentClasspath
    var files = urls.map(url => try {new File(url.toURI)} catch {case urie: URIException => new File(url.getPath)})
    files.mkString(File.pathSeparator)
  }
}
