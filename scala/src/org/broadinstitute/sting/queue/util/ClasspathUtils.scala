package org.broadinstitute.sting.queue.util

import collection.JavaConversions._
import org.reflections.util.ManifestAwareClasspathHelper
import java.io.File
import javax.print.URIException
import java.net.{URL, URLClassLoader}

/**
 * Builds the correct class path by examining the manifests
 */
object ClasspathUtils {

  /**
   * Returns a list of files that build up the classpath, taking into account jar file manifests.
   * @return List[File] that build up the current classpath.
   */
  def manifestAwareClassPath = {
    var urls = ManifestAwareClasspathHelper.getUrlsForManifestCurrentClasspath
    urls.map(url => try {new File(url.toURI)} catch {case urie: URIException => new File(url.getPath)})
  }

  /**
   * Adds the directory to the system class loader classpath using reflection.
   * HACK: Uses reflection to modify the class path, and assumes loader is a URLClassLoader
   * @param path Directory to add to the system class loader classpath.
   */
  def addClasspath(path: File): Unit = {
    val url = path.toURI.toURL
    val method = classOf[URLClassLoader].getDeclaredMethod("addURL", classOf[URL]);
    if (!method.isAccessible)
      method.setAccessible(true);
    method.invoke(ClassLoader.getSystemClassLoader(), url);
  }
}
