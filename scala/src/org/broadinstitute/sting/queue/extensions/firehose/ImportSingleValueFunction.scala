package org.broadinstitute.sting.queue.extensions.firehose

import org.broadinstitute.sting.queue.function.JarCommandLineFunction
import org.broadinstitute.sting.commandline.{Input, Argument}
import java.io.File

/**
 * Runs the Firehose ImportSingleValue jar file.
 */
class ImportSingleValueFunction extends JarCommandLineFunction {
  @Argument(doc="firehose host")
  var host: String = _

  @Argument(doc="firehose port")
  var port: Int = _

  @Argument(doc="firehose domain")
  var domain: String = _

  @Argument(doc="firehose entity type")
  var entityType: String = _

  @Argument(doc="firehose entity id")
  var entityID: String = _

  @Argument(doc="firehose annotation type name", shortName="bamFHAnn", required=false)
  var annotationTypeName: String = _

  @Argument(doc="clean bam firehose security token", shortName="bamFHToken", required=false)
  var securityToken: String = _

  @Input(doc="imports the path to this file", exclusiveOf="importValueInFile")
  var importValue: File = _

  @Input(doc="imports the value contained in the file", exclusiveOf="importValue")
  var importValueInFile: File = _

  override def commandLine = super.commandLine + ("" +
          " PORT=%s HOST=%s DOMAIN=%s ENTITY_TYPE=%s" +
          " ENTITY_ID=%s ANNOTATION_TYPE_NAME=%s SECURITY_TOKEN=%s" +
          "%s%s"
          ).format(
    port, host, domain, entityType, entityID, annotationTypeName, securityToken,
    optional(" VALUE=", importValue), optional(" VALUE_FILE=", importValueInFile))
}
