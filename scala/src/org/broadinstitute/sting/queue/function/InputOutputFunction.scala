package org.broadinstitute.sting.queue.function

import java.lang.reflect.Field
import org.broadinstitute.sting.queue.util._

/**
 * A function with @Inputs and @Outputs tagging fields that can be set by the user in a QScript
 */
trait InputOutputFunction extends QFunction with Cloneable {
  def getFieldValue(field: Field) = ReflectionUtils.getValue(this, field)
  def setFieldValue(field: Field, value: Any) = ReflectionUtils.setValue(this, field, value)

  def functionFields: List[Field] = inputFields ::: outputFields ::: internalFields
  def inputFields = ReflectionUtils.filterFields(fields, classOf[Input])
  def outputFields = ReflectionUtils.filterFields(fields, classOf[Output])
  def internalFields = ReflectionUtils.filterFields(fields, classOf[Internal])

  private lazy val fields = ReflectionUtils.getAllFields(this.getClass)
  def inputs = ReflectionUtils.getFieldNamesValues(this, inputFields)
  def outputs = ReflectionUtils.getFieldNamesValues(this, outputFields)
  def internals = ReflectionUtils.getFieldNamesValues(this, internalFields)

  /**
   * Sets a field value using the name of the field.
   * Field must be annotated with @Input, @Output, or @Internal
   * @returns true if the value was found and set
   */
  def addOrUpdateWithStringValue(name: String, value: String) = {
    fields.find(_.getName == name) match {
      case Some(field) =>
        val isInput = ReflectionUtils.hasAnnotation(field, classOf[Input])
        val isOutput = ReflectionUtils.hasAnnotation(field, classOf[Output])
        val isInternal = ReflectionUtils.hasAnnotation(field, classOf[Internal])
        if (isInput || isOutput || isInternal) {
          ReflectionUtils.addOrUpdateWithStringValue(this, field, value)
        }
        true
      case None => false
    }
  }

  def cloneFunction() = clone.asInstanceOf[this.type]
  // explicitly overriden so that trait function cloneFunction can use this.clone
  override protected def clone = super.clone

  /**
   * As the function is frozen, changes all fields to their canonical forms.
   */
  override def freeze = {
    for (field <- this.functionFields)
      mapField(field, canon)
    super.freeze
  }

  def mapField(field: Field, f: Any => Any): Any = {
    var fieldValue = this.getFieldValue(field)
    fieldValue = CollectionUtils.updated(fieldValue, f).asInstanceOf[AnyRef]
    this.setFieldValue(field, fieldValue)
    fieldValue
  }

  /**
   * Set value to a uniform value across functions.
   * The biggest example is file paths relative to the command directory in DispatchFunction
   */
  protected def canon(value: Any): Any = value
}
