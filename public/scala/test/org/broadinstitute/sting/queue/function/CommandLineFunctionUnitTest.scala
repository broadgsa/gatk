package org.broadinstitute.sting.queue.function

import org.testng.Assert
import org.testng.annotations.{DataProvider, Test}

// Since "protected" in Scala is subclass-only and doesn't allow package-level access, we need to
// extend a class if we want to test it
class CommandLineFunctionUnitTest extends CommandLineFunction {
  def commandLine = ""

  @DataProvider( name="formatArgumentTestData" )
  def formatArgumentDataProvider = {
    Array(Array("", "argvalue", "", true, true, "'argvalue'"),
          Array("", "argvalue", "", true, false, "argvalue"),
          Array("", "argvalue", "", false, true, "'argvalue'"),
          Array("", "argvalue", "", false, false, "argvalue"),
          Array("-arg", "argvalue", "", true, true, "'-arg' 'argvalue'"),
          Array("-arg", "argvalue", "", true, false, "-arg argvalue"),
          Array("ARGNAME=", "ARGVALUE", "", false, true, "'ARGNAME=ARGVALUE'"),
          Array("ARGNAME=", "ARGVALUE", "", false, false, "ARGNAME=ARGVALUE"),
          Array("-Xmx", "4", "G", true, true, "'-Xmx' '4' 'G'"),
          Array("-Xmx", "4", "G", true, false, "-Xmx 4 G"),
          Array("-Xmx", "4", "G", false, true, "'-Xmx4G'"),
          Array("-Xmx", "4", "G", false, false, "-Xmx4G"),
          Array("", "", "", true, true, "''"),
          Array("", "", "", true, false, ""),
          Array("", "", "", false, true, "''"),
          Array("", "", "", false, false, ""),
          Array("", null, "", true, true, ""),
          Array("", Nil, "", true, true, ""),
          Array("", None, "", true, true, ""),
          Array(null, null, null, true, true, ""),
          Array("", Some("argvalue"), "", true, true, "'argvalue'")
         )
  }

  @Test( dataProvider="formatArgumentTestData" )
  def testFormatArgument( prefix: String, param: Any, suffix: String, spaceSeparated: Boolean, escape: Boolean, expectedReturnValue: String ) {
    Assert.assertEquals(formatArgument(prefix, param, suffix, spaceSeparated, escape, "%s"),
                        expectedReturnValue)
  }

  @Test
  def testFormatArgumentCustomFormatString() {
    Assert.assertEquals(formatArgument("", "argvalue", "", true, true, "%.3s"), "'arg'")
  }

  @DataProvider( name = "requiredTestData" )
  def requiredDataProvider = {
    Array(Array("", "argvalue", "", true, true, " 'argvalue' "),
          Array("", "argvalue", "", true, false, " argvalue "),
          Array("", "argvalue", "", false, true, " 'argvalue' "),
          Array("", "argvalue", "", false, false, " argvalue "),
          Array("-arg", "argvalue", "", true, true, " '-arg' 'argvalue' "),
          Array("-arg", "argvalue", "", true, false, " -arg argvalue "),
          Array("ARGNAME=", "ARGVALUE", "", false, true, " 'ARGNAME=ARGVALUE' "),
          Array("ARGNAME=", "ARGVALUE", "", false, false, " ARGNAME=ARGVALUE "),
          Array("-Xmx", "4", "G", true, true, " '-Xmx' '4' 'G' "),
          Array("-Xmx", "4", "G", true, false, " -Xmx 4 G "),
          Array("-Xmx", "4", "G", false, true, " '-Xmx4G' "),
          Array("-Xmx", "4", "G", false, false, " -Xmx4G "),
          Array("", "", "", true, true, " '' "),
          Array("", "", "", true, false, "  "),
          Array("", "", "", false, true, " '' "),
          Array("", "", "", false, false, "  "),
          Array("", null, "", true, true, "  "),
          Array("", Nil, "", true, true, "  "),
          Array("", None, "", true, true, "  ")
         )
  }

  @Test( dataProvider="requiredTestData" )
  def testRequired( prefix: String, param: Any, suffix: String, spaceSeparated: Boolean, escape: Boolean, expectedReturnValue: String ) {
    Assert.assertEquals(required(prefix, param, suffix, spaceSeparated, escape),
                        expectedReturnValue)
  }

  @DataProvider( name = "optionalTestData" )
  def optionalDataProvider = {
    Array(Array("-arg", "argvalue", "", true, true, " '-arg' 'argvalue' "),
          Array("-arg", null, "", true, true, ""),
          Array("-arg", Nil, "", true, true, ""),
          Array("-arg", None, "", true, true, ""),
          Array("-arg", "", "", true, true, " '-arg' ")
         )
  }

  @Test( dataProvider="optionalTestData" )
  def testOptional( prefix: String, param: Any, suffix: String, spaceSeparated: Boolean, escape: Boolean, expectedReturnValue: String ) {
    Assert.assertEquals(optional(prefix, param, suffix, spaceSeparated, escape),
                        expectedReturnValue)
  }

  @DataProvider( name = "conditionalTestData" )
  def conditionalDataProvider = {
    Array(Array(true, "-FLAG", true, " '-FLAG' "),
          Array(true, "-FLAG", false, " -FLAG "),
          Array(false, "-FLAG", true, ""),
          Array(false, "-FLAG", false, ""),
          Array(true, null, true, "  "),
          Array(true, Nil, true, "  "),
          Array(true, None, true, "  "),
          Array(false, null, true, ""),
          Array(false, Nil, true, ""),
          Array(false, None, true, "")
         )
  }

  @Test( dataProvider="conditionalTestData" )
  def testConditional( condition: Boolean, param: Any, escape: Boolean, expectedReturnValue: String ) {
    Assert.assertEquals(conditional(condition, param, escape),
                        expectedReturnValue)
  }

  @DataProvider( name = "repeatTestData" )
  def repeatDataProvider = {
    Array(Array("", List("a", "bc", "d"), "", " ", true, true, " 'a' 'bc' 'd' "),
          Array("", List("a", "bc", "d"), "", " ", true, false, " a bc d "),
          Array("", List("a", "bc", "d"), "", "", true, true, " 'a''bc''d' "),
          Array("", List("a", "bc", "d"), "", "", true, false, " abcd "),
          Array("-f", List("file1", "file2", "file3"), "", " ", true, true, " '-f' 'file1' '-f' 'file2' '-f' 'file3' "),
          Array("-f", List("file1", "file2", "file3"), "", " ", true, false, " -f file1 -f file2 -f file3 "),
          Array("-f", List("file1", "file2", "file3"), "", " ", false, true, " '-ffile1' '-ffile2' '-ffile3' "),
          Array("-f", List("file1", "file2", "file3"), "", " ", false, false, " -ffile1 -ffile2 -ffile3 "),
          Array("-f", List("file1", "file2", "file3"), "", "", false, true, " '-ffile1''-ffile2''-ffile3' "),
          Array("-f", List("file1", "file2", "file3"), "", "", false, false, " -ffile1-ffile2-ffile3 "),
          Array("-f", List("file1", "file2", "file3"), "suffix", " ", true, true, " '-f' 'file1' 'suffix' '-f' 'file2' 'suffix' '-f' 'file3' 'suffix' "),
          Array("-f", List("file1", "file2", "file3"), "suffix", " ", true, false, " -f file1 suffix -f file2 suffix -f file3 suffix "),
          Array("-f", List("file1", "file2", "file3"), "suffix", " ", false, true, " '-ffile1suffix' '-ffile2suffix' '-ffile3suffix' "),
          Array("-f", List("file1", "file2", "file3"), "suffix", " ", false, false, " -ffile1suffix -ffile2suffix -ffile3suffix "),
          Array("-f", null, "", " ", true, true, ""),
          Array("-f", Nil, "", " ", true, true, "")
         )
  }

  @Test( dataProvider="repeatTestData" )
  def testRepeat( prefix: String, params: Traversable[_], suffix: String, separator: String,
                  spaceSeparated: Boolean, escape: Boolean, expectedReturnValue: String ) {
    Assert.assertEquals(repeat(prefix, params, suffix, separator, spaceSeparated, escape),
                        expectedReturnValue)
  }

  // Need to test None separately due to implicit conversion issues when using None in a TestNG data provider
  @Test
  def testRepeatNone() {
    testRepeat("", None, "", " ", true, true, "")
  }

  @DataProvider( name = "repeatWithPrefixFormattingTestData" )
  def repeatWithPrefixFormattingDataProvider = {
    Array(Array("-f", List("file1", "file2", "file3"), "", " ", true, true, (prefix: String, value: Any) => "%s:tag%s".format(prefix, value),
                " '-f:tagfile1' 'file1' '-f:tagfile2' 'file2' '-f:tagfile3' 'file3' "),
          Array("-f", List("file1", "file2", "file3"), "", " ", true, false, (prefix: String, value: Any) => "%s:tag%s".format(prefix, value),
                " -f:tagfile1 file1 -f:tagfile2 file2 -f:tagfile3 file3 "),
          Array("", List("file1", "file2", "file3"), "", " ", true, true, (prefix: String, value: Any) => "-%s".format(value),
                " '-file1' 'file1' '-file2' 'file2' '-file3' 'file3' "),
          Array("-f", null, "", " ", true, true, (prefix: String, value: Any) => "%s:tag%s".format(prefix, value),
                ""),
          Array("-f", Nil, "", " ", true, true, (prefix: String, value: Any) => "%s:tag%s".format(prefix, value),
                "")
         )
  }

  @Test( dataProvider = "repeatWithPrefixFormattingTestData" )
  def testRepeatWithPrefixFormatting( prefix: String, params: Traversable[_], suffix: String, separator: String,
                                      spaceSeparated: Boolean, escape: Boolean, formatPrefix: (String, Any) => String,
                                      expectedReturnValue: String ) {
    Assert.assertEquals(repeat(prefix, params, suffix, separator, spaceSeparated, escape, "%s", formatPrefix),
                        expectedReturnValue)
  }
}