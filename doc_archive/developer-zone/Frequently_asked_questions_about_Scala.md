## Frequently asked questions about Scala

http://gatkforums.broadinstitute.org/gatk/discussion/1315/frequently-asked-questions-about-scala

<h3>1. What is Scala?</h3>
<p>Scala is a combination of an object oriented framework and a functional programming language. For a good introduction see the free online book <a href="http://programming-scala.labs.oreilly.com/">Programming Scala</a>.</p>
<p>The following are extremely brief answers to frequently asked questions about Scala which often pop up when first viewing or editing QScripts. For more information on Scala there a multitude of resources available around the web including the <a href="http://www.scala-lang.org/">Scala home page</a> and the online <a href="http://www.scala-lang.org/api/2.8.1/index.html">Scala Doc</a>.</p>
<h3>2. Where do I learn more about Scala?</h3>
<ul>
<li><a href="http://www.scala-lang.org">http://www.scala-lang.org</a></li>
<li><a href="http://programming-scala.labs.oreilly.com">http://programming-scala.labs.oreilly.com</a></li>
<li><a href="http://www.scala-lang.org/docu/files/ScalaByExample.pdf">http://www.scala-lang.org/docu/files/ScalaByExample.pdf</a></li>
<li><a href="http://devcheatsheet.com/tag/scala/">http://devcheatsheet.com/tag/scala/</a></li>
<li><a href="http://davetron5000.github.com/scala-style/index.html">http://davetron5000.github.com/scala-style/index.html</a></li>
</ul>
<h3>3. What is the difference between <code>var</code> and <code>val</code>?</h3>
<p><code>var</code> is a value you can later modify, while <code>val</code> is similar to <code>final</code> in Java.</p>
<h3>4. What is the difference between Scala collections and Java collections? / Why do I get the error: type mismatch?</h3>
<p>Because the GATK and Queue are a mix of Scala and Java sometimes you'll run into problems when you need a Scala collection and instead a Java collection is returned.</p>
<pre><code class="pre_md">   MyQScript.scala:39: error: type mismatch;
     found   : java.util.List[java.lang.String]
     required: scala.List[String]
        val wrapped: List[String] = TextFormattingUtils.wordWrap(text, width)</code class="pre_md"></pre>
<p>Use the implicit definitions in <code>JavaConversions</code> to automatically convert the basic Java collections to and from Scala collections.</p>
<pre><code class="pre_md">import collection.JavaConversions._</code class="pre_md"></pre>
<p>Scala has a very rich collections framework which you should take the time to enjoy. One of the first things you'll notice is that the default Scala collections are immutable, which means you should treat them as you would a String. When you want to 'modify' an immutable collection you need to capture the result of the operation, often assigning the result back to the original variable.</p>
<pre><code class="pre_md">var str = "A"
str + "B"
println(str) // prints: A
str += "C"
println(str) // prints: AC

var set = Set("A")
set + "B"
println(set) // prints: Set(A)
set += "C"
println(set) // prints: Set(A, C)</code class="pre_md"></pre>
<h3>5. How do I append to a list?</h3>
<p>Use the <code>:+</code> operator for a single value.</p>
<pre><code class="pre_md">  var myList = List.empty[String]
  myList :+= "a"
  myList :+= "b"
  myList :+= "c"</code class="pre_md"></pre>
<p>Use <code>++</code> for appending a list.</p>
<pre><code class="pre_md">  var myList = List.empty[String]
  myList ++= List("a", "b", "c")</code class="pre_md"></pre>
<h3>6. How do I add to a set?</h3>
<p>Use the <code>+</code> operator.</p>
<pre><code class="pre_md">  var mySet = Set.empty[String]
  mySet += "a"
  mySet += "b"
  mySet += "c"</code class="pre_md"></pre>
<h3>7. How do I add to a map?</h3>
<p>Use the <code>+</code> and <code>-&gt;</code> operators.</p>
<pre><code class="pre_md">  var myMap = Map.empty[String,Int]
  myMap += "a" -&gt; 1
  myMap += "b" -&gt; 2
  myMap += "c" -&gt; 3</code class="pre_md"></pre>
<h3>8. What are Option, Some, and None?</h3>
<p>Option is a Scala generic type that can either be some generic value or <code>None</code>. Queue often uses it to represent primitives that may be null.</p>
<pre><code class="pre_md">  var myNullableInt1: Option[Int] = Some(1)
  var myNullableInt2: Option[Int] = None</code class="pre_md"></pre>
<h3>9. What is _ / What is the underscore?</h3>
<p><a href="http://blog.normation.com/2010/07/01/scala-dreaded-underscore-psug/">François Armand</a>'s slide deck is a good introduction: <a href="http://www.slideshare.net/normation/scala-dreaded">http://www.slideshare.net/normation/scala-dreaded</a></p>
<p>To quote from his slides:</p>
<pre><code class="pre_md">Give me a variable name but
- I don't care of what it is
- and/or
- don't want to pollute my namespace with it</code class="pre_md"></pre>
<h3>10. How do I format a String?</h3>
<p>Use the <code>.format()</code> method.</p>
<p>This Java snippet:</p>
<pre><code class="pre_md">String formatted = String.format("%s %i", myString, myInt);</code class="pre_md"></pre>
<p>In Scala would be: </p>
<pre><code class="pre_md">val formatted = "%s %i".format(myString, myInt)</code class="pre_md"></pre>
<h3>11. Can I use Scala Enumerations as QScript @Arguments?</h3>
<p>No. Currently Scala's <code>Enumeration</code> class does not interact with the Java reflection API in a way that could be used for Queue command line arguments. You can use Java <code>enum</code>s if for example you are importing a Java based walker's <code>enum</code> type.</p>
<p>If/when we find a workaround for Queue we'll update this entry. In the meantime try using a String.</p>