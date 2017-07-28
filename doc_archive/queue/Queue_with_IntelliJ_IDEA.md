## Queue with IntelliJ IDEA

http://gatkforums.broadinstitute.org/gatk/discussion/1309/queue-with-intellij-idea

<p>We have found it that Queue works best with <a href="http://www.jetbrains.com/idea/download">IntelliJ IDEA</a> Community Edition (free) or Ultimate Edition installed with the Scala Plugin enabled. Once you have downloaded IntelliJ IDEA, follow the instructions below to setup a Sting project with Queue and the Scala Plugin.</p>
<p>[[File:sting_project_libraries.png|300px|thumb|right|Project Libraries]]
[[File:sting_module_sources.png|300px|thumb|right|Module Sources]]
[[File:sting_module_dependencies.png|300px|thumb|right|Module Dependencies]]
[[File:sting_module_scala_facet.png|300px|thumb|right|Scala Facet]]</p>
<h3>1. Build Queue on the Command Line</h3>
<p>Build Queue from source from the command line with <code>ant queue</code>, so that:</p>
<ul>
<li>The lib folder is initialized including the scala jars</li>
<li>The <code>queue-extensions</code> for the GATK are generated to the build folder</li>
</ul>
<h3>2. Add the scala plugin</h3>
<ul>
<li>In IntelliJ, open the menu <code>File</code> > <code>Settings</code></li>
<li>Under the <code>IDE Settings</code> in the left navigation list select <code>Plugins</code></li>
<li>Click on the <code>Available</code> tab under plugins</li>
<li>Scroll down in the list of available plugins and install the <code>scala</code> plugin</li>
<li>If asked to retrieve dependencies, click <code>No</code>.  The correct scala libraries and compiler are already available in the lib folder from when you built Queue from the command line</li>
<li><code>Restart IntelliJ</code> to load the scala plugin</li>
</ul>
<h3>3. Creating a new Sting Project including Queue</h3>
<ul>
<li>
<p>Select the menu <code>File...</code> > <code>New Project...</code></p>
</li>
<li>
<p>On the first page of &quot;New Project&quot; select
<code>Create project from scratch</code>
Click <code>Next &gt;</code></p>
</li>
<li>
<p>On the second page of &quot;New Project&quot; select
Set the project <code>Name:</code> to <code>Sting</code>
Set the <code>Project files location:</code> to the directory where you checked out the Sting git repository, for example <code>/Users/jamie/src/Sting</code>
Uncheck <code>Create Module</code>
Click <code>Finish</code></p>
</li>
<li>
<p>The &quot;Project Structure&quot; window should open.  If not open it via the menu <code>File</code> > <code>Project Structure</code></p>
</li>
<li>
<p>Under the <code>Project Settings</code> in the left panel of &quot;Project Structure&quot; select <code>Project</code>
Make sure that <code>Project SDK</code> is set to a build of <code>1.6</code>
If the Project SDK only lists <code>&lt;No SDK&gt;</code> add a <code>New</code> > <code>JSDK</code> pointing to <code>/System/Library/Frameworks/JavaVM.framework/Versions/1.6</code></p>
</li>
<li>
<p>Under the <code>Project Settings</code> in the left panel of &quot;Project Structure&quot; select <code>Libraries</code>
Click the plus (+) to create a new Project Library
Set the <code>Name:</code> to <code>Sting/lib</code>
Select <code>Attach Jar Directories</code>
Select the path to <code>lib</code> folder under your SVN checkout</p>
</li>
<li>
<p>Under the <code>Project Settings</code> in the left panel of &quot;Project Structure&quot; select <code>Modules</code></p>
</li>
<li>
<p>Click on the <code>+</code> box to add a new module</p>
</li>
<li>
<p>On the first page of &quot;Add Module&quot; select
<code>Create module from scratch</code>
Click <code>Next \&gt;</code></p>
</li>
<li>
<p>On the second page of &quot;Add Module&quot; select
Set the module <code>Name:</code> to <code>Sting</code>
Change the <code>Content root</code> to: <code>&lt;directory where you checked out the Sting SVN repository&gt;</code>
Click <code>Next \&gt;</code></p>
</li>
<li>
<p>On the third page
Uncheck all of the other source directories only leaving the <code>java/src</code> directory checked
Click <code>Next \&gt;</code></p>
</li>
<li>
<p>On fourth page click <code>Finish</code></p>
</li>
<li>
<p>Back in the <code>Project Structure</code> window, under the <code>Module 'Sting'</code>, on the <code>Sources</code> tab make sure the following folders are selected</p>
<ul>
<li><code>Source Folders</code> (in blue):
<code>public/java/src</code>
<code>public/scala/src</code>
<code>private/java/src</code> (Broad only)
<code>private/scala/src</code> (Broad only)
<code>build/queue-extensions/src</code></li>
<li><code>Test Source Folders</code> (in green):
<code>public/java/test</code>
<code>public/scala/test</code>
<code>private/java/test</code> (Broad only)
<code>private/scala/test</code> (Broad only)</li>
</ul>
</li>
<li>
<p>In the <code>Project Structure</code> window, under the <code>Module 'Sting'</code>, on the <code>Module Dependencies</code> tab select
Click on the button <code>Add...</code>
Select the popup menu <code>Library...</code>
Select the <code>Sting/lib</code> library
Click <code>Add selected</code></p>
</li>
<li>
<p>Refresh the Project Structure window so that it becomes aware of the Scala library in <code>Sting/lib</code>
Click the <code>OK</code> button
Reopen Project Structure via the menu <code>File</code> > <code>Project Structure</code></p>
</li>
<li>In the second panel, click on the <code>Sting</code> module
Click on the plus (+) button above the second panel module
In the popup menu under <code>Facet</code> select <code>Scala</code>
On the right under <code>Facet 'Scala'</code> set the <code>Compiler library:</code> to <code>Sting/lib</code>
Click <code>OK</code></li>
</ul>
<h3>4. Enable annotation processing</h3>
<ul>
<li>Open the menu <code>File</code> > <code>Settings</code></li>
<li>Under <code>Project Settings [Sting]</code> in the left navigation list select <code>Compiler</code> then <code>Annotation Processors</code></li>
<li>Click to enable the checkbox <code>Enable annotation processing</code></li>
<li>Leave the radio button <code>obtain processors from the classpath</code> selected</li>
<li>Click <code>OK</code></li>
</ul>
<h3>5. Debugging Queue</h3>
<h4>Adding a Remote Configuration</h4>
<p>[[File:queue_debug.png|300px|thumb|right|Queue Remote Debug]]</p>
<ul>
<li>
<p>In IntelliJ 10 open the menu <code>Run</code> > <code>Edit Configurations</code>.</p>
</li>
<li>
<p>Click the gold <code>[+]</code> button at the upper left to open the <code>Add New Configuration</code> popup menu.</p>
</li>
<li>
<p>Select <code>Remote</code> from the popup menu.</p>
</li>
<li>
<p>With the new configuration selected on the left, change the configuration name from 'Unnamed' to something like 'Queue Remote Debug'.</p>
</li>
<li>
<p>Set the <code>Host</code> to the hostname of your server, and the <code>Port</code> to an unused port. You can try the default port of 5005.</p>
</li>
<li>
<p>From the <code>Use the following command line arguments for running remote JVM</code>, copy the argument string.</p>
</li>
<li>
<p>On the server, paste / modify your command line to run with the previously copied text, for example <code>java -Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=5005 Queue.jar -S myscript.scala ...</code>.</p>
</li>
<li>
<p>If you would like the program to wait for you to attach the debugger before running, change <code>suspend=n</code> to <code>suspend=y</code>.</p>
</li>
<li>Back in IntelliJ, click <code>OK</code> to save your changes.</li>
</ul>
<h4>Running with the Remote Configuration</h4>
<ul>
<li>Ensure <code>Queue Remote Debug</code> is selected via the configuration drop down or <code>Run</code> > <code>Edit Configurations</code>.</li>
<li>Set your breakpoints as you normally would in IntelliJ.</li>
<li>Start your program by running the full java path (with the above -Xdebug -Xrunjdwp ...) on the server.</li>
<li>In IntelliJ go to the <code>Run</code> > <code>Debug</code>.</li>
</ul>
<h3>6. Binding javadocs and source</h3>
<p>From <a href="http://stackoverflow.com/questions/4145734/jdk-documentation-in-intellij-idea-on-mac-os-x">Stack overflow</a>:</p>
<h4>Add javadocs:</h4>
<p>Point IntelliJ to <a href="http://download.oracle.com/javase/6/docs/api/">http://download.oracle.com/javase/6/docs/api/</a>.<br />
Go to File -&gt; Project Structure -&gt; SDKs -&gt; Apple 1.x -&gt; DocumentationPaths, and the click specify URL.</p>
<h4>Add sources:</h4>
<p>In IntelliJ, open File -&gt; Project Structure.
Click on &quot;SDKs&quot; under &quot;Platform Settings&quot;.
Add the following path under the Sourcepath tab:
/Library/Java/JavaVirtualMachines/1.6.0_29-b11-402.jdk/Contents/Home/src.jar!/src</p>