## (howto) Set up remote debugging in IntelliJ

http://gatkforums.broadinstitute.org/gatk/discussion/4712/howto-set-up-remote-debugging-in-intellij

<p>Remote debugging is a powerful tool but requires a little bit of setup. Here is the 3-step process to an easier life.</p>
<h3>1. Set up the remote config in IntelliJ</h3>
<p>Do the following in IntelliJ: </p>
<ul>
<li>
<p>Run -&gt; Edit Configurations -&gt; Add new configuration (+ symbol top left) -&gt; Remote</p>
</li>
<li>Fill in the appropriate host (<code>gsa[_machine#_].broadinstitute.org</code>) and port number (<em>XXXXX</em>), where <em>xxxxx</em> is a 5-digit port number you make up to avoid accidentally connecting to someone else's debug session. Press OK. Add breakpoint(s) where you want them in the code.</li>
</ul>
<h3>2. Run the tool on gsa machine</h3>
<p>Run the GATK command from the server with </p>
<pre>
java -agentlib:jdwp=transport=dt_socket,server=y,suspend=y,address=<i>5-digit_port_number</i> \
     -jar <i>_toolName_</i> \
     <i>args</i>
</pre>
<p>GATK will wait for IntelliJ to actually start running.</p>
<h3>3. Chase bug(s) in IntelliJ</h3>
<p>Go to IntelliJ</p>
<ul>
<li>Run -&gt; Debug -&gt; Select the configuration you just created.</li>
</ul>
<p>Now chase.</p>
<p>You can also add the <code>agentlib</code> business as an alias in your <code>.profile</code> or <code>.my.bashrc</code> on the server like I did. Boom.</p>