## The 10+ Queuemandents

http://gatkforums.broadinstitute.org/gatk/discussion/8027/the-10-queuemandents

<h4>In no particular order:</h4>
<ul>
<li>Thou shalt never run multiple Queue jobs in the same directory</li>
<li>Thou shalt not run Queue in thy home directory, lest thy quota be overrun</li>
<li>Thou shalt wait patiently, and without whining while Queue compiles (unless thou useth <code>mvn -Ddisable.shadepackage verify</code> to disable the shade package build, and/or <code>-P !queue</code> to disable building Queue -- but why wouldst thou want to do that if thou hopeth to run Queue?)</li>
<li>Thou shalt use <code>val</code> over <code>var</code> whenever possible</li>
<li>Thou shalt not change the scatter count in mid-stream</li>
<li>Thou shalt use a proper IDE like IntelliJ for writing Queue scripts, as having Scala syntax checking shall greatly increase thy coding efficiency and reduce the number of questions thou shalt need to ask the team</li>
<li>Thou shalt recall that changing a Queue script and rerunning on the same inputs shall require <code>-startFromScratch</code></li>
<li>Thou shalt clean up thy <code>/tmp/</code> directory for finished jobs (I'm looking at thee, Ami)</li>
<li>Thou shalt consult thy <code>.queue/scatterGather/...*.out</code> files for job-specific logs</li>
<li>Thou shalt manufacture <code>.*.done</code> files with great caution</li>
<li>Thou shalt consult the names of thy <code>.*.fail</code> files for more information on failing jobs</li>
<li>Thou shalt recall that increasing the scatter number increases graph complexity and  thusly graph building time (200 samples x 140 scatter taketh about 20 minutes to build)</li>
<li>Thou shalt ensure the <code>queue-extensions</code> exist in thy repo by running a full <code>mvn verify</code> before loading the project for the first time in thy IDE</li>
</ul>