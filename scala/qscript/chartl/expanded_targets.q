import org.broadinstitute.sting.commandline.ArgumentCollection
import org.broadinstitute.sting.queue.library.ipf.ExpandIntervals
import org.broadinstitute.sting.queue.pipeline.PipelineArgumentCollection
import org.broadinstitute.sting.queue.QScript

class expanded_targets extends QScript {
  @ArgumentCollection var args : PipelineArgumentCollection = new PipelineArgumentCollection
  @Argument(shortName="bait",doc="The list of baits associated with the target list",required=false) var baitFile : File = _

  def script = {

    val intervalExpands : List[ExpandIntervals] = (new Range(0,40,1)).toList.map( u => {
      new ExpandIntervals(args.projectIntervals,1+5*u,5,new File(System.getProperty("user.dir")+"/"+args.projectName+"_expanded_%d_%d.interval_list".format(1+5*u,6+5*u)),args.projectRef,"BED")
    })

    addAll(intervalExpands)

    //  intervalExpands.map(_.outList).foreach(makeCalls(_))

  }
}