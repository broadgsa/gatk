### Description

_Please explain the changes you made here._
_Explain the **motivation** for making this change. What existing problem does the pull request solve?_
_Mention any issues fixed, addressed or otherwise related to this pull request, including issue numbers or hard links for issues in other repos._

----

### Checklist (never delete this)

Never delete this, it is our record that procedure was followed. If you find that for whatever reason one of the checklist points doesn't apply to your PR, you can leave it unchecked but please add an explanation below. 

#### Content
- [ ] Added or modified tests to cover changes and any new functionality
- [ ] All tests passing on Bamboo (including license tests!)

#### Form
- [ ] Rebased, squashed and reworded to produce a single commit (exemptions may be made -- sparingly) that is _uncluttered by excess lines_ (lines that just say "addressed review", "fixed" etc left over from squashing)
- [ ] Both the PR and the final rebased commit have a **concise** yet **descriptive** title (they are used when we compile version notes)
- [ ] Edited the README / documentation accordingly and got sign-off from support team (private code is exempt)

#### Review
- [ ] Suggest a reviewer or ask your team lead to suggest one
- [ ] Final (thumbsup) from the reviewer(s)

#### GATK4 Implications
- [ ] If your fix/change is applicable to GATK4 as well, and is reasonably small and self-contained (< ~50 lines or so), port the change to GATK4 and open a PR against https://github.com/broadinstitute/gatk or https://github.com/broadinstitute/gatk-protected as appropriate, or at least make a "best effort" attempt to do so.
- [ ] If your fix/change cannot yet be ported to GATK4 because the tool in question hasn't been ported yet, or has only been partially ported, or it would be difficult/burdensome to port the change, or you tried to port the change and failed, then add the ticket to our list of [GATK3 PRs to be eventually ported to GATK4](https://docs.google.com/document/d/1DjEHw57k5h0i8MZRGYPlQA3InvURKwQ7pCoi_Eigc4M/edit)

Once everything is checked off, you can go ahead and merge the PR. Don't forget to also delete the branch.
