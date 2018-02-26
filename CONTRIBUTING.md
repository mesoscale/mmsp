# Contributing to [MMSP: The Mesoscale Microstructure Simulation Project][_mmsp]

Thank you for spending some time with MMSP. We sincerely appreciate your
interest, and hope you will consider contributing to the documentation,
source code, and community at large. This code is open-source, and we really
appreciate input from the community. There are several ways you can contribute,
including editing or adding documentation, writing tutorials, submitting bug
reports, requesting new features, and writing code to incorporate. If you'd like
to help out, please do!

If you're looking for support or guidance, please visit our [Gitter chat][_gitter]
to talk with the developers and community members, or send an e-mail to
<trevor.keller@gmail.com>.

## Ground Rules

- Help us to maintain a considerate and supportive community by reading and
  enforcing the [Code of Conduct][_conduct]. Be nice to newcomers.
- Promote clean code by ensuring your code and documentation build with neither
  warnings nor errors using compiler flags equivalent to [GCC's][_gcc]
  `-Wall -pedantic`.
- Engage with the community by creating [issues][_issue] for changes or
  enhancements you'd like to make. Discuss ideas the open and get feedback.
- Maximize reusability by keeping code simple. Program in C if possible, and
  generally follow the [Google C++ Style Guide][_goog]. Avoid creating new
  classes if possible.
- Document new functions and classes with [Doxygen][_doxy]-compatible comments
  in the source code.

## Getting Started

Interested in helping, but unsure where to start? Consider proofreading the PDF
documentation and reporting or fixing errors! There's bound to be a typo in
there. Or look through the existing [issues][_issue] for fixes or enhancements
you're interested in. The number of comments on an issue can indicate its level
of difficulty, as well as the impact closing it will have.

If you're brand new to open source, welcome! You might benefit from the
[GitHub Help Pages][_ghhelp] and a [tutorial][_tut].

## Branching Workflow

[MMSP][_mmsp] uses [git][_git] version control with a [branching workflow][_branch],
summarized below. Following this workflow ensures a pristine `master` branch
that remains synchronized, and allows collaborators to quickly and easily review
contributions for content and quality.

The core mantra is ***never commit to `master`***. On [mesoscale/mmsp][_mmsp],
`master` is the rolling release branch, and code must only be merged into it
following [continuous integration][_ci] testing, code review, and discussion.

Changes are made on *branches*, using `master` as the trunk. Ideally, a branch
is created based on an [issue][_issue], which captures a *feature* or *hotfix*
to MMSP. To paraphrase illustrations of more involved workflows, a simple
branching model for MMSP looks like the following, with `B` indicating a branch
named `___`, `M` a commit with message `"___"`, `?` a pull request, and `D`
deletion of a branch. Note that branches are numbered sequentially, using the
number of the [issue][_issue] outlining the missing or broken functionality to
be patched.

```
            mesoscale/ origin/ origin/
              master   hotfix  feature  hotfix         feature
                :
                |
                +-----------------------------------------B issue4
                |                                         |
                +--------------------------B bug5         M "create function"
                |                          M "missing ;"  M "update docs"
                |         ?----------------+              |
 "merge bug5"   M---------+                |              |
                |         D                D              |
                +-----------------------------------------M "merge master"
                |                 ?-----------------------+
                |                 |                       M "address comments"
                |                 M-----------------------+
 "merge issue4" M-----------------+                       |
                |                 D                       D
                : 
```

## Pull Requests

*Please neither commit directly to the `master` branch, nor push directly to
mesoscale/mmsp.*

1. Create a fork of [MMSP][_mmsp] on your personal GitHub account.
2. For obvious changes, such as typos or edits to `.gitignore`, you can edit
   your fork directly in the browser, then file a [pull request][_pr].
3. Most changes begin with an [Issue][_issue]. If one does not exist, create it.
4. On your local machine, use the command line to work on your local fork of
   [MMSP][_mmsp]. Create a working branch off of `master`. If you're working,
   for example, on issue #42, *function template documentation*:
   - [Fork][_ghhelp] mesoscale/mmsp and [clone][_ghhelp] it to your local machine.
     Your fork is a `git remote` with the alias `origin` by default.
   - `git pull https://github.com/mesoscale/mmsp.git master` to update your
     local `master` branch from the main repository.
   - `git checkout -b issue42_function-template-documentation master` to
     create a local feature branch off of `master`.
   - Write your edits, updates, or feature code.
   - `git status`, `git add -u`, `git status`, and `git commit`,
     in that order, to add your changes. If something is amiss, follow
     the terminal guidance to fix it.
   - Write a concise commit message (first line), then in-depth commentary below
     using the keywords ["Addresses" or "Closes"][_ghkey] where appropriate.
   - `git push origin issue42_summarize-usage-in-pseudocode` to push the
     *working branch* to your fork on GitHub.
5. Visit GitHub and make the pull request official. Please fill out the pull
   request template and assign a reviewer &emdash; `@tkphd` if noone else suits.

Obvious fixes will be merged quickly. Enhancements may undergo a code review,
which we typically conduct using [reviewable][_review].

## Happy Coding!

[_mmsp]:    https://github.com/mesoscale/mmsp
[_branch]:  http://nvie.com/posts/a-successful-git-branching-model/
[_ci]:      https://docs.travis-ci.com/
[_conduct]: https://github.com/mesoscale/mmsp/blob/master/CODE_OF_CONDUCT.md
[_doxy]:    https://www.stack.nl/~dimitri/doxygen/manual/docblocks.html
[_gcc]:     https://gcc.gnu.org/
[_ghhelp]:  https://help.github.com/
[_ghkey]:   https://help.github.com/articles/closing-issues-using-keywords/
[_git]:     https://git-scm.com/
[_gitter]:  https://gitter.im/mesoscale/mmsp
[_goog]:    https://google.github.io/styleguide/cppguide.html
[_issue]:   https://github.com/mesoscale/mmsp/issues
[_pr]:      https://help.github.com/articles/about-pull-requests/
[_review]:  https://reviewable.io/reviews/mesoscale/mmsp
[_tut]:     https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github
