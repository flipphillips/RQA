# RQA

![logo](RQA/Documentation/icon.png)

v 0.3.0dev

***

A set Mathematica tools for RQA — Recurrence quantification analysis.

## Intro

Note that I posted an early version of this to GitHub for sharing with the 2019 Wolfram Summer School. It is, as such, incomplete, undocumented, and other fun un-things. I will continue to update and upgrade, of course. I invite any and all enhancements.

## Installation

This package uses a subset of the Paclet mechanism but can be installed manually as a subdirectory to the

```mathematica
FileNameJoin[{$UserBaseDirectory, "Applications"}]
```

directory. Basically just copy the sub-directory `RQA` from the `RQA` repository to `FileNameJoin[{$UserBaseDirectory, "Applications"}]`

Alternately, there is a Wolfram Language script, `install.wls` that will do the best it can to make a symbolic link. No guarantees, of course. To try it:

* Check out distribution to a local GitHub location
* Run the enclosed `install.wls` via `wolframscript -f install.wls`

The nice thing about this method is that, when you update/pull fresh versions of `RQA` from github you'll automatically get updates installed into Mathematica. **But** I'm fairly certain this method will only work for macOS and Linux for now, since it uses symbolic links. I'd be happy to accept some more platform-independent solution.

## Repository Notebook Viewer

This lets you view notebooks in this package using the WolframCloud

* [![View notebooks](https://wolfr.am/Etv7EZ90)](https://wolfr.am/FVNv9Yfe) - master branch
* [![View notebooks](https://wolfr.am/Etv7EZ90)](https://wolfr.am/FFDrp9F5) - dev branch

## Documentation

![Documentation how-to](RQA/Documentation/dochowto.png)

## Version History

The basic code is pretty straight forward, but there are some idiosyncrasies. Back when Dave and I started this, we used some meta-references to the Webber work. It turns out that, not all of those analyses were properly described. I'll leave those out for right now to avoid both my and their embarrassment.

The most recent work is based on the book recently published by Webber et al. As we use different methods we will note the references here.

* Version 0.1 - Based on a stack of papers that were not the Webber originals, but rather the interpretative work of others.
* Version 0.2 - Based on the paper in the dynsys NSF book. Webber Jr, C. L., & Zbilut, J. P. (2005). Recurrence quantification analysis of nonlinear dynamical systems. _Tutorials in contemporary nonlinear methods for the behavioral sciences_, 94, 26-94. <https://www.nsf.gov/pubs/2005/nsf05057/nmbs/chap2.pdf>
* Version 0.3+ - Based on Webber Jr, C. L., & Marwan, N. (2015). _Recurrence quantification analysis. Theory and Best Practices._ Springer.

## People

<https://flipphillips.com>

## References

* Webber Jr, C. L., & Zbilut, J. P. (2005). Recurrence quantification analysis of nonlinear dynamical systems. _Tutorials in contemporary nonlinear methods for the behavioral sciences_, 94, 26-94. <https://www.nsf.gov/pubs/2005/nsf05057/nmbs/chap2.pdf>
* Webber Jr, C. L., & Marwan, N. (2015). _Recurrence quantification analysis. Theory and Best Practices._ Springer.
