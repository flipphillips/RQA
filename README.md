# RQA

![logo](RQA/Documentation/icon.png)

v 0.2.1dev

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

[![View notebooks](https://wolfr.am/Etv7EZ90)](https://wolfr.am/FFDrp9F5)

## Documentation

![dochowto](RQA/Documentation/dochowto.png)

## People

[flipphillips.com][https://flipphillips.com]
