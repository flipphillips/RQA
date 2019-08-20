# Change Log

## The RQA Package

![logo](RQA/Documentation/icon.png)

<!--
## Types of changes

- `Added` for new features.
- `Changed` for changes in existing functionality.
- `Deprecated` for soon-to-be removed features.
- `Removed` for now removed features.
- `Fixed` for any bug fixes.
- `Security` in case of vulnerabilities. 
-->

## [Unreleased changes]

A major revision that uses the definitions from the Webber book, sparse arrays, and other tools instead of the brute force stuff we were using before.

### Changed

- `README.md` has installation instructions.
- `README.md` has documentation png.
- Routines refactored again for parameter estimation.
- `RQARecurrence` has two versions - one that operates on the whole system and another with RR_star that operates on a single diagonal.
- Helpers renamed / beginning of refactor for upper-triangular vs whole-array calculations
- Symmetric options (see above)
- `dmin` / segment-length threshold moved to an option rather than a positional parameter
- `LOI` option added to `xxxMap` routines
- **Important** `RQARecurrenceMap` now can take _either_ a timeseries _or_ an already-computed distance map.
- **Important** `RQARecurrenceMap` returns a `SparseArray` of distances rather than binary 0/1
- Significant change in determining 'occupancy' of a given location in the RM.
- semicolons after stuff that creates output in RQA.wl because I hate that.

### Added

- `RQARemoveLOI`
- `RQARatio`
- `RQADApply`
- `RQADmean`
- `RQANonzeroValues`
- added diagonal-only computations (like RR_star) so that we could calculate other difficult-to-do things. These are `RQAXxx[rm,k,opts]`

### Fixed

- `RQATrend`

***

## [0.2.0] - 2019-08-18

Released up-to-date version with changed documentation and calculation methods.

### Changed

- `DimensionalityEstimation` works but with old 'false neighbor' threshold.
- Moved `icon.png` to `Documentation` tree
- Redesigned `EstimateDimensionality` a little to handle errors better. Still not happy with it.
- Significant refactoring of parameter estimation code based on new functions added for fundamental computations.

### Added

- `EmbeddedLastTime` gives last moment in an embedded TS.
- `Dmax` and `Vmax` as per <https://github.com/alisono2>
- `PacletInfo`
- Warning message for irregular sampling
- Documentation path
- Some usage messages
- Data from RQA151 added to `Development` for debugging
- `HX` generated from RQC.EXE
- `RQA Comparison.nb` attempting to rectify differences between Webber's approach and mine.
- `VerticalLengths`, `DiagonalLengths`, etc added since they act as foundations for computations.
- Helpers for array shape, line and point finders.

### Removed

- Unnecessary `Unprotect`

***

## [0.1.1 - dev] - 2019-06-29

Semi-released to Alison Ord for playing around with.

### Fixed

- Better time series support

### Changed

- Many more options, methods, etc

### Added

- Readme, Todo, etc.

***

## [0.1.0] - 2019-06-26

Resurrected as part of the Wolfram Summer School 2019.

### Added

- initial check-in based on start of package developed with Dave Jacobs c.a. 2016.
