## Test environments
* local Debian 11 (bullseye) install, R 4.1.1
* Rhub
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Debian Linux, R-devel, GCC, no long double
* win-builder, R-devel and R-release
* Github actions
  - Windows Server x64, R 4.1.1
  - macOS Catalina 10.15.7, R 4.1.1
  - Ubuntu 20.04.2 LTS, R 4.1.1
  - Ubuntu 20.04.2 LTS, R-devel (2021-08-25 r80817)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Michal Kolesár <kolesarmi@googlemail.com>’

  New submission

  Possibly misspelled words in DESCRIPTION:
    Kolesár (22:58)
    M�ller (23:14)
    Plagborg (23:5)

  Found the following (possibly) invalid URLs:
    URL: https://opportunityinsights.org/data/?paper_id=599
      From: man/cz.Rd
      Status: 403
      Message: Forbidden


  The misspelled words are the package authors' names. The URL is valid, I
  double-checked manually.

## Downstream dependencies
There are currently no downstream dependencies for this package
