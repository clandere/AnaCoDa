## Test environments
* ubuntu 16.04 LTS, R 3.4.1
* ubuntu 14.04 LTS (on Travis-CI), R 3.4.1
* win-builder (devel, release, and oldrelease)

## R CMD check results
There were no ERRORs or WARNINGs. 

NOTEs:

On Ubuntu

* checking installed package size ... NOTE
  installed size is 23.8Mb
  sub-directories of 1Mb or more:
    libs  23.6Mb
 
  Response: The package contains mostly compiled code to improve computation speed. The compiled objects can not be reduced in size.

* checking top-level files ... NOTE
  Non-standard file/directory found at top level:
    'cran-comments.md'

  Response: cran-comments.md included for cran maintainers as this is a new submission.


On Windows

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Cedric Landerer <cedric.landerer@gmail.com>'

  New submission

  Possibly mis-spelled words in DESCRIPTION:
    Codon (3:20)
    Stationarity (3:37)
    checkpointing (21:139)
    codon (21:64)

  Response: All words are correclty spelled and Maintainer is Cedric Landerer

* checking top-level files ... NOTE
  Non-standard file/directory found at top level:
    'cran-comments.md'

  Response: cran-comments.md included for cran maintainers as this is a new submission.


## Downstream dependencies
There are no downstream depencies as this is a first submission.
