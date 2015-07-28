# calin #

__Calin__ is a calibration and analysis package for data from arrays
  of imaging atmospheric Cherenkov detectors (IACTs). Its main focus
  is on computational accuracy and robustness. Secondary considerations
  are simplicity, speed, end-user ease-of-use and maintainability.

### The name ###

__Calin__ is a (French) word that my three-year-old daughter says a
  lot. Here it isn't meant to mean anything specific, but to evoke the
  word "calibration". Since acronyms abound in astronomy it could at a
  stretch stand for "calibration and analysis for large IACT
  networks", but that wasn't the original motivation.

### Design requirements ###

**Computational accuracy and robustness:** It seems self-evident that
  this should be the primary requirement of any data analysis product, no
  matter what the domain of study. In __calin__, this is the primary
  design requirement and we expect that significant amount of the
  development time will go in this direction. This goal demands that
  algorithms used must be well understood and tested, probably with
  specific Monte Carlo studies. Wherever possible multiple algorithms
  should be applied to the data and their results compared to identify
  problems and understand systematics. It is proper that __calin__ be
  based on the algorithms that have been used by existing IACT
  experiments for decades, but our specific focus on accuracy suggests
  that they be (re)studied in some detail before being included. There
  may be room to improve even long-established methods. Since accuracy
  is more important than speed, __calin__ should adopt
  computentionally expensive methods should they provide improvements
  over simpler methods.

**Simplicity:** The software should be as simple as possible given
  that it must accomplish the main task of (accurately) analyzing IACT
  data. Needless complexity should be avoided where possible. This
  requirement is hard to achieve in practice, for two main
  reasons. First, it requires assessing the needs of the end product
  in advance, and it can be tempting to engineer the design of the
  system for a worst-case of what __might__ be needed, to avoid
  refactoring at a later stage. However refactoring is something that
  must avoided at all costs, rather its probability and ultimate cost
  should be balanced against the wasted effort of producing an
  over-engineered system that is ultimately not needed. Certainly the
  design should accomodate the types of analysis that are done now,
  and the increase in data rate that will accompany next generation
  IACT networks. Care must be taken going much beyond that.  We
  attempt to write justifications for major design decisions to help
  clarify whether they are really necessary.

**Speed:** 

**End-user ease-of-use** 

**Modularity** 


