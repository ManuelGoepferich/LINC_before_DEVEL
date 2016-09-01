# LINC_CLASSES

## CLASS DEFINITION
LINCmatrix <- setClass("LINCmatrix",
                       slots       = list(
                         results     = "list", 
                         assignment  = "vector",
                         correlation = "list",
                         expression  = "matrix",
                         history     = "environment",
                         linCenvir   = "environment")
)

LINCcluster <- setClass("LINCcluster",
                        slots       = list(
                          results     = "list", 
                          assignment  = "vector",
                          correlation = "list",
                          expression  = "matrix",
                          history     = "environment",
                          linCenvir   = "environment"),
                          contains    = "LINCmatrix"
)

## CLASS DEFINITION
LINCsingle <- setClass("LINCsingle",
                       slots       = list(
                         results     = "list", 
                         assignment  = "vector",
                         correlation = "list",
                         expression  = "matrix",
                         history     = "environment",
                         linCenvir   = "environment"),
                         contains    = "LINCmatrix"
)

## CLASS DEFINITION
LINCbio <- setClass("LINCbio",
                    slots       = list(
                      results     = "list", 
                      assignment  = "vector",
                      correlation = "list",
                      expression  = "matrix",
                      history     = "environment",
                      linCenvir   = "environment"),
                      contains = "LINCmatrix"
)

## CLASS DEFINITION
LINCfeature <- setClass("LINCfeature",
                        slots        = list(
                          customID     = "character",
                          customCol    = "character",
                          setLevel     = "character",
                          showLevels   = "logical")
)

