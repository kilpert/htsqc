library(knitr)

args = commandArgs(TRUE)
rmdfile = args[1]
tsvfile = args[2]
htmlfile = args[3]

knit2html( rmdfile,
           output=htmlfile,
           envir=parent.frame()
          )

# knit2html( rmdfile,
#            output=htmlfile,
#            envir=parent.frame()
# )
