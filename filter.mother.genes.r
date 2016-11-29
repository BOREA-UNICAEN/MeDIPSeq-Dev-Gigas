library(xlsx)
expr <- read.xlsx2("expressions.xlsx", 1, header = TRUE, colClasses =
                  c("character", rep("numeric",12)))

anal <- logical(7)
keep <- logical(nrow(expr))

for (i in 1:nrow(expr))
{
 for (j in 3:9)
 {
  if (expr[i,j] < expr[i,j-1]) anal[j-2]<-FALSE else anal[j-2]<-TRUE
 }
 if (any(anal)) keep[i]<-TRUE else keep[i]<-FALSE
}
