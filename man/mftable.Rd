\name{mftable}
        \alias{mftable}
        \alias{xts}
        \alias{xte}
        \alias{xtr}
        \alias{xtc}
        \alias{xtp}
        \alias{xta}
        \alias{xtl}
	    \alias{xth}
	    \alias{ltx}
        \alias{expected}
        \alias{subset}
        \alias{\%regex\%}
        \title{CrossTabulation - extending ftable}
        \description{
          Crosstabulate summaries from \code{ftable} object (i.e. a flat contingency table).
          \code{xts, xtr, xtc} gives counts, row and column proportions.
          \code{xte} expected values
          \code{xtp} cell proportion of total
          \code{xta} prints a summary crosstabulation
          \code{glm} pretty overload stats::glm prettyfying output
	      \code{xtl} \code{ftable} to LaTex
	      \code{xtl.default} \code{matrix} to LaTex
	      \code{xth} \code{ftable} to HTML
	      \code{subset} allow arg: drop.unused.levels
	      \code{\%regex\%} match regular expressions within string array
          \tabular{r}{
                            N\cr
                     Expected\cr
                N / Row Total\cr
                N / Col Total\cr
                N / Table Total
          }
          \code{expected} calculate expected cell frequencies from a array, matrix or ftable.
        }
        \usage{
        mftable(z, ...)
        xta(z, dig=3, expected=FALSE, prop.c=TRUE, prop.r=TRUE, prop.t=TRUE, correct=NA)
        xtr(z, dig=3)
        xtc(z, dig=3)
        xtp(z, rcsum=TRUE, dig=3)
        xte(z, dig=3)
        xtl(z, size='tiny', align='center', caption=NULL, loc='!htbp', hline=FALSE)
        xtl.default(m, size='tiny', align='center', caption=NULL, loc='!htbp', prn=TRUE)
        xth(z, ...)
        xts(z, dig=3)
        ltx(z, title="", caption="", label="", loc="!htbp")
        expected(m)
        subset(data, ..., drop.unused.levels=TRUE)
        regex_array \%regex\% str_array
        }
        \arguments{
          \item{z}{ A \code{ftable} or \code{mftable} object.}
          \item{dig}{ Number of digits.}
          \item{expected}{ Boolen, show expected values.}
          \item{prop.r}{ Boolen, show row proportions.}
          \item{prop.c}{ Boolen, show column proportions.}
          \item{prop.t}{ Boolen, show cell proportion of total.}
          \item{correct}{ Boolen, if given show chisq.test with continuity correction as specified.}
          \item{rcsum}{ Boolen, show row and column totals.}
          \item{size}{ Font size, i.e. any valid LaTeX specification.}
          \item{align}{'left', 'center', 'right'.}
          \item{caption}{ Subtitel of LaTeX or HTML tabel.}
          \item{loc}{LaTeX table placement.}
          \item{hline}{ Boolean, draw horizontal lines.}
          \item{title}{ String, title of table.}
          \item{prn}{ Boolean, show rownames.}
          \item{m}{ 2D matrix or mftable object}
	      \item{...}{ Extra arguments for mftable, R2HTML::HTML, base::subset}
	      \item{drop.unused.levels}{ Boolen whether or not to refactorize}
	      \item{regex_array}{ Array of one or more regular expressions}
	      \item{str_array}{ Array of one or more strings}
	      \item{data}{ DataFrame}
        }
        \value{
	      \code{mftable} A mftable object.
          \code{ltx} void.
          \code{expected} matrix of expected cell frequecies.
        }
        \author{Soren Merser \email{soren.merser@gmail.com}}
        \seealso{
          \code{\link[stats]{ftable}},
          \code{\link[R2HTML]{HTML}}
        }
        \examples{
        # library(mftable)
        # A simple matrix
        # m <- matrix(c(41, 216, 64, 180), byrow=TRUE, nrow=2)
        # expected(m)
        # xta(m, prop.c=FALSE, prop.r=FALSE, prop.t=FALSE, ex=TRUE)

        # A ftable object
        # x <- ftable(Titanic, row.vars = 1:3)
        # xta(x, dig = 2, correct = TRUE)


        # htmlHead  <- '<html>
        #                   <head>
        #                       <meta charset="utf-8" />
        #                       <link rel=stylesheet type=text/css href=xtab.css>
        #                   </head>
        #               <body>'
        # htmlTail  <- "</body></html>"
        # htmlTbl   <-  xth(x, caption='Class, Sex and Age vs. Survived')
        # write(paste(
        #            htmlHead,
        #            htmlTbl,
        #            htmlTail
        #               ),
        #       file='c://Programmer/Apache Group/Apache2/htdocs/test.html'
        #       )
        }
        \keyword{Utilities}
