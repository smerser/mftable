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
	    \alias{ncr}
	    \alias{as.numeric}
	    \alias{as.double}
        \alias{expected}
        \alias{subset}
        \alias{\%regex\%}

        \title{CrossTabulation - extending ftable}

        \description{
            Crosstabulate summaries from \code{ftable} object (i.e. a flat contingency table). Functions for LaTex or HTML output.
            \tabular{llr}{

            \code{xts, xtr, xtc} \tab Gives counts, row and column proportions \tab \cr

            \code{xte}  \tab Expected values \tab \cr

            \code{xtp}  \tab Cell proportion of total \tab \cr

            \code{xtl}  \tab \code{ftable} or \code{matrix} to LaTex \tab \cr

    	    \code{xth}  \tab \code{ftable} to HTML \tab \cr

            \code{xta}  \tab Prints a summary. Cell entries are: \tab \cr
                        \tab    \tab                N   \cr
                        \tab    \tab         Expected   \cr
                        \tab    \tab    N / Row Total   \cr
                        \tab    \tab    N / Col Total   \cr
                        \tab    \tab    N / Table Total \cr

            \code{ltx}  \tab LateX output. OBS: z must be a \code{ftable} object!
            }

        Additional auxillary functions:

            \code{as.numeric} and \code{as.double} returns the numerical matrix of a \code{mftable} object

            \code{ncr} return row and column indexes af the numeric information

            \code{expected} calculate expected cell frequencies from a array, matrix or ftable.

            \code{subset} has arg: drop.unused.levels.

            \code{\%regex\%} match one or more regular expressions with multiple string array.

        }

        \usage{
        mftable(z, ...)
        xta(z, dig=3, expected=FALSE, prop.c=TRUE, prop.r=TRUE, prop.t=TRUE, correct=NA)
        xtr(z, dig=3)
        xtc(z, dig=3)
        xtp(z, rcsum=TRUE, dig=3)
        xte(z, dig=3)
        xtl(z, size='tiny', align='center', caption=NULL, loc='!htbp', hline=FALSE)
        xtl(m, size='tiny', align='center', caption=NULL, loc='!htbp', prn=TRUE)
        xth(z, ...)
        xts(z, dig=3)
        ltx(z, title="", caption="", label="", loc="!htbp")
        as.numeric(z)
        as.double(z)
        ncr(z)
        expected(m)
        subset(data, ..., drop.unused.levels=TRUE)
        regex_array \%regex\% str_array
        }
        \arguments{
            \item{z}{ A \code{ftable} or \code{mftable} object.}
            \item{dig}{ Number of digits.}
            \item{expected}{ Boolen, show expected values.}
            \item{prop.r}{ Boolen, show cell proportion of row.}
            \item{prop.c}{ Boolen, show cell proportion of column.}
            \item{prop.t}{ Boolen, show cell proportion of total.}
            \item{correct}{ Boolen, if given show chisq.test with continuity correction as specified.}
            \item{rcsum}{ Boolen, show row and column totals.}
            \item{size}{ Font size, i.e. any valid LaTeX specification.}
            \item{align}{ 'left', 'center', 'right'.}
            \item{caption}{ Subtitel of LaTeX or HTML tabel.}
            \item{loc}{ LaTeX table placement.}
            \item{hline}{ Boolean, draw horizontal lines.}
            \item{label}{ String, subtitle.}
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
          \code{ltx} void; i.e. just print
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
        \keyword{Utilities, ftable, Crosstab, Pivot, Contingency table, Regular expression, subset, factor}
