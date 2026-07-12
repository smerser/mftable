require(methods)
require(R2HTML)

"%regex%" <- function(regex_array, str_array) {
    b <- F
    for (i in regex_array) {
        a <- regexpr(i, str_array, perl=T) > 0
        b <- a | b
    }
    return(b)
}

subset <- function(data, ..., drop.unused.levels=TRUE) {
    data <- base::subset(data, ...)
    if (drop.unused.levels) {
        if (is.list(data)) {
            data[] <- lapply(data, function(x) if (is.factor(x)) factor(x) else x)
        } else {
            data <- sapply(data, function(x) if (is.factor(x)) factor(x) else x)
        }
    }
    return(data)
}

# ── S4 Class ──────────────────────────────────────────────────────────────────
setClass('mftable', representation(t='matrix', nr='integer', nc='integer', di='vector'))

# ── Generics — defined BEFORE methods ─────────────────────────────────────────
xtl <- function(z, size="tiny", align="center", caption=NULL, loc='!htbp', hline=FALSE) {
    UseMethod("xtl")
}
xts <- function(z, dig=3) { UseMethod("xts") }
xtp <- function(z, rcsum=TRUE, dig=3) { UseMethod("xtp") }
xtr <- function(z, dig=3) { UseMethod("xtr") }
xtc <- function(z, dig=3) { UseMethod("xtc") }
xta <- function(z, dig=3, expected=FALSE, prop.c=TRUE, prop.r=TRUE, prop.t=TRUE, correct=NA) {
    UseMethod("xta")
}
xte <- function(z, dig=3) { UseMethod("xte") }
xth <- function(z, ...) { UseMethod("xth") }
xtm <- function(z, data, y_col, dig=3) { UseMethod("xtm") }
pNr <- function(z, ...) { UseMethod("pNr") }

# ── Constructor ────────────────────────────────────────────────────────────────
mftable <- function(z, ...) {

    bs <- 1
    if (inherits(z, 'mftable')) {
        return(z)
    } else if (inherits(z, c('matrix', 'table'))) {
        bs <- 2
        x <- ftable(z, ...)
    } else if (inherits(z, 'formula')) {
        x <- ftable(z, ...)
    } else if (inherits(z, 'ftable')) {
        x <- z
    } else {
        stop("Arguments not castable to 'mftable'")
    }

    while (TRUE) {
        lenc <- NULL
        for (i in attr(x, 'col.vars'))
            lenc <- c(lenc, length(unlist(i)))
        ngc  <- length(lenc)
        ccp  <- rev(cumprod(rev(lenc)))
        cv   <- names(attr(x, 'col.vars'))

        lenr <- NULL
        for (i in attr(x, 'row.vars'))
            lenr <- c(lenr, length(unlist(i)))
        ngr  <- length(lenr)
        rcp  <- rev(cumprod(rev(lenr)))
        rv   <- names(attr(x, 'row.vars'))

        if (ngc == 1 && ngr == 1) {
            z <- NULL
            for (i in 1:nrow(x))
                z <- rbind(z, as.character(x[i, ]))
            z <- cbind('', z)
            z <- cbind(unlist(attr(x, 'row.vars')), z)
            y <- c(unlist(rv), rep('', ccp[1] + bs))
            z <- rbind(y, z)
            y <- c(rep('', bs), unlist(cv), unlist(attr(x, 'col.vars')))
            z <- rbind(y, z)
            break
        }

        y <- unlist(attr(x, 'row.vars')[ngr])
        z <- cbind(y, x)
        y <- NULL
        for (i in (ngr - 1):(if (ngr == 1) 0 else 1)) {
            for (j in unlist(attr(x, 'row.vars')[i])) {
                r <- c(j, rep('', rcp[i + 1] - 1))
                y <- c(y, r)
            }
            z <- cbind(y, z)
            y <- NULL
        }

        len <- dim(z)[2]
        z   <- cbind(z[, 1:ngr], '', z[, (ngr + 1):len])
        len <- len + 1

        y <- c(rv, rep('', len - ngr))
        z <- rbind(y, z)

        if (ngc == 1) {
            y <- c(rep('', ngr), cv[ngc], unlist(attr(x, 'col.vars')[ngc]))
            z <- rbind(y, z)
            break
        }

        y <- c(rep('', ngr), cv[ngc],
               rep(unlist(attr(x, 'col.vars')[ngc]), cumprod(lenc)[ngc - 1]))
        z <- rbind(y, z)

        for (i in (ngc - 1):(if (ngc == 1) 0 else 1)) {
            y <- c(rep('', ngr), cv[i])
            for (h in 1:(cumprod(lenc)[i + 1] / (lenc[i] * lenc[i + 1])))
                for (j in unlist(attr(x, 'col.vars')[i])) {
                    r <- c(j, rep('', ccp[i + 1] - 1))
                    y <- c(y, r)
                }
            z <- rbind(y, z)
        }
        break
    }

    rownames(z) <- z[, 1]
    colnames(z) <- z[1, ]
    z <- z[-1, -1]

    z[, 1:ngr] <- format(z[, 1:ngr], justify='left')
    z <- new('mftable', t=z, nr=as.integer(ngr), nc=as.integer(ngc), di=attr(x, 'dim'))
    z
}

# ── Print methods ──────────────────────────────────────────────────────────────
print.mftable <- function(z) print(z@t, quote=FALSE, right=TRUE)

setMethod("show", "mftable",
    function(object) print(object@t, quote=FALSE, right=TRUE)
)

# ── Coercion ───────────────────────────────────────────────────────────────────
as.double.mftable <- function(z) {
    d <- dim(z@t)
    y <- z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]]
    y <- as.numeric(y)
    dim(y) <- z@di
    y
}

as.numeric.mftable <- function(z) as.double.mftable(z)

# ── Numeric row/col indices ────────────────────────────────────────────────────
ncr <- function(z) {
    d <- dim(z@t)
    list(row=(z@nc + 1):d[1], col=(z@nr + 1):d[2])
}

# ── Expected values ────────────────────────────────────────────────────────────
expected <- function(m) {
    fi.  <- m %*% rep(1., ncol(m))
    f.j  <- rep(1., nrow(m)) %*% m
    (fi. %*% f.j) / sum(fi.)
}

xte.default <- function(z, dig=3) round(expected(z), dig)

xte.mftable <- function(z, dig=3) {
    d    <- dim(z@t)
    y    <- z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]]
    y    <- as.numeric(y)
    dim(y) <- z@di
    e    <- expected(y)
    z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]] <- round(e, dig)
    z
}

xte.ftable <- function(z, dig=3) xte(mftable(z), dig)

# ── xtm: cell means of a numeric variable ─────────────────────────────────────
xtm.ftable <- function(z, data, y_col, dig=3) xtm(mftable(z), data, y_col, dig)

xtm.mftable <- function(z, data, y_col, dig=3) {
    d       <- dim(z@t)
    row_idx <- (z@nc + 1):d[1]
    col_idx <- (z@nr + 1):d[2]

    row_varnames <- trimws(rownames(z@t)[1:z@nc])
    col_varnames <- trimws(colnames(z@t)[1:z@nr])

    rn <- trimws(rownames(z@t)[row_idx])
    cn <- trimws(colnames(z@t)[col_idx])

    means_mat <- matrix('', nrow=length(row_idx), ncol=length(col_idx),
                        dimnames=list(rn, cn))

    for (i in seq_along(rn)) {
        for (j in seq_along(cn)) {
            mask <- rep(TRUE, nrow(data))
            for (rv in row_varnames)
                if (rv %in% colnames(data))
                    mask <- mask & trimws(as.character(data[[rv]])) == rn[i]
            for (cv in col_varnames)
                if (cv %in% colnames(data))
                    mask <- mask & trimws(as.character(data[[cv]])) == cn[j]
            vals <- data[[y_col]][mask]
            if (length(vals) > 0)
                means_mat[i, j] <- round(mean(vals, na.rm=TRUE), dig)
        }
    }

    z@t[row_idx, col_idx] <- means_mat
    z
}

# ── xts: counts + totals ───────────────────────────────────────────────────────
xts.default <- function(z, dig=3) {
    z <- rbind(z, colSums(z))
    round(cbind(z, rowSums(z)), dig)
}

xts.mftable <- function(z, dig=3) {
    d   <- dim(z@t)
    y   <- z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]]
    y   <- as.numeric(y)
    dim(y) <- z@di
    res <- xts(y, dig)
    rn  <- rownames(z@t)
    cn  <- colnames(z@t)
    z@t <- cbind(rbind(z@t, ''), '')
    z@t[(z@nc + 1):(d[1] + 1), (z@nr + 1):(d[2] + 1)] <- res
    rownames(z@t) <- c(rn, 'total')
    colnames(z@t) <- c(cn, 'total')
    z@di <- as.integer(c(z@di[1] + 1, z@di[2] + 1))
    z
}

xts.ftable <- function(z, dig=3) xts(mftable(z), dig)

xts.matrix <- function(z, dig=3) {
    rn <- if (is.null(rownames(z)))
        c(paste('[', seq(1, dim(z)[1]), ',]', sep=''), 'total')
    else
        c(rownames(z), 'total')
    cn <- if (is.null(colnames(z)))
        c(paste('[,', seq(1, dim(z)[2]), ']', sep=''), 'total')
    else
        c(colnames(z), 'total')
    z  <- rbind(z, colSums(z))
    z  <- cbind(z, rowSums(z))
    colnames(z) <- cn
    rownames(z) <- rn
    z
}

# ── xtr: row proportions ───────────────────────────────────────────────────────
xtr.default <- function(z, dig=3) round(prop.table(z, 1), dig)

xtr.mftable <- function(z, dig=3) {
    cn  <- colnames(z@t)
    d   <- dim(z@t)
    y   <- z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]]
    y   <- as.numeric(y)
    dim(y) <- z@di
    y   <- round(cbind(xtr(y, dig), margin.table(y, 1) / sum(y)), dig)
    z@t <- cbind(z@t, '')
    z@t[(z@nc + 1):d[1], (z@nr + 1):(d[2] + 1)] <- y
    colnames(z@t) <- c(cn, 'c.prop')
    z@di <- as.integer(c(z@di[1], z@di[2] + 1))
    z
}

xtr.ftable <- function(z, dig=3) xtr(mftable(z), dig)

# ── xtc: column proportions ────────────────────────────────────────────────────
xtc.default <- function(z, dig=3) round(prop.table(z, 2), dig)

xtc.mftable <- function(z, dig=3) {
    rn  <- rownames(z@t)
    d   <- dim(z@t)
    y   <- z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]]
    y   <- as.numeric(y)
    dim(y) <- z@di
    y   <- round(rbind(xtc(y, dig), margin.table(y, 2) / sum(y)), dig)
    z@t <- rbind(z@t, '')
    z@t[(z@nc + 1):(d[1] + 1), (z@nr + 1):d[2]] <- y
    rownames(z@t) <- c(rn, 'r.prop')
    z@di <- as.integer(c(z@di[1] + 1, z@di[2]))
    z
}

xtc.ftable <- function(z, dig=3) xtc(mftable(z), dig)

# ── xtp: proportion of total ───────────────────────────────────────────────────
xtp.default <- function(z, rcsum=TRUE, dig=3) round(xts(z) / sum(z), dig)

xtp.mftable <- function(z, rcsum=FALSE, dig=3) {
    d   <- dim(z@t)
    y   <- z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]]
    y   <- as.numeric(y)
    dim(y) <- z@di
    res <- xtp(y, dig)
    if (!rcsum) {
        rd  <- dim(res)
        res <- res[-rd[1], -rd[2]]
        z@t[(z@nc + 1):d[1], (z@nr + 1):d[2]] <- res
        return(z)
    }
    rn  <- rownames(z@t)
    cn  <- colnames(z@t)
    z@t <- cbind(rbind(z@t, ''), '')
    z@t[(z@nc + 1):(d[1] + 1), (z@nr + 1):(d[2] + 1)] <- res
    rownames(z@t) <- c(rn, 'total')
    colnames(z@t) <- c(cn, 'total')
    z@di <- as.integer(c(z@di[1] + 1, z@di[2] + 1))
    z
}

xtp.ftable <- function(z, rcsum=TRUE, dig=3) xtp(mftable(z), rcsum, dig)

# ── xta: full summary ──────────────────────────────────────────────────────────
xta.mftable <- function(z, dig=3, expected=FALSE, prop.c=TRUE, prop.r=TRUE,
                        prop.t=TRUE, correct=NA) {
    s <- xts(z, dig)

    cat("   Cell Contents\n|-----------------|\n|               N |\n")
    if (expected) cat("|        expected |\n")
    if (prop.r)   cat("|   N / Row Total |\n")
    if (prop.c)   cat("|   N / Col Total |\n")
    if (prop.t)   cat("| N / Table Total |\n")
    cat("|-----------------|\n\n\n")

    if (!expected && !prop.c && !prop.r && !prop.t) {
        wdt  <- max(nchar(s@t), nchar(rownames(s@t)), nchar(colnames(s@t)))
        s@t  <- formatC(s@t, width=wdt)
        y    <- rbind(colnames(s@t), s@t)
        y    <- cbind(rownames(y), y)
        for (i in 1:dim(y)[1]) {
            cat(formatC(y[i, ], width=wdt, flag='-'), sep='|', collapse='\n')
            cat(paste(rep('-', dim(y)[2] * (wdt + 1))), '\n', sep='')
        }
    } else {
        e   <- xts(xte(z, dig))
        r   <- xtr(z, dig)
        cc  <- xtc(z, dig)
        wdt <- max(if (expected) nchar(e@t) else nchar(r@t),
                   nchar(rownames(e@t)), nchar(colnames(e@t))) + 1
        cn  <- colnames(s@t)
        cat(formatC(c(' ', cn[1:z@nr]), width=wdt + 1, flag='-'))
        cat('\n', rep('-', (wdt + 2) * ncol(e@t)), '\n', sep='')
        cn[z@nr] <- ''
        colnames(s@t) <- cn

        if (prop.t && !prop.r && !prop.c && !expected) {
            p   <- xtp(z, rcsum=TRUE, dig)
            res <- list(s=s@t, e=e@t, p=p@t, c=cbind(cc@t, ''), r=rbind(r@t, ''))
        } else {
            p   <- xtp(z, rcsum=FALSE, dig)
            res <- list(s=s@t, e=e@t, p=cbind(rbind(p@t, ''), ''),
                        c=cbind(cc@t, ''), r=rbind(r@t, ''))
        }

        d <- dim(s@t)
        for (i in 1:d[1]) {
            for (j in 1:d[2]) {
                if (j == 1)
                    cat(formatC(rownames(s@t)[i], width=wdt, flag='-'), collapse='|')
                else if (i == 1)
                    cat(formatC(colnames(s@t)[j], width=wdt, flag='-'), collapse='|')
                else
                    cat(formatC(res$s[i, j], width=wdt), collapse='|')
            }
            cat('\n')
            if (expected) {
                for (a in 1:d[2])
                    cat(formatC(paste(res$e[i, a]), width=wdt), collapse='|')
                cat('\n')
            }
            if (prop.r) {
                for (a in 1:d[2])
                    cat(formatC(paste(res$r[i, a]), width=wdt), collapse='|')
                cat('\n')
            }
            if (prop.c) {
                for (a in 1:d[2])
                    cat(formatC(paste(res$c[i, a]), width=wdt), collapse='|')
                cat('\n')
            }
            if (prop.t) {
                for (a in 1:d[2])
                    cat(formatC(paste(res$p[i, a]), width=wdt), collapse='|')
                cat('\n')
            }
            cat(rep('-', (wdt + 2) * ncol(e@t)), '\n', sep='')
        }
    }

    if (!is.na(correct)) {
        correct <- as.logical(correct)
        print(chisq.test(as.double(z), correct=correct))
    }
}

xta.ftable <- function(z, dig=3, expected=FALSE, prop.c=TRUE, prop.r=TRUE,
                       prop.t=TRUE, correct=NA) {
    xta(mftable(z), dig, expected, prop.c, prop.r, prop.t, correct)
}

xta.table <- function(z, dig=3, expected=FALSE, prop.c=TRUE, prop.r=TRUE,
                      prop.t=TRUE, correct=NA) {
    xta(mftable(z), dig, expected, prop.c, prop.r, prop.t, correct)
}

xta.matrix <- function(z, dig=3, expected=FALSE, prop.c=TRUE, prop.r=TRUE,
                       prop.t=TRUE, correct=NA) {
    d  <- dim(z)
    dn <- dimnames(z)
    z  <- ftable(z)
    attr(z, 'col.vars') <- list(col.vars=if (is.null(dn)) letters[1:d[2]] else dn[2])
    attr(z, 'row.vars') <- list(row.vars=if (is.null(dn)) letters[(d[2]+1):(d[2]+d[1])] else dn[1])
    xta(z, dig, expected, prop.c, prop.r, prop.t, correct)
}

# ── xth: HTML output ───────────────────────────────────────────────────────────
xth.ftable  <- function(z, ...) xth(mftable(z), ...)

xth.mftable <- function(z, caption='', ...) {
    a <- cbind(attr(z@t, 'dimnames')[[1]], z@t)
    a <- rbind(c('', attr(z@t, 'dimnames')[[2]]), a)
    rownames(a) <- NULL
    colnames(a) <- NULL
    a[1:z@nc, z@nr + 1] <- paste('<b>', a[1:z@nc, z@nr + 1], '</b>', sep='')
    a[z@nc + 1, 1:z@nr] <- paste('<b>', a[z@nc + 1, 1:z@nr], '</b>', sep='')
    z@t <- a
    a   <- capture.output(HTML(z@t, file='', caption=caption))
    a   <- paste(a, collapse=' ')
    a   <- gsub('\t', '', a)
    a   <- gsub('</tr>', '</tr><tr>', a)
    a   <- gsub('total', '<b>total</b>', a)
    d   <- dim(z@t)
    i   <- 0
    while (i < z@nc * d[2]) {
        a <- sub('<td class=cellinside>', '<td class=_frow>', a)
        i <- i + 1
    }
    sub(' class=dataframe> ',
        paste0(' class=dataframe><col span="', z@nr + 1, '" />', sep=''), a)
}

# ── xtl: LaTeX output ─────────────────────────────────────────────────────────
xtl.mftable <- function(z, size='tiny', align='center', caption=NULL,
                        loc='!htbp', hline=FALSE) {
    d  <- dim(z@t)
    rn <- rownames(z@t)
    cat('\\begin{table}[', loc, ']\n\\begin{', align, '}\n\\begin{', size, '}\n', sep='')
    cat('\\begin{tabular}{', c(rep("l", z@nr + 1), rep("r", d[2])), '}\n', sep='')
    cat(c('', colnames(z@t)), sep='&')
    cat('\\\\\n')
    for (i in 1:d[1]) {
        cat(rn[i], z@t[i, ], sep=c(rep("&", d[2]), ''))
        cat('\\\\\n')
        if (i == 1 && hline) cat('\\hline ')
    }
    cat('\\end{tabular}\n')
    if (!is.null(caption)) cat('\\caption{', caption, '}\n', sep='')
    cat('\\end{', size, '}\n\\end{', align, '}\n\\end{table}\n', sep='')
}

xtl.matrix <- function(z, size="tiny", align="center", caption=NULL,
                       loc='!htbp', prn=TRUE) {
    d  <- dim(z)[2]
    nr <- nrow(z)
    nc <- ncol(z)
    rn <- rownames(z)
    cn <- colnames(z)
    cat('\\begin{table}[', loc, ']\n\\begin{', align, '}\n\\begin{', size, '}\n', sep='')
    cat("\\begin{tabular}{l", rep("r", d), "}\n", sep="")
    if (!is.null(cn)) {
        cat('', cn, sep='&')
        cat("\\\\\n")
    }
    for (r in 1:nr) {
        if (!is.null(rn) && prn) cat(rn[r])
        for (cc in 1:nc)
            cat('&', z[r, cc], sep='')
        cat("\\\\\n")
    }
    cat("\\end{tabular}\n")
    if (!is.null(caption)) cat("\\caption{", caption, "}\n", sep="")
    cat("\\end{", size, "}\n\\end{", align, "}\n\\end{table}\n", sep="")
}

xtl.ftable <- function(z, size="tiny", align="center", caption=NULL,
                       loc='!htbp', hline=FALSE) {
    xtl(mftable(z), size, align, caption, loc, hline)
}

# ── pNr: print wide tables in multiple rows ────────────────────────────────────
pNr.mftable <- function(z, nrow, wdt, vspace='5mm') {
    if (nrow > 1) {
        len   <- dim(z@t)[2]
        delta <- len %/% nrow
        begin <- 2
        end   <- delta
        for (i in 1:(nrow - 1)) {
            print(formatC(z@t[, c(1, begin:end)], width=wdt), quote=FALSE, right=TRUE)
            cat('\\vspace{', vspace, '}\n', sep='')
            begin <- begin + delta - 1
            end   <- end   + delta - 1
        }
        print(formatC(z@t[, c(1, begin:len)], width=wdt), quote=FALSE, right=TRUE)
    } else {
        print(z@t)
    }
}

# ── ltx: alternative LaTeX output (Dieter Menne) ──────────────────────────────
ltx <- function(z, title="", caption="", label="", loc="!htbp") {
    ft       <- format(z, quote=FALSE)
    cv       <- attr(z, "col.vars")
    rv       <- attr(z, "row.vars")
    nr       <- nrow(ft)
    nc       <- ncol(ft)
    ncolvars <- length(cv)
    nrowvars <- length(rv)

    ft[ncolvars, 1:nrowvars] <- ft[ncolvars + 1, 1:nrowvars]

    align1 <- paste(rep("l", ncolvars), collapse="")
    align2 <- paste(rep("r", nc - ncolvars), collapse="")

    cat("\\ctable[ caption={", caption, "}, label=", label, ",pos=", loc,
        ", botcap]{", align1, align2, "}{} \n{\\FL\n", sep="")

    for (i in 1:ncolvars) {
        head <- paste("\\multicolumn{1}{c}{", ft[i, ], "}", collapse="&\n", sep="")
        if (i == ncolvars) cat(head, "\n\\ML\n") else cat(head, "\n\\NN\n")
    }

    for (i in (ncolvars + 2):nr) {
        cat(paste(ft[i, ], collapse="&"))
        if (i != nr) {
            if (substr(ft[i, 1], 1, 1) == ' ' && substr(ft[i + 1, 1], 1, 1) != ' ')
                cat("\\ML\n")
            else
                cat("\\NN\n")
        }
    }
    cat("\n\\LL\n}\n")
}
