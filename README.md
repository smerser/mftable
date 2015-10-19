#R package.
## library(mftable): easy print ftable object to LaTex or HTML

```
Format ftable object (flat contingency table, crosstabulation, pivot) into a LaTex or HTML rapport
```

Cell content:

> count (default) 

> row- or column proportion

> expected values

Optional include row- and column totals.

Two extra functions included:

1. `%regex%` operator: multiple regular expressions matching on string arrays

2. `subset` function: automatic refactoring excluding non-used levels in the subset

# Install package
```
git clone https://github.com/smerser/mftable.git

R CMD INSTALL --with-keep.source mftable
```
