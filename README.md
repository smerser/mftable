#  R package

## mftable

## Easy print ftable object to LaTex or HTML

```
Format ftable object (flat contingency table, crosstabulation, pivot) into a LaTex or HTML rapport
```

Cell content:

> count (default) 

> row- and column proportion

> expected values

> row- and column totals

Two extra functions included:

1. `%regex%` operator: multiple regular expressions matching on multiple string arrays; i.e. returns strings matching one or more of the regex's.

2. `subset` function: automatic refactoring all factored variables in the subset, i.e. excludes non-used levels 

# Install package
```
git clone https://github.com/smerser/mftable.git

R CMD INSTALL --with-keep.source mftable
```
