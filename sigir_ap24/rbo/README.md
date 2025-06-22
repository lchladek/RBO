# Rank-Biased Overlap implementation in R

- [`rbo.R`](./rbo.R) is a verbatim copy of the [Rank-Biased Overlap implementation](https://github.com/julian-urbano/sigir2024-rbo/blob/main/rbo/R/) from the following paper:  M. Corsi and J. Urbano, "[The Treatment of Ties in Rank-Biased Overlap](http://julian-urbano.info/files/publications/068-treatment-ties-rank-biased-overlap.pdf)", *International ACM SIGIR Conference on Research and Development in Information Retrieval*, 2024. This can calculate all relevant scores ($RBO_{EXT}$, $RBO_{MIN}$, $RBO_{MAX}$ and $RBO_{RES}$) in all three tie-aware variants ($RBO^w$, $RBO^a$ and $RBO^b$). There is also a [Python implementation](https://github.com/julian-urbano/sigir2024-rbo/tree/main/rbo/Python). 
- [`low.R`](./low.R) and [`hig.R`](./hig.R) contain functions to compute the permutations leading to the lowest and highest possible scores, from which one can compute $RBO^{low}$ and $RBO^{hig}$ variants, as defined in the following paper: M. Corsi and J. Urbano, "[How do Ties Affect the Uncertainty in Rank-Biased Overlap?](https://julian-urbano.info/files/publications/069-how-ties-affect-uncertainty-rank-biased-overlap.pdf)", *International ACM SIGIR Conference on Information Retrieval in the Asia Pacific*, 2024.
- [`avg.R`](./avg.R) contains functions to compute the $RBO^{avg}$ variants, as defined in the same paper.

If you use this software, please cite these papers.

## How to create rankings

Items are represented with strings, using only alphanumeric characters. A ranking is a list of items. Items that are tied are simply concatenated under the same vector. For example:

```r
x <- list("red", c("blue", "green"), "yellow", "pink")
y <- list(c("blue", "red"), "white", c("yellow", "black", "purple"), "green")
```

As this can become tedious, we provide three helper functions to create rankings:

- `from_string` allows you to create a ranking from a plain text representation. For example, ranking `x` above can be created as follows:

  ```r
  > from_string("red (blue green) yellow pink")
  # [[1]]
  # [1] "red"
  # 
  # [[2]]
  # [1] "blue"  "green"
  # 
  # [[3]]
  # [1] "yellow"
  # 
  # [[4]]
  # [1] "pink"
  ```
  
- `to_string` is the inverse operation; it creates a plain text representation from a ranking:

  ```r
  > to_string(x)
  # [1] "red (blue green) yellow pink"
  ```
  
- `extract_ranking` creates a ranking given a vector of `items` and their `scores`, sorting in descending order:

  ```r
  > extract_ranking(c("green", "yellow", "blue", "pink", "red"),
                    c(4, 2, 4, 1, 6))
  # [[1]]
  # [1] "red"
  # 
  # [[2]]
  # [1] "blue"  "green"
  # 
  # [[3]]
  # [1] "yellow"
  # 
  # [[4]]
  # [1] "pink"
  ```
  
  This is useful for example to extract the rankings as represented in a retrieval run.

## How to compute RBO

Given two rankings `x` and `y`, RBO can be easily computed with a `p` of your choice:

```r
> rbo(x, y, p = .95)
#       ext       min       max       res
# 0.6922853 0.3310519 0.8930692 0.5620173
```

All four scores will be computed by default. If only some of them are necessary, you can specify which through argument `score`:

```r
> rbo(x, y, p = .95, score = c("ext", "res"))
#       ext       res
# 0.6922853 0.5620173
```

The $RBO^a$ variant is computed by default, but this can be changed through argument `ties`:

```r
> rbo(x, y, p = .95, ties = "w")
#       ext       min       max       res
# 0.7068257 0.3429684 0.9049857 0.5620173

> rbo(x, y, p = .95, ties = "b")
#       ext       min       max       res
# 0.7207131 0.3509163 0.9129336 0.5620173
```

If you want to compute the $RBO^{avg}$ variant, you can use the appropriate function from `avg.R`:

```r
> rbo_avg_ext(x, y, .95)
# [1] 0.6922853

> rbo_avg_min(x, y, .95)
# [1] 0.3310519

> rbo_avg_max(x, y, .95)
# [1] 0.8930692
```

## How to compute bounds for the seen part? 

Given two rankings `x` and `y`, the arrangements of tied items that lead to the lowest and highest possible RBO scores are computed respectively as follows:

```r
> l <- low_rankings(x, y)
> l
# $x
# [1] "red"    "green"  "blue"   "yellow" "pink"  
# $y
# [1] "blue"   "red"    "white"  "black"  "purple" "yellow" "green" 

> h <- hig_rankings(x, y)
> h
# $x
# [1] "red"    "blue"   "white"  "yellow" "black"  "purple" "green" 
# $y
# [1] "red"    "blue"   "green"  "yellow" "pink"  
```

With these, we can use `rbo` to compute $RBO_{EXT}$, $RBO_{MIN}$ and $RBO_{MAX}$ as usual. Specifically, the residuals, as defined in the paper, are:

```r
# RES_U
> rbo_avg_max(x, y, .95) - rbo_avg_min(x, y, .95)
# [1] 0.5620173

# RES_S
> rbo(h$x, h$y, .95, score = "ext") - rbo(l$x, l$y, .95, score = "ext")
# [1] 0.1359071

# RES_SU
> rbo(h$x, h$y, .95, score = "max") - rbo(l$x, l$y, .95, score = "min")
# [1] 0.6546295
```

## What RBO variant should you use? 

- When a tie represents equality, so that tied items *really* occur at the same rank, you should compute $RBO^w$.
- When a tie represents uncertainty, so that it is not known which item appears first:

  - Ties should not be broken deterministically, such as by doc ID, because it inflates $RBO$ scores.
  - Ties should not be broken at random because it introduces noise. $RBO^{avg}$ should be used instead, as it precisely computes the expected $RBO$ when breaking ties at random.
  - If the measured overlap should be corrected by the amount of information lost due to ties, $RBO^b$ should be used. This ensures $RBO^b(X,X)=1$, and implies $RBO^a\leq RBO^b$.
  - Uncertainty in the $RBO$ scores should always be computed and reported when appropriate, for example via residuals $RES_U$ (unseen part), $RES_S$ (seen part) or $RES_{SU}$ (seen and unseen parts).

For more details, please refer to Section 1.1 of the SIGIR paper and Section 5 of the SIGIR-AP paper.


