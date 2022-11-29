Coding for slope direction
================

Taking the secreted peptidase results, we want to code them for slope
direction (so negative slopes have one color, and positive slopes
another).

``` r
library(tidyverse)
theme_set(theme_classic())
d <- read_csv("../data/new_summary_of_all_results.csv")
# give your results files more meaningful names! You will have more results later!
glimpse(d)
```

    Rows: 80
    Columns: 25
    $ ...1              <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15…
    $ seq_id            <chr> "SRR7066492_k141_800055_1", "SRR7066492_k141_712724_…
    $ paired            <chr> "SRR7066492_k141_800055_1_VS_SRR7066493_k141_1516210…
    $ depth             <chr> "15 mbsf", "15 mbsf", "15 mbsf", "15 mbsf", "15 mbsf…
    $ DDG               <chr> "\"Negative\"", "\"Negative\"", "\"Negative\"", "\"N…
    $ total_energy      <dbl> 68.06490, 6.33556, 11.69100, 54.64490, 38.04470, 118…
    $ BackHbond         <dbl> -168.0380, -111.3070, -94.9517, -266.5790, -97.1110,…
    $ SideHbond         <dbl> -50.1989, -32.1791, -24.9081, -97.2998, -22.9511, -1…
    $ Energy_VdW        <dbl> -309.560, -175.676, -176.133, -531.509, -168.491, -8…
    $ Electro           <dbl> -12.27580, -5.39436, -3.85522, -35.50950, -6.01924, …
    $ Energy_SolvP      <dbl> 423.167, 240.810, 229.703, 729.536, 233.353, 1155.93…
    $ Energy_SolvH      <dbl> -403.904, -233.328, -230.770, -685.529, -217.887, -1…
    $ Energy_vdwclash   <dbl> 24.72020, 9.11277, 7.91306, 55.51110, 11.35000, 91.5…
    $ energy_torsion    <dbl> 5.06507, 1.43978, 1.19308, 4.76428, 1.74607, 8.88756…
    $ backbone_vdwclash <dbl> 120.4640, 124.3500, 83.2963, 259.9240, 84.4464, 433.…
    $ Entropy_sidec     <dbl> 140.8650, 86.0928, 80.8627, 246.1220, 78.6843, 410.5…
    $ Entropy_mainc     <dbl> 416.100, 230.592, 227.956, 640.298, 225.855, 988.038…
    $ water_bonds       <chr> "./ranked_0.pdb", "./ranked_0.pdb", "./ranked_0.pdb"…
    $ helix_dipole      <dbl> -2.13474, -4.18676, -3.22516, -7.58027, -1.40288, -9…
    $ loop_entropy      <chr> "./ranked_0.pdb", "./ranked_0.pdb", "./ranked_0.pdb"…
    $ cis_bond          <dbl> -2.13474, -4.18676, -3.22516, -7.58027, -1.40288, -9…
    $ disulfide         <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    $ kn_electrostatic  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    $ Energy_Ionisation <dbl> 0.887267, 0.358804, 0.211391, 1.295270, 0.918279, 3.…
    $ Entropy_Complex   <chr> "./ranked_0.pdb", "./ranked_0.pdb", "./ranked_0.pdb"…

Let’s recode `DDG` as to whether the delta G (`total_energy`) is greater
or smaller in the deep sample compared to the surface sample. We’ll use
the tidyverse `group_by()` function to separate our data frame into
different groups, for each unique value of `paired`, and perform an
operation on each group.

I’ll write a little function that will return a character vector of
length 2, containing “negative”, “positive”, or “zero” depending on the
relative value of the folding energy in deep vs surface sediments.

``` r
code_slope <- function(depths, energies) {
  # Input length should be exactly 2: one surface and one deep protein
  if(length(depths) != 2) {
    warning(paste("length of depth input vector is ", length(depths), "; expected length is 2"))
  }
  
  # Calculate "delta delta G" and then test its sign
  # I imagine there is a nice way to abstract out the depth names, 
  # so that you can pass depths as an argument
  DDG <- energies[depths=="15 mbsf"] - energies[depths == "0.25 mbsf"]
  
  # There's probably a more elegant way to do this, but I can't think of it
  if(DDG < 0) {
    slope <- "negative"
  } else if(DDG > 0) {
    slope <- "positive"
  } else if(DDG == 0) {
    slope <- "zero"
  } else {
    slope <- NA
  }
  # The return needs to be of length 2
  slope <- rep(slope, 2)
  slope # In R, stating the return value of a function causes it to auto-print
}
```

This function relies on there being exactly two proteins for each
“pair”, so let’s check to make sure that’s true:

``` r
unpaired <- d %>%
  select(paired) %>% # no need for all those other columns
  group_by(paired) %>%
  summarise(n = n()) %>%
  filter(n != 2) 
head(unpaired)
```

    # A tibble: 2 × 2
      paired                                                    n
      <chr>                                                 <int>
    1 SRR7066492_k141_177316_1_VS_SRR7066493_k141_2351188_1     1
    2 SRR7066492_k141_440424_2_VS_SRR7066493_k141_3197379_9     1

Hmm. Not sure why there are two unpaired values, but for now we can just
remove them and apply our function.

``` r
df <- d %>%
  anti_join(unpaired, by="paired") %>% # removes all rows of d that exist in unpaired
  group_by(paired) %>%
  mutate(DDG = code_slope(depths=depth, energies=total_energy))
```

Now we can plug this into the plotting function that I wrote previously
(and which I have now hived off into its own function. I feel like it
would make sense to build a package for some of these operations).

``` r
source("../R/ddg_plot.R")
ddg_plot(df, legend.pos = "right")
```

![](slope_direction_tutorial_files/figure-gfm/unnamed-chunk-4-1.png)

This plot is a little more ambiguous than previous ones, so let’s do
some more analysis.

``` r
direction_table <- df %>%
  filter(depth=="0.25 mbsf") %>% # Just so we're not double-counting
  group_by(DDG) %>%
  summarise(n = n())
knitr::kable(direction_table)
```

| DDG      |   n |
|:---------|----:|
| negative |  29 |
| positive |  10 |

So more negatives than positives. Let’s do the paired t test. For that,
we’ll want to calculate the differences between the folding energies by
depth.

``` r
depth_difs <- df %>%
  group_by(paired, DDG) %>%
  summarise(DDG.numeric = total_energy[depth=="15 mbsf"] - total_energy[depth=="0.25 mbsf"])
ggplot(depth_difs, aes(x=DDG.numeric)) + 
  geom_density() + 
  geom_rug(aes(color=DDG)) + 
  geom_vline(xintercept=0, colour="gray50") + 
  scale_color_manual(values = c("#517C96", "#8D2048"))
```

![](slope_direction_tutorial_files/figure-gfm/unnamed-chunk-6-1.png)

These are distributed normally (ish), so it is appropriate (ish) to do a
parametric paired t-test.

``` r
df <- df %>% arrange(paired)
surf <- df %>%
  filter(depth == "0.25 mbsf") %>%
  pull(total_energy)
deep <- df %>% 
  filter(depth == "15 mbsf") %>%
  pull(total_energy)
t.test(surf, deep, paired=TRUE, var.equal = TRUE)
```


        Paired t-test

    data:  surf and deep
    t = 0.45675, df = 38, p-value = 0.6505
    alternative hypothesis: true mean difference is not equal to 0
    95 percent confidence interval:
     -28.79625  45.57632
    sample estimates:
    mean difference 
           8.390034 

Hmm. I’m not sure this is correct: almost 3 times more differences are
negative than positive, which seems reasonably unlikely if there is no
true difference between surface and deep. However it is about time to
put my kids to bed, so I can’t look into this further.
