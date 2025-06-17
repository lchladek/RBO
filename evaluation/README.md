## DelftBlue

When running locally or installing packages, run this first to load R:

```
module load 2024r1
module load r/4.3.0 r-purrr r-dplyr r-ggplot2 r-glue
```

Then, install the following missing packages in the R interpreter:
```
```
install.packages("combinat")
install.packages("extraDistr")
install.packages("future.apply")
install.packages("latex2exp")
install.packages("mvtnorm")
install.packages("rio")
install.packages("retry")
install.packages("plyr")
install.packages("rlang")
install.packages("foreach")
install.packages("doParallel")
install.packages("DBI")
install.packages("RSQLite")
```

## Requirements

Run the following commands in the R interpreter to install the required packages:

```
install.packages("combinat")
install.packages("purrr")
install.packages("dplyr")
install.packages("extraDistr")
install.packages("future.apply")
install.packages("ggplot2")
install.packages("glue")
install.packages("latex2exp")
install.packages("mvtnorm")
install.packages("rio")
install.packages("retry")
install.packages("plyr")
install.packages("rlang")
install.packages("foreach")
install.packages("doParallel")
install.packages("DBI")
install.packages("RSQLite")
```

If the global package folder is unwritable (for example, you are running in DelftBlue), choose 'yes' at the prompt to install to your user home folder.
