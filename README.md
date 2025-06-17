# A Probabilistic Model of Uncertainty in Rank-Biased Overlap for Tied Data

This repository contains the R source code for reproducing the results in the article,
'A Probabilistic Model of Uncertainty in Rank-Biased Overlap for Tied Data' (Lukáš Chládek, 2025).

# Usage

## Using the R libraries

R is required. 

```r
# Libraries required by Corsi and Urbano (2024)
install.packages("dplyr")
install.packages("extraDistr")
install.packages("future.apply")
install.packages("ggplot2")
install.packages("glue")
install.packages("latex2exp")
install.packages("mvtnorm")
install.packages("rio")

# Additional libraries required for this evaluation
install.packages("combinat")
install.packages("purrr")
install.packages("argparse")
install.packages("rlang")
install.packages("progressr")
install.packages("beepr")
```


## Reproducing the evaluation locally

This requires R and [Gourd](https://github.com/ConSol-Lab/gourd/).

Please note that R's random number generator may give different results depending on the software version.
The results in the paper use R version 4.3.0.

1. **Install R and Gourd**

   If these are not yet installed, use the instructions at [https://www.r-project.org/]() and [https://github.com/ConSol-Lab/gourd]().
2. **Install R dependencies**

   Run the installation commands from "Using the R libraries" above.
2. **Clone this repository**

   `git clone https://github.com/lchladek/RBO && cd RBO`
3. **Enter the `evaluation` directory**

   `cd evaluation`
4. **Run the evaluation**
   
   `gourd run local`

## Reproducing the evaluation on a SLURM cluster

This requires an R environment as well as [Gourd](https://github.com/ConSol-Lab/gourd/) installed on a SLURM-capable cluster.

Please note that R's random number generator may give different results depending on the software version.
The results in the paper use R version 4.3.0.


1. **Install R and Gourd**

If your cluster does not use the `module` command to enable R (not DelftBlue):
install Gourd normally using the instructions at [https://github.com/ConSol-Lab/gourd]().
You may also need to install R - consult the user manual for your HPC facility.
Then **continue to step 2**.

If your cluster uses the `module` command (DelftBlue):   please build the the Gourd version with a `modules-support` tag.
R is preinstalled.
```bash
# install the rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# clone the Gourd repository
git clone --branch modules-support --single-branch https://github.com/ConSol-Lab/gourd
cd gourd
# build
cargo build --release
# copy the program to your PATH
# on most clusters, you will need to change this to a directory in your /home
cp ./target/release/{gourd, gourd-wrapper} /usr/bin/
```

2. **Install R dependencies**

   Run the installation commands from "Using the R libraries" above in an interactive R shell. 
 
   If your cluster  uses the `module` command, first load the modules required for R. For DelftBlue:
   `module load 2024r1 r/4.3.0`. During installation, you will be asked to create a
   personal R package library: do so with the default path.
2. **Clone this repository**

   `git clone https://github.com/lchladek/RBO && cd RBO`
3. **Enter the `evaluation` directory**

   `cd evaluation`
4. **Configure the SLURM parameters**

   Edit the `gourd.toml` according to the accounting parameters of your cluster and the desired output file paths.
   See the example below with 'lchladek' as the username.

   If your cluster uses the `module` command, uncomment the 'modules' line and if necessary change the module names.
   
```
# Put all outputs and metrics in the same folder
output_path = "/scratch/lchladek/experiments"
metrics_path = "/scratch/lchladek/experiments"
experiments_folder = "/scratch/lchladek/experiments"

[slurm]
experiment_name = "RBO experimental evaluation"
output_folder = "/scratch/lchladek/slurm-output"
partition = "compute-p2"
account = "education-eemcs-courses-cse3000"
mail_type = "ALL"
#mail_user = "lchladek@tudelft.nl"
#modules = ["2024r1", "r/4.3.0" ]
```
5. **Run the evaluation**

   Use `gourd run slurm`. Check on the status using `gourd status`.
   If something does not work, refer to the  to the [Gourd documentation](https://andreats.com/gourd/).


## License

All source code is distributed under the terms of the MIT license. Different copyright applies to each subdirectory.

### `/src` and `/evaluation`

Copyright 2025 Lukáš Chládek <l@chla.cz>

The MIT license is located at `/src/LICENSE.txt` and `/evaluation/LICENSE.txt`.

### `/sigir_ap24`

This subdirectory is a modified copy of the [SIGIR-AP24](https://github.com/matteo-corsi/sigir_ap24) repository by Corsi and Urbano, 2024.
Please cite their work: [https://doi.org/10.1145/3673791.369842]()

Copyright 2024 Matteo Corsi [<corsi.matteo95@gmail.com>](mailto:corsi.matteo95@gmail.com)
and Julián Urbano [<urbano.julian@gmail.com>](mailto:urbano.julian@@gmail.com)

The MIT license is located at `/sigir_ap24/LICENSE.txt`.

