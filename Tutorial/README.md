# Tutorial of "Zero-Inflated Poisson Models with Measurement Error in Response"

2021-07-20

Tutorial of "Zero-Inflated Poisson Models with Measurement Error in Response"
by [Qihuang Zhang](https://qihuangzhang.com) and [Grace Y. Yi](https://www.uwo.ca/stats/people/bios/Yi,%20Grace.html).







## Steps of Implementing Example
* Step 1: Install the R package

Before implementing the code, R package [GeneErrorMis](https://github.com/QihuangZhang/GeneErrorMis) is needed to be installed:

``` r
# install devtools if necessary
install.packages('devtools')

# install the GeneErrorMis package
devtools::install_github('QihuangZhang/GeneErrorMis')

# load
library(GeneErrorMis)
```

* Step 2: Configure the global parameters
Modify the following global parameter in "main.R"

seed: the seed of generating simulation data

* Step 3: Implement method in "main.R"



## File Structure
* [Main](https://github.com/QihuangZhang/ZeroInf/blob/master/Tutorial/Main.R): Main R file in implementating the tutorial exmaple
* [Auxiliary Functions](https://github.com/QihuangZhang/ZeroInf/blob/master/Tutorial/Functions.R): Auxiliary Functions that are needed during the implementation
* [Example Resulting Data](https://github.com/QihuangZhang/ZeroInf/blob/master/Tutorial/ExampleResults.RData): An example of the generated results
