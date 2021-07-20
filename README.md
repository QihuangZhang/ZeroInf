# Zero-Inflated Poisson Models with Measurement Error in Response

2021-07-20

Data and code for reproducing the results in the article of "Zero-Inflated Poisson Models with Measurement Error in Response"
by [Qihuang Zhang](https://qihuangzhang.com) and [Grace Y. Yi](https://www.uwo.ca/stats/people/bios/Yi,%20Grace.html).






Before implementing the code, R package [GeneErrorMis](https://github.com/QihuangZhang/GeneErrorMis) is needed to be installed:

``` r
# install devtools if necessary
install.packages('devtools')

# install the GeneErrorMis package
devtools::install_github('QihuangZhang/GeneErrorMis')

# load
library(GeneErrorMis)
```


## File Structure
* [Naive Method](https://github.com/QihuangZhang/ZeroInf/blob/master/code/Simulation/Simulation0.R)
* [Section 6.1： Simulation 1](https://github.com/QihuangZhang/ZeroInf/blob/master/code/Simulation/Simulation1.R)
* [Section 6.2： Simulation 2](https://github.com/QihuangZhang/ZeroInf/blob/master/code//Simulation/Simulation2.R)
* Section 6.3： Simulation 3
  * [Proposed Method](https://github.com/QihuangZhang/ZeroInf/blob/master/code/Simulation/Simulation3.R)
  * [External Validation](https://github.com/QihuangZhang/ZeroInf/blob/master/code/Simulation/Simulation4.R): X_i is not available but W_i is available in validation data
  * [Internal Validation](https://github.com/QihuangZhang/ZeroInf/blob/master/code/Simulation/Simulation5.R): Both X_i and W_i is available in validation data
