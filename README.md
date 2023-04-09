# modern_stat_4_modern_biology
A worksheet for Modern Statistics for Modern Biology book

## 1  Generative Models for Discrete Data

### 1.2 A real example

```R
dpois(x = 3, lambda = 5)
#[1] 0.1403739

#to generate the probabilities of all values from 0 to 12,
dpois(x = 0:12, lambda = 5)

#[1] 0.006737947 0.033689735 0.084224337 0.140373896 0.175467370 0.175467370 0.146222808 0.104444863
#[9] 0.065278039 0.036265577 0.018132789 0.008242177 0.003434240

# Plot
barplot(dpois(0:12, 5), names.arg = 0:12, col = "red", main = "Probabilities of seeing 0,1,2,…,12 mutations, as modeled by the Poisson(5) distribution")
```
![Fig_1](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_1.png) 

### 1.3 Using discrete probability models

```R
genotype = c("AA","AO","BB","AO","OO","AO","AA","BO","BO",
             "AO","BB","AO","BO","AB","OO","AB","BB","AO","AO")
table(genotype)
#genotype
#AA AB AO BB BO OO 
# 2  2  7  3  3  2 
genotypeF = factor(genotype)
levels(genotypeF)
#[1] "AA" "AB" "AO" "BB" "BO" "OO"
```
### 1.3.1 Bernoulli trials
Suppose we want to simulate a sequence of 15 fair coin tosses. To get the outcome of 15 Bernoulli trials with a probability of success equal to 0.5 (a fair coin), we write
```R
rbinom(15, prob = 0.5, size = 1)
#[1] 0 0 1 1 1 0 1 0 1 0 1 1 0 0 1
```
To simulate twelve trials of throwing a ball into the two boxes  with probability of falling in the right-hand box 2/3 in the left-hand box 1/3, we write:

```R
rbinom(12, prob = 2/3, size = 1)
#[1] 1 1 1 1 1 1 1 0 1 1 1 1
```
### 1.3.2 Binomial success counts

```R
rbinom(1, prob = 2/3, size = 12)
```
The number of successes in 15 Bernoulli trials with a probability of success of 0.3 is called a binomial random variable or a random variable that follows the 
 distribution.
 ```R
 set.seed(235569515)
rbinom(1, prob = 0.3, size = 15)
```
The complete **probability mass distribution** is available by typing:
```R
probabilities = dbinom(0:15, prob = 0.3, size = 15)
round(probabilities, 2)
# [1] 0.00 0.03 0.09 0.17 0.22 0.21 0.15 0.08 0.03 0.01 0.00 0.00 0.00 0.00 0.00 0.00

barplot(probabilities, names.arg = 0:15, col = "red", main = "Theoretical distribution of B(15, 0.3)")
```
![Fig_2](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_1_2.png)

### 1.3.3 Poisson distributions
When the probability of success *p* is small and the number of trials *n* large, the binomial distribution *B(p,n)* can be faithfully approximated by a simpler distribution, the Poisson distribution with rate parameter *&lambda; = np*
We have seen this in HIV example:
- mutation rate 5*10^-4 per nts per replication cycle
- HIV genome size = 10^4
```R
rbinom(1, prob = 5e-4, size = 10000)
```
Simulate a mutation process along 10,000 positions with a mutation rate of 5*10^-4 and count the number of mutations.Repeat this many times and plot the distribution
```R
simulations = rbinom(n = 300000, prob = 5e-4, size = 10000)
barplot(table(simulations), col = "lavender", main = 'Simulated distribution of B(10000, 10^-4) for 300000 simulations')
```
![Fig_3](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_1_3.png)

### 1.3.4 A generative model for epitope detection

High level story: Detection of allergic reaction with ELISA for 50 patients and testing for 100 epitope with an FDR 1%

We’re going to study the data for all 50 patients tallied at each of the 100 positions. If there are no allergic reactions, the false positive rate means that for one patient, each individual position has a probability of 1 in 100 of being a 1. So, after tallying 50 patients, we expect at any given position the sum of the 50 observed (0,1) variables to have a Poisson distribution with parameter 0.5. 
```R
load("data/e100.RData")

barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
  names.arg = seq(along = e100), col = "darkolivegreen", main = "Output of the ELISA array results for 50 patients in the 100 positions")
```
![Fig_4](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_1_4.png)

The spike in the above figure is striking. What are the chances of seeing a value as large as 7, if no epitope is present?
If we look for the probability of seeing a number as big as 7 (or larger) when considering one Poisson(0.5) random variable, the answer can be calculated in closed form as:

<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">
  <mi>P</mi>
  <mo stretchy="false">(</mo>
  <mi>X</mi>
  <mo>&#x2265;</mo>
  <mn>7</mn>
  <mo stretchy="false">)</mo>
  <mo>=</mo>
  <munderover>
    <mo data-mjx-texclass="OP">&#x2211;</mo>
    <mrow data-mjx-texclass="ORD">
      <mi>k</mi>
      <mo>=</mo>
      <mn>7</mn>
    </mrow>
    <mi mathvariant="normal">&#x221E;</mi>
  </munderover>
  <mi>P</mi>
  <mo stretchy="false">(</mo>
  <mi>X</mi>
  <mo>=</mo>
  <mi>k</mi>
  <mo stretchy="false">)</mo>
  <mo>.</mo>
</math>

This is, of course, the same as <math xmlns="http://www.w3.org/1998/Math/MathML">
  <mn>1</mn>
  <mo>&#x2212;</mo>
  <mi>P</mi>
  <mo stretchy="false">(</mo>
  <mi>X</mi>
  <mo>&#x2264;</mo>
  <mn>6</mn>
  <mo stretchy="false">)</mo>
</math>
The probability  P(X <= 6) is the so-called **cumulative distribution function** at 6, and R has the function ```ppois``` for computing it, which we can use in either of the following two ways:
```R
1 - ppois(6, 0.5)
#[1] 1.00238e-06
# OR 
ppois(6, 0.5, lower.tail = FALSE)
#[1] 1.00238e-06
```
So <math xmlns="http://www.w3.org/1998/Math/MathML" display="block">
  <mi>&#x3F5;</mi>
  <mo>=</mo>
  <mi>P</mi>
  <mo stretchy="false">(</mo>
  <mi>X</mi>
  <mo>&#x2265;</mo>
  <mn>7</mn>
  <mo stretchy="false">)</mo>
  <mo>=</mo>
  <mn>1</mn>
  <mo>&#x2212;</mo>
  <mi>P</mi>
  <mo stretchy="false">(</mo>
  <mi>X</mi>
  <mo>&#x2264;</mo>
  <mn>6</mn>
  <mo stretchy="false">)</mo>
  <mo>&#x2243;</mo>
  <msup>
    <mn>10</mn>
    <mrow data-mjx-texclass="ORD">
      <mo>&#x2212;</mo>
      <mn>6</mn>
    </mrow>
  </msup>
  <mo>.</mo>
</math>

The above calculation is not the correct computation in this case.We looked at all 100 positions, looked for the largest value and found that it was 7. Due to this selection, a value as large as 7 is more likely to occur than if we only looked at one position.
So instead of asking what the chances are of seeing a Poisson(0.5) as large as 7, we should ask ourselves, what are the chances that the maximum of 100 Poisson(0.5) trials is as large as 7? If we follow the extereme value analysis (code and calculation not shown here), the probability would be 10^-4 (vs 10^-6 above).

