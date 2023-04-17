## 2  Statistical Modeling

statistical modeling starts with the available data and aims to build a model that can explain the observed patterns in the data. This process involves making assumptions about the underlying processes and relationships between variables, and then fitting the model to the data to estimate the unknown parameters. The estimated parameters can then be used to make predictions or draw conclusions about the population of interest.

However, in many cases, the true model and parameters are unknown, and researchers need to use the available data to estimate them. This is where statistical inference comes in, which involves using the available data to make inferences about the unknown parameters or the model itself.

### 2.3 A simple example of statistical modeling

```R
load("data/e100.RData")
# exclude tricky outlier 
e99 = e100[-which.max(e100)]
```
#### Goodness-of-fit : visual evaluation
```R
barplot(table(e99), space = 0.8, col = "chartreuse4")
```
However, it is hard to decide which theoretical distribution fits the data best without using a comparison. One visual goodness-of-fit diagram is known as the **rootogram**; it hangs the bars with the observed counts from the theoretical red points. If the counts correspond exactly to their theoretical values, the bottom of the boxes will align exactly with the horizontal axis.

```R
library("vcd")
gf1 = goodfit( e99, "poisson")
rootogram(gf1, xlab = "", rect_gp = gpar(fill = "chartreuse4"), main = "rootgram for goodness-of-fit visual evaluation")
```
![Fig2_1](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_2_1.png)

Rootogram showing the square root of the theoretical values as red dots and the square root of the observed frequencies as drop down rectangles. We see that the rootogram for e99 seems to fit the Poisson model reasonably well. But remember, to make this happen we removed the outlier. The Poisson is completely determined by one parameter, often called the Poisson mean &lambda;. In most cases where we can guess the data follows a Poisson distribution, we will need to estimate the Poisson parameter from the data.

The most common way of estimating &lambda; , is to choose the value &#955;&#770; that makes the observed data the most likely. This is called the **maximum likelihood estimator (MLE)**. 

The shortcut to calculate Poisson parameter for ```e100``` data is :
```R
gf  =  goodfit(e100, "poisson")
names(gf)
#[1] "observed" "count"    "fitted"   "type"     "method"   "df"       "par"     

gf$par
#$lambda
#[1] 0.55
```
### Binomial distributions and maximum likelihood

Suppose we take a sample of n=120 males and test them for red-green colorblindness. We can code the data as 0 if the subject is not colorblind and 1 if he is. We summarize the data by the table:
 
```R
cb <- c(rep(0, 110), rep(1, 10))

table(cb)
#cb
#  0   1 
#110  10 

mean(cb)
#[1] 0.08333333
```
In the case of the Poisson, if we compute the likelihood for many possible ```p``` , we can plot it and see where its maximum falls:

```R
probs  =  seq(0, 0.3, by = 0.005)
likelihood = dbinom(sum(cb), prob = probs, size = length(cb))
plot(probs, likelihood, pch = 16, xlab = "probability of success",
       ylab = "likelihood", cex=0.6)
probs[which.max(likelihood)]
#[1] 0.085
```
![Fig2_2](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_2_2.png)

##### Likelihood for the binomial distribution
Suppose ```n=300```, and we observe ```y=40```successes. Then, for the binomial distribution:
```R
loglikelihood = function(p, n = 300, y = 40) {
  log(choose(n, y)) + y * log(p) + (n - y) * log(1 - p)
}
p_seq = seq(0, 1, by = 0.001)
plot(p_seq, loglikelihood(p_seq), xlab = "p", ylab = "log f(p|y)", type = "l", main = "Plot of the log likelihood function for n=300, y=40")
```
![Fig2_3](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_2_3.png)

#### multinomial data
##### DNA count modeling: base pairs
This section combines estimation and testing by simulation in a real example. Data from one strand of DNA for the genes of Staphylococcus aureus bacterium are available in a fasta file ```staphsequence.ffn.txt```, which we can read with a function from the Bioconductor package:
```R
library("Biostrings")
staph = readDNAStringSet("data/staphsequence.ffn.txt", "fasta")
#Let’s look at the first gene:
staph[1]
#DNAStringSet object of length 1:
#    width seq                                                                       names            #   
#[1]  1362 ATGTCGGAAAAAGAAATTTGGGAAAAAGTGCTTGA...TAGAGAATCTTGAAAAAGAAATAAGAAATGTATAA #lcl|NC_002952.2_c...

letterFrequency(staph[[1]], letters = "ACGT", OR = 0)
#  A   C   G   T 
#522 219 229 392 
```
Due to their different physical properties, evolutionary selection can act on the nucleotide frequencies. So we can ask whether, say, the first ten genes from these data come from the same multinomial.
```R
# calculate letter freq for genes in staph
letterFrq = vapply(staph, letterFrequency, FUN.VALUE = numeric(4),
         letters = "ACGT", OR = 0)
# set column names
colnames(letterFrq) = paste0("gene", seq(along = staph))
tab10 = letterFrq[, 1:10]
computeProportions = function(x) { x/sum(x) }
prop10 = apply(tab10, 2, computeProportions)
round(prop10, digits = 2)

#  gene1 gene2 gene3 gene4 gene5 gene6 gene7 gene8 gene9 gene10
#A  0.38  0.36  0.35  0.37  0.35  0.33  0.33  0.34  0.38   0.27
#C  0.16  0.16  0.13  0.15  0.15  0.15  0.16  0.16  0.14   0.16
#G  0.17  0.17  0.23  0.19  0.22  0.22  0.20  0.21  0.20   0.20
#T  0.29  0.31  0.30  0.29  0.27  0.30  0.30  0.29  0.28   0.36

p0 = rowMeans(prop10)
#p0
#        A         C         G         T 
#0.3470531 0.1518313 0.2011442 0.2999714 
```
So let’s suppose `p0` is the vector of **multinomial probabilities** for all the ten genes and use a *Monte Carlo simulation* to test whether the departures between the observed letter frequencies and expected values under this supposition are within a plausible range.

We compute the expected counts by taking *the outer product* of the vector of probabilities p0 with the sums of nucleotide counts from each of the 10 columns, `cs`.
```R
cs = colSums(tab10)
#cs
# gene1  gene2  gene3  gene4  gene5  gene6  gene7  gene8  gene9 gene10 
#  1362   1134    246   1113   1932   2661    831   1515   1287    696 

expectedtab10 = outer(p0, cs, FUN = "*")
round(expectedtab10)

#  gene1 gene2 gene3 gene4 gene5 gene6 gene7 gene8 gene9 gene10
#A   473   394    85   386   671   924   288   526   447    242
#C   207   172    37   169   293   404   126   230   195    106
#G   274   228    49   224   389   535   167   305   259    140
#T   409   340    74   334   580   798   249   454   386    209
```
We can now create a random table with the correct column sums using the `rmultinom` function. This table is generated according to the null hypothesis that the true proportions are given by p0.
```R
randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) } )
all(colSums(randomtab10) == cs)

#[1] TRUE

```
Now we repeat this `B = 1000` times. For each table we compute our test statistic from the function `stat` and store the results in the vector `simulstat`. Together, these values constitute our null distribution, as they were generated under the null hypothesis that p0 is the vector of multinomial proportions for each of the 10 genes.

```R
stat = function(obsvd, exptd) {
   sum((obsvd - exptd)^2 / exptd)
}
B = 1000
simulstat = replicate(B, {
  randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) })
  stat(randomtab10, expectedtab10)
})
S1 = stat(tab10, expectedtab10)
sum(simulstat >= S1)

#[1] 0
hist(simulstat, col = "lavender", breaks = seq(0, 75, length.out=50))
abline(v = S1, col = "red")
abline(v = quantile(simulstat, probs = c(0.95, 0.99)),
       col = c("darkgreen", "blue"), lty = 2)
```
![Fig2_4](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_2_4.png)

The above is a histogram of simulstat. The value of S1 is marked by the vertical red line, those of the 0.95 and 0.99 quantiles  by the dotted lines.
We see that the probability of seeing a value as large as S1=70.1 is very small under the null model. It happened 0 times in our 1000 simulations that a value as big as S1 occurred. Thus the ten genes do not seem to come from the same multinomial model.



