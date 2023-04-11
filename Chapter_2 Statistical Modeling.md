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
