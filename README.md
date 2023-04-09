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
[Fig_1](https://github.com/hamidghaedi/modern_stat_4_modern_biology/blob/main/figs/Fig_1.png)