# Entropy Measures for Analysis of Biomedical Signals

[Entropy Measures for Analysis of Biomedical Signals paper](https://www.researchgate.net/publication/290219771_Entropy-Based_Algorithms_in_the_Analysis_of_Biomedical_Signals)

![](https://i.imgur.com/xMF2y5B.png)

Other References:
[Exploring Entropy Measurements to Identify Multi-Occupancy in Activities of Daily Living](http://web.a.ebscohost.com.libproxy1.nus.edu.sg/plink?key=10.83.8.64_8000_1334801804&site=ehost&scope=site&db=a9h&AN=136174470&msid=-419426216)

[Characterization of Surface EMG Signal Based on Fuzzy Entropy
](https://ieeexplore-ieee-org.libproxy1.nus.edu.sg/document/4237165)

[Permutation Entropy: A Natural Complexity Measure for Time Series](https://journals-aps-org.libproxy1.nus.edu.sg/prl/abstract/10.1103/PhysRevLett.88.174102)

[Multiscale Entropy Analysis: A New Method to Detect Determinism in a Time Series](https://arxiv.org/pdf/physics/0604040.pdf)

[Dispersion Entropy: A Measure for Time-Series Analysis](https://ieeexplore-ieee-org.libproxy1.nus.edu.sg/document/7434608)

## What is Entropy and Why is it used

Entropy has become an appropriate measure to study time series from biological systems, which are characterised by complex dynamics. 

Entropy has emerged as a suitable complexity measure for the amount of disorder or uncertainty in a system or time-series data.

### Shannon Entropy 

In communication theory, Shannon et al. (1949) introduced Shannon Entropy. Information entropy S of a random variable X that takes values $x_1, x_2,...,x_N$ from a set of values $\theta$ and probability mass function $p(x_i) = P_r\{X=x_i\}, x_i \in \theta$ is defined as:

$$S_{en} = \sum^{n}_{i=1} p(x_i) log_a{\frac{1}{p(x_i)}} = \sum^{n}_{i=1} p(x_i) log_a{p(x_i)} , a >1$$

Shannon entropy denotes the degree of uncertainty associated with the occurrence of the result. A higher value of entropy gives a more uncertain outcome and is more difficult to predict.

### Approximate Entropy 
[Approximate entropy as a measure of system complexity - 1991
](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC51218/?tool=pmcentrez&report=abstract)

#### Rationale

Numerous entropy algorithms have been introduced to quantify the irregularity of signals. The computations, however, are limited with these challenges:
- Insufficient data points
- Recorded data contaminated with noise

Since real-world recorded physiological signals are often short and noisy, Approximate Entropy (ApEn) was proposed to avert the challenges in finite length of a time series. 

- High regularity and low randomness -> smaller entropy values
- Less regularity -> higher entropy values

#### Formula

- Original work: Kolmogorov-Sinai (KS) entropy measures that mean rate of creation of information, or in other words, the decrease in uncertainty at a receiver by knowing the current state of the system given the past history. However, this formula is of limited use to estimate entropy of finite length "real-world" time series.
- Grass- berger and Procaccia proposed a fomula to estimate KS entropy with reasonable precision, characterizing chaotic signals by calculating $K_2$ entropy which is a lower bound of KS entropy.

- Let {X_i} = {x_1, ..., x_i, ..., x_N} represents a time series of length N. Consider $m$-length vectors: $u_m(i) = \{x_i, x_{i+1}, ... , x_{i+m-1}\}, 1 \leq i \leq N-m+1$. Let $n^m_i(r)$ represents the number of vectors $u_m(j)$ that are close to vector $u_m(i)$, which means the number of vectors that satisfy $d[u_m(i), u_m(j)] \leq r$, where $d$ is the Euclidean distance. $C^m_i(r) = n^m_i(r)/(N-m+1)$ represents the probability that any vector $u_m(j)$ is close to the vector $u_m(i)$. The average of $C^m_i$, 
$$C^m(r) = 1/(N-m+1) \sum^{(N-m+1)}_{i=1}C^m_i(r)$$ represents the probability that any two vectors are within r of each other. $K_2$ is defined as: $$K_2 = \lim_{N \to \infty} \lim_{m \to \infty} \lim_{r \to 0} -ln[C^{m+1}(r) - C^m(r)]$$

- Eckmann and Ruelle further defined the function $\phi^m(r) = \frac{1}{N-m+1} \sum^{N-m+1}_{i=1} ln C^m_i(r)$, considering the distance between two vectors as maximum absolute difference between their components: $d[u^m(i), u^m(j)] = max\{[x(i+k)-x(j+k)]:0 \leq k \leq m-1 \}$. As such, $\phi^{m+1}(r) - \phi^m(r) \approx \sum^{N-m+1}_{i=1} ln[\frac{C^m_i(r)}{C^{m+1}_i(r)}]$ represents the **average of the natural logarithm of the conditional probability that sequences that are close to each other for $m$ consecutive data points will be close to each other when one more data point is known**. Therefore, KS entropy is estimated as: $$H_{ER}= \lim_{N \to \infty} \lim_{m \to \infty} \lim_{r \to 0} [\phi^{m}(r) - \phi^{m+1}(r)]$$

- Finally, Pincus et al. introduced a family of measures termed approximate entropy, $A_E(m,r)$ defined as: $$A_E(m,r) = \lim_{N \to \infty}[\phi^{m}(r) - \phi^{m+1}(r)]$$ and $A_E$ is estimated by the statistics $$A_E(m,r,N) = \phi^{m}(r) - \phi^{m+1}(r)$$

#### Usage
$A_E$ is not intended as an approximate value of ER entropy. It is instead a **regularity statistics**. It can be applied to "real-world" time series.

- Lower $A_E$ values are assigned to more regular time series while higher $A_E$ are assigned to more irregular, less predictable time series.

#### Choice of paramerters $m , r , N$
Usually, $r$ = 20% of the standard deviation of the amplitude values and $m = 2$

#### Disadvantages

- Strongly dependent on the record length
- Is often lower than expected for short records
- Lacks relative consistency 

#### Python Implementation
[EntroPy](https://github.com/raphaelvallat/entropy)

### Sample Entropy

[Physiological time-series analysis using approximate entropy and sample entropy - 2000](https://journals-physiology-org.libproxy1.nus.edu.sg/doi/full/10.1152/ajpheart.2000.278.6.H2039)

#### Rationale

Highlight the limitations of Approximate Entropy: ApEn includes self-matches (count each sequence as matching itself) to avoid the occurrence of $ln(0)$ which leads to **bias** that causes lack of two important properties:
1. Heavily dependent on record length and Is uniformly lower than expected for short records
2. Lacks relative consistency: if ApEn of one data set is higher than that of another, it should but does not remain higher for all conditions tested

- To reduce the bias, introduce a new family of statistics, Sample Entropy (SampEn) that does not count self-matches.

    - $SampEn(m, r, N)$ is also the negative natural logarithm of the conditional probability that two sequences similar for $m$ points remain similar at the next point, where self-matches are not included in calculating probability. Thus, lower SampEn value also indicates more self-similarity in the series and SampEn only requires about half time to calculate.
    
#### Properties

SampEn is:
- Largely independent of record length
- Display relative consistency under circumstances where ApEn does not, over a broader range of possible $r,m,N$ values
    
#### Formula

We also need to define $\phi^{m}(r)$ and $\phi^{m+1}(r)$. The probability $\phi^{m}(r)$ that two sequences match for $m$ points is computed by counting the average number of vector pairs for which the distance is lower than the tolerance $r$: $$\phi^{m}(r) = \frac{1}{N-m} \sum^{N-m}_{i=1}C^m_i(r)$$ Similary, $\phi^{m+1}(r)$ is defined for an embedding dimension $m+1$. $SampEn$ is calculated as: $$Samp_{en} = ln\frac{\phi^{m}(r)}{\phi^{m+1}(r)}$$

More precisely, starting from definition of $K_2$ entropy, $$S_E(m,r) = \lim_{N \to \infty} -ln\frac{U^{m+1}(r)}{U^m(r)}$$ which is estimated by statistic: $$S_E(m,r,N) = -ln\frac{U^{m+1}(r)}{U^m(r)}$$

$U^{m}(r)$ is different from $C^{m}(r)$ in these:
- Distance between 2 vectors are maximum absolute difference between their components
- Exclude self-matches (vectors not compared to themselves)
- Given a time series with N points, only the first $N-m$ vectors of length $m$, $u_m(i)$ are considered so that for $1 \leq i \leq N-m$, the vector $u_{m+1}(i)$ of length $m+1$ is also defined. 

#### Limitations
Both $S_E$ and $A_E$ measure the degree of randomness of a time series. However, there is no straightforward relationship between regularity, measured by entropy-based metrics, and complexity. 

**An increase in entropy is usually not always associated with an increase in complexity**

Entropy-based metrics are maximised for random sequences, although both **perfectly ordered** and **maximally disordered** systems are agreed to possess no complex structures. So a more meaningful physiologic complexity measure should **vanish** for these two extreme states.

Hence, Costa M., Goldberger A. L., Peng C. -K. introduced **Multiscale entropy**.

#### Python implementation
[EntroPy](https://github.com/raphaelvallat/entropy)

### Corrected Conditional Entropy

[Measuring regularity by means of a corrected conditional entropy in sympathetic outflow - 1997](https://link-springer-com.libproxy1.nus.edu.sg/article/10.1007%2Fs004220050414)

#### Rationale

ApEn is limited: Whatever process (even white noise) is considered, a short data sequence generates zero values of entropy rate for a large (w.r.t total amount of data) length of pattern $(m)$. Hence, we cannot quantify the regularity of the time series. 

In other words, embedding the dynamics into a too highly dimensional phase space produces entropy rate estimates equal to 0. This forces Pincus et al. (1993) to fix the pattern lengh (i.e. the embedding dimension m) at very small arbitraty values in order to obtain a reliable entropy rate estimate.

- Aim: Propose a new method for quantification of regularity of a process from short data sequences, based on defining a new CCE function and search for its minimum in respect of the pattern length. 

- Attempt to **distinguish** the **decrease of entropy rate** related to **presence of recurrences** from that related to the **shortness of the data sequence**

#### Formula

1. Conditional Entropy

To take advantage of stationary, normalize the series $\{X(i)\}$ to a process with 0 mean and unitary variance $$x(i) = \frac{X(i) - av[X]}{std[X]}$$

From the normalized series, a reconstructed $L-dimensional$ phase space with a delay of reconstruction equal to 1 is obtained by considering $N-L+1$ vectors $x_L(i) = (x(i), x(i-1), ..., x(i-L+1)$. Each vector $x_L(i)$ represents a pattern of L consecutive samples.

CE is defined: $$E(L/L-1) = - \sum_{L-1} p_{L-1} \sum_{L/L-1} p_{L/L-1} \log p_{L/L-1}$$

where $p_{L-1}$ denotes the joint probability of the pattern $x_{L-1}(i)$ and $p_{L/L-1}$ symbolizes conditional probability of the Lth sample of the pattern $x_L(i)$ given the previous $L-1$ ones.

CE can be directly derived from definition of Shannon Entropy (SE) of $x_L(i)$ $$E(L) = - \sum_L p_L \log p_L$$

CE can be obtained as the variation of the SE w.r.t L: $$E(L/L-1) = E(L) - E(L - 1)$$

- SE represents the amount of information needed to specify the point $x_L(i)$ in a L-dimensional phase space.

- CE quantifies the variation of information necessary to specify a new state in a one-dimension incremental phase space. Small CE values are obtained when a length L pattern can be almost completely predicted by a length L - 1 pattern.

2. Conditional entropy estimate

$$E(L/L-1) = E(L) - E(L-1)$$ where $E(L)$ and $E(L-1)$ represent estimate of SE in a L-dimensional and (L-1)-dimensional phase space. $E(L)$ can be estimated by approximating probabilities in $E(L) = - \sum_L p_L \log p_L$ 