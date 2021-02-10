import numpy as np

def quantize(signal, partitions, codebook):
    indices = []
    quanta = []
    for datum in signal:
        index = 0
        while index < len(partitions) and datum > partitions[index]:
            index += 1
        indices.append(index)
        quanta.append(codebook[index])
    return np.array(indices), np.array(quanta)

def shannon_entropy_with_comments(series, L, num_int):
    """
    Function which computes the Shannon Entropy (SE) of a time series of length
    'N' using an embedding dimension 'L' and 'Num_int' uniform intervals of
    quantification. The algoritm presented by Porta et al. at "Measuring
    regularity by means of a corrected conditional entropy in sympathetic
    outflow" (PMID: 9485587) has been followed.

    Parameters:
        - series: the time series.
        - L: the embedding dimension.
        - num_int: the number of uniform intervals used in the quantification
        of the series.

    Returns:
        - SE: the SE value.
        - unique: the number of patterns which have appeared only once. This
        output is only useful for computing other more complex entropy measures
        such as Conditional Entorpy or Corrected Conditional Entropy.

    Credits:
        - Jesús Monge Álvarez. Research Master in signal theory and bioengineering - University of Valladolid
    """
    print("Original series:", series)
    series = (series - np.mean(series)) / np.std(series, ddof=1)
    print("Normalized series:", series)
    # parameters required for quantification
    epsilon = (np.max(series) - np.min(series)) / num_int
    print("epsilon:", epsilon)
    partition = np.arange(start=np.min(series), step=epsilon, stop=np.max(series)+0.1) # partition the series
    print("partition:", partition)
    codebook = np.arange(start=-1, step=1, stop=num_int+1) # which common value to assign to inputs that fall into each range of the partition
    print("codebook:", codebook)
    # uniform quantification of the time series
    indices, quants = quantize(series, partition, codebook) # slightly different
    print("Uniform quantification of the time series:")
    print("quantizations before:", quants)
    quants[quants == -1] = 0
    print("quantizations after:", quants)

    N = len(quants)
    X = quants[:]
    for j in range(1, L):
        new_row = np.append(quants[j:], np.zeros(j))
        X = np.vstack((X, new_row))

    print("Compose patterns of length 'L':")
    if X.ndim == 1:
        X = X[np.newaxis, :]
        X = X.astype(float)

    print("X:\n", X)
    print("Eliminate last 'L-1' columns of 'X' since they are not real patterns")
    X = X[:, :N - L + 1]
    print("X after\n:", X)

    num = np.ones(N-L+1); # This vector will contain the repetition of each pattern

    print("Get the number of repetitions of each pattern by going through columns of 'X':")
    for j in range(N-L+1):
        for i in range(j+1, N-L+1):
            tmp = ~np.isnan(X[:,j])
            print("col j of X:", X[:,j])
            print("col i (j+1 onwards) of X:", X[:,i])
            if np.all(tmp) and np.array_equal(X[:,j], X[:,i]):
                print("2 columns are equal, set col i to nan")
                num[j] += 1
                nan_array = np.empty((L))
                nan_array[:] = np.nan
                X[:,i] = nan_array
            tmp = np.nan

    print("Number of patterns:", num)
    # Get patterns which are not nan
    aux = ~np.isnan(X[0,:])
    # print(X[0,:])
    # print("aux:", aux)

    # print("num of patterns:", num)

    # compute number of different patterns
    new_num = num[aux]

    # number of patterns which have appeared only once
    unique = np.sum(new_num[new_num == 1])
    print("Number of patterns which have appeared only once:", unique)

    # probability of each pattern
    print("Probability of each pattern")
    p_i = new_num/(N-L+1);
    print("p_i", p_i)

    SE = np.squeeze((-1) * np.matmul((p_i)[np.newaxis,:], np.transpose(np.log(p_i)[np.newaxis])))
    print("Compute Shannon Entropy:", SE)

    return SE, unique

def shannon_entropy(series, L, num_int):
    # normalize input time series
    # print(np.mean(series), np.std(series))
    series = (series - np.mean(series)) / np.std(series, ddof=1)
    # parameters required for quantification
    epsilon = (np.max(series) - np.min(series)) / num_int
    partition = np.arange(start=np.min(series), step=epsilon, stop=np.max(series)+0.1) # partition the series
    codebook = np.arange(start=-1, step=1, stop=num_int+1) # which common value to assign to inputs that fall into each range of the partition
    # uniform quantification of the time series
    indices, quants = quantize(series, partition, codebook) # slightly different
    quants[quants == -1] = 0

    N = len(quants)
    X = quants[:]
    for j in range(1, L):
        new_row = np.append(quants[j:], np.zeros(j))
        X = np.vstack((X, new_row))

    if X.ndim == 1:
        X = X[np.newaxis, :]
        X = X.astype(float)

    X = X[:, :N - L + 1]

    num = np.ones(N-L+1); # This vector will contain the repetition of each pattern

    for j in range(N-L+1):
        for i in range(j+1, N-L+1):
            tmp = ~np.isnan(X[:,j])
            if np.all(tmp) and np.array_equal(X[:,j], X[:,i]):
                num[j] += 1
                nan_array = np.empty((L))
                nan_array[:] = np.nan
                X[:,i] = nan_array
            tmp = np.nan

    # Get patterns which are not nan
    aux = ~np.isnan(X[0,:])

    # compute number of different patterns
    new_num = num[aux]

    # number of patterns which have appeared only once
    unique = np.sum(new_num[new_num == 1])

    # probability of each pattern
    p_i = new_num/(N-L+1);

    SE = np.squeeze((-1) * np.matmul((p_i)[np.newaxis,:], np.transpose(np.log(p_i)[np.newaxis])))

    return SE, unique

def conditional_entropy_with_comments(series, L, num_int):
    """
        Function which computes the Conditional Entropy (CE) of a time series of
        length 'N' using an embedding dimension 'L' and 'Num_int' uniform intervals
        of quantification. The algoritm presented by Porta et al. at "Measuring
        regularity by means of a corrected conditional entropy in sympathetic
        outflow" (PMID: 9485587) has been followed.

        Parameters:
            - series: the time series.
            - L: the embedding dimension.
            - num_int: the number of uniform intervals used in the quantification
            of the series.

        Returns:
            - CE: the CE value.
            - unique: the number of patterns which have appeared only once. This
            output is only useful for computing other more complex entropy
            measures such as Corrected Conditional Entorpy.

        Credits:
            - Jesús Monge Álvarez. Research Master in signal theory and bioengineering - University of Valladolid
    """
    # L as embedding dimension
    SE, unique = shannon_entropy(series, L, num_int)
    print("Calculate Shanon entropy for series with embedding dimension 'L':", SE, ", number of unique patterns:", unique)

    # L - 1 as embedding dimension
    SE_l, unique_l = shannon_entropy(series, L-1, num_int)
    print("Calculate Shanon entropy for series with embedding dimension 'L-1':", SE_l, ", number of unique patterns:", unique)

    print("Conditional entropy:", SE - SE_l)
    return SE - SE_l, unique

def conditional_entropy(series, L, num_int):
    """
        Function which computes the Conditional Entropy (CE) of a time series of
        length 'N' using an embedding dimension 'L' and 'Num_int' uniform intervals
        of quantification. The algoritm presented by Porta et al. at "Measuring
        regularity by means of a corrected conditional entropy in sympathetic
        outflow" (PMID: 9485587) has been followed.

        Parameters:
            - series: the time series.
            - L: the embedding dimension.
            - num_int: the number of uniform intervals used in the quantification
            of the series.

        Returns:
            - CE: the CE value.
            - unique: the number of patterns which have appeared only once. This
            output is only useful for computing other more complex entropy
            measures such as Corrected Conditional Entorpy.

        Credits:
            - Jesús Monge Álvarez. Research Master in signal theory and bioengineering - University of Valladolid
    """
    # L as embedding dimension
    SE, unique = shannon_entropy(series, L, num_int)

    # L - 1 as embedding dimension
    SE_l, unique_l = shannon_entropy(series, L-1, num_int)

    return SE - SE_l, unique

def corrected_conditional_entropy_with_comments(series, Lmax, num_int):
    """
        Function which computes the Corrected Conditional Entropy (CCE) of a time series of
        length 'N' using an embedding dimension 'L' and 'Num_int' uniform intervals
        of quantification. The algoritm presented by Porta et al. at "Measuring
        regularity by means of a corrected conditional entropy in sympathetic
        outflow" (PMID: 9485587) has been followed.

        Parameters:
            - series: the time series.
            - L: the embedding dimension.
            - num_int: the number of uniform intervals used in the quantification
            of the series.

        Returns:
            - CCE_min: the CCE value. The best estimation of the CCE is the
            minimum value of all the CCE that have been computed.

        Credits:
            - Jesús Monge Álvarez. Research Master in signal theory and bioengineering - University of Valladolid
    """
    N = len(series)

    E_est_1, unique = shannon_entropy(series, 1, num_int)
    print("Calculate Ê(1) with 'L=1' which will be used to calculate the corrective term:", E_est_1)

    CCE = np.empty((Lmax))
    CCE[:] = np.nan

    CCE[0] = 100
    # print("CCE:", CCE)
    CE = np.empty((Lmax))
    CE[:] = np.nan

    uniques = np.empty((Lmax))
    uniques[:] = np.nan

    correc_term = np.empty((Lmax))
    correc_term[:] = np.nan

    # print("original values:", CCE, CE, uniques, correc_term)
    print("CCE will be a vector that will contain several CCE values computed:")
    print("We loop for different value of embedding dimensions 'L':")
    for L in range(2, Lmax + 1):
        print("\nL:", L, "\n")
        print("First, compute CE for the current embedding dimension:")
        CE[L-1], uniques[L-1] = conditional_entropy(series, L, num_int)
        print("CE:", CE)
        print("uniques:", uniques)

        print("Second, compute the percentage of patterns which are not repeated")
        perc_L = uniques[L-1]/(N-L+1)
        print("perc_L:", perc_L)
        print("Ê(1):", E_est_1)
        correc_term[L-1] = perc_L * E_est_1

        print("CCE is CE + corrective term:")
        print("correct_term:", correc_term)
        CCE[L-1] = CE[L-1] + correc_term[L-1]
        print("CCE:", CCE)

    print("\nFinal CCE:", CCE)
    print("Get min CCE value:", np.min(CCE))
    return np.min(CCE)

def corrected_conditional_entropy(series, Lmax, num_int):
    """
            Function which computes the Corrected Conditional Entropy (CCE) of a time series of
            length 'N' using an embedding dimension 'L' and 'Num_int' uniform intervals
            of quantification. The algoritm presented by Porta et al. at "Measuring
            regularity by means of a corrected conditional entropy in sympathetic
            outflow" (PMID: 9485587) has been followed.

            Parameters:
                - series: the time series.
                - L: the embedding dimension.
                - num_int: the number of uniform intervals used in the quantification
                of the series.

            Returns:
                - CCE_min: the CCE value. The best estimation of the CCE is the
                minimum value of all the CCE that have been computed.

            Credits:
                - Jesús Monge Álvarez. Research Master in signal theory and bioengineering - University of Valladolid
        """
    N = len(series)

    E_est_1, unique = shannon_entropy(series, 1, num_int)

    CCE = np.empty((Lmax))
    CCE[:] = np.nan

    CCE[0] = 100
    CE = np.empty((Lmax))
    CE[:] = np.nan

    uniques = np.empty((Lmax))
    uniques[:] = np.nan

    correc_term = np.empty((Lmax))
    correc_term[:] = np.nan

    for L in range(2, Lmax + 1):
        CE[L - 1], uniques[L - 1] = conditional_entropy(series, L, num_int)

        perc_L = uniques[L - 1] / (N - L + 1)
        correc_term[L - 1] = perc_L * E_est_1

        CCE[L - 1] = CE[L - 1] + correc_term[L - 1]

    return np.min(CCE)