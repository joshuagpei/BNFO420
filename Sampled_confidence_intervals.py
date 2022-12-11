from math import log2
import numpy as np
import matplotlib.pyplot as plt
import time
from pprint import pprint

def multidimensional_shifting(num_samples, sample_size, elements, probabilities):
    # replicate probabilities as many times as `num_samples`
    replicated_probabilities = np.tile(probabilities, (num_samples, 1))    # get random shifting numbers & scale them correctly
    random_shifts = np.random.random(replicated_probabilities.shape)
    random_shifts /= random_shifts.sum(axis=1)[:, np.newaxis]    # shift by numbers & find largest (by finding the smallest of the negative)
    shifted_probabilities = random_shifts - replicated_probabilities
    return np.argpartition(shifted_probabilities, sample_size, axis=1)[:, :sample_size]

def confidence_interval(count_a, count_b, count_ab, num_sequences, number_of_samples, sample_size=None, confidence=0.95,
                        draw_result=False):
    if not sample_size:  # It is common to resample using the same size as the original sample
        sample_size = float(num_sequences)  # ..but that won't always work well for rare mutations
    freqa = count_a / float(num_sequences)   # convert those counts into frequencies
    freqb = count_b / float(num_sequences)

    o_freq_ab = count_ab / num_sequences  # This is "O" -- the actual observed frequency of the mutations co-occuring

    e_freq_ab = freqa * freqb  # calculate the expected frequency of the mutations occuring together by chance

    p_a_alone = (count_a - count_ab) / num_sequences  # the A and B counts above include sequences where both mutations
    p_b_alone = (count_b - count_ab) / num_sequences  # occur so we calculate the frequency that mutation occurs just by
    p_ab_together = o_freq_ab                         # itself.
    p_none = 1 - (p_a_alone + p_b_alone + p_ab_together)

    print('Observed AB frequency', o_freq_ab)
    print('Expected AB frequency', e_freq_ab)
    print('Calculated AB O over E is:', o_freq_ab / e_freq_ab)
    print('Calculated log2(AB O over E) is:', log2(o_freq_ab / e_freq_ab))

    #  print('Sum of all probabilities (should be one): ', p_a_alone, '+', p_b_alone, '+', p_ab_together, '+', p_none,
    #      '=', p_a_alone + p_b_alone + p_ab_together + p_none)

    sampled_log_odds = []  # we will store the values for each sampled O / E value

    for i in range(number_of_samples):
        if not i % 10:
            print('Working on sample', i)
            

        elements = ['A', 'B', 'X', '_']  # X symbol here represents both mutations together, _ represents no mutations.
        probabilities = [p_a_alone, p_b_alone, p_ab_together, p_none]
        # np.random.choice is rate-limiting step, replacing w/ multidimension shifting/other method
        startTime = time.perf_counter()
        
        #Attempting alternative random choice with replacement
        #sample = multidimensional_shifting(num_samples, sample_size, elements, probabilities)
        
        
        sample = np.random.choice(elements, sample_size, p=probabilities)

        #sample = np.random.Generator.choice(elements, sample_size, replace=True, p=probabilities, axis=0, shuffle=False)

        
        startTime = time.perf_counter()
        sample_count_ab = (sample == 'X').sum()
        sample_count_a = (sample == 'A').sum() + sample_count_ab
        sample_count_b = (sample == 'B').sum() + sample_count_ab

        sample_freq_a = sample_count_a / sample_size
        sample_freq_b = sample_count_b / sample_size
        sample_freq_ab = sample_count_ab / sample_size

        
        try: 
            startTime = time.perf_counter()
            sampled_log_odds.append(log2(sample_freq_ab / (sample_freq_a * sample_freq_b)))


        except ValueError:
            print(sample_freq_ab, sample_freq_a, sample_freq_b)
            continue
            #quit('Numerical error, try increasing samples size') #Remove comment once done testing

    # print(sorted(sampled_log_odds))  # comment this out once you are happy with what is happening

    mu = np.mean(sampled_log_odds)
    sigma = np.std(sampled_log_odds)
    partial = (1 - confidence) / 2
    lower, upper = np.quantile(sampled_log_odds, [partial, confidence + partial])
    median = np.median(sampled_log_odds)

    if draw_result:

        # the histogram of the data
        num_bins = 25
        fig, axs = plt.subplots(1, 3) 
        n, bins, patches = axs[0].hist(sampled_log_odds, num_bins, density=True)

        # add a 'best fit' line

        y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
             np.exp(-0.5 * (1 / sigma * (bins - mu))**2))

        axs[0].plot(bins, y, '--')
        axs[0].set_xlabel('Log2 Odds')
        # axs[0].set_ylabel('(not actually) Probability Density')
        axs[0].set_title(r'Sampled Odds')
        axs[0].set_yticks([])
        axs[0].set_ylabel('Counts')

        # Draw a matplotlib violin plot
        axs[1].set_title('Violin plot')
        axs[1].violinplot(sampled_log_odds, positions=[1], showmeans=True, showmedians=False)
        axs[1].set_xticks([1])
        axs[1].set_xticklabels(['Pair 1'])
        axs[1].set_xlabel('Mutant Pairs')


        # Matplotlib error bar ... note, you can have LOTS of error bars in one plot if you like
        axs[2].set_title('Error bar plot')
        axs[2].set_xticks([1, 2])
        axs[2].set_xticklabels(['Pair 1', 'Pair 2'])
        axs[2].set_xlabel('Mutant Pairs')
        axs[2].errorbar(1, median, yerr=(upper-lower)/2, capsize=10, fmt='o')
        axs[2].errorbar(2, median, yerr=(upper-lower)/2, capsize=10, fmt='o')
        fig.tight_layout()
        plt.show()

    return lower, median, upper

# Example of using this function
# Note that results return are by default Log2 transformed
"""
sequences = 100000      # These numbers will all come from the COVIDCG database in practice
A_count = 1500           # count of sequences with mutation A  -- these given here are for example only
B_count = 2400           # count of sequences with mutation B
AB_count = 3            # count of sequences with both mutations

samp_size = sequences * 10  # The size of each sample, often this will be the same size as the original dataset
# However, because some mutants may be so rare that some samples will cuase a numerical
# error when calculating the odds, this might need to be larger.

samples = 100  # The total number of samples that will be taken. Ideally this would be 1000 or greater

low, med, high = confidence_interval(A_count, B_count, AB_count, sequences, samples, samp_size,
                                     confidence=0.95, draw_result=True)

print('Log2 lower:', low, 'log2 median:', med, 'Log2 upper:', high)
print('Real lower:', 2 ** low, 'log2 median:', 2 ** med, 'Log2 upper:', 2 ** high)
"""