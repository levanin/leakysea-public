from collections import Counter
from numpy import mean
import csv

# Attacking Seasign modifications in https://ieeexplore.ieee.org/document/10416853/
# Which involve presampling ephemeral exponents in a biased manner instead of rejection sampling
# Author: Shai Levin

# Parameters
#================================================================================================

# Bound B which defines the interval [-B,B] for sampling secret exponents
B = 5
# Dimension n for CSIDH lattice
n = 74

# Parameters:
# delta: Bound delta which defines the interval [-(delta+1)B, (delta+1)B] for sampling ephemeral exponents
# t: Number of repetitions t in 1 execution

# Parameter set I from https://ieeexplore.ieee.org/document/10416853/
delta = 9472 # parameter set I
t = 128

# Parameter set II from https://ieeexplore.ieee.org/document/10416853/
# delta = 114
# t = 337

# Number of known signatures where you'd expect to see a large enough sample size from the distribution
optimal_sigs=floor(4*delta*B/t)

# Fixed key for testing
FIXED_KEY = vector(ZZ, [5, 1, -3, 0, 5, -4, 5, -1, -4, 3, 1, 3, 2, 2, -3, -4, 5, 0, 3, -5, -1, -4, 0, 2, -2, -5, 1, -1, -1, 3, -2, 2, 2, 2, 4, 3, -5, -1, -3, -4, -5, -3, -5, -1, 4, 2, 2, 5, -5, 4, -4, 1, 3, 3, 4, -3, 4, 0, -3, 4, -4, 5, -5, -2, -4, -4, 2, 3, 4, 2, 0, 0, 5, -1])
#================================================================================================


# Sampling function - might want to make this cryptographically secure vs just using random.randint
def sample_exponent_vector(bound):
    return vector(ZZ,[randint(-bound, bound) for _ in range(n)])

# Generate secret exponent in [-B,B]^n
def gen_secret():
    return sample_exponent_vector(B)

# Generate ephemeral exponent in [-(delta+1)B, (delta+1)B]^n
def gen_ephemeral():
    return sample_exponent_vector((delta+1)*B)

# Generate biased ephemeral exponent (which is conditioned on secret) in [-(delta+1)B, (delta+1)B]^n
def gen_biased_ephemeral(secret):
    while True:
        ephemeral = gen_ephemeral()
        test = ephemeral - secret
        if min(test) >= -delta*B and max(test) <= delta*B:
            return ephemeral
        
def attack(biased_ephemerals):
    # Given biased ephemeral exponents, attempt to guess the secret exponent
    guess = []
    for i in range(n):
        # Determines the i-th entry of the secret exponent vector

        # Collect the i-th entries of all the vectors
        i_entries = [v[i] for v in biased_ephemerals]
        # Determine the minimum and maximum values of the i-th entries
        minval, maxval = min(i_entries), max(i_entries)
        # Shift by the bound to get the relative distance from the bound.
        b,a=maxval - delta*B, minval + delta*B
        if a == b:
        #We are certain of this guess!
            guess.append(a)
        else:
        #The guess lies somewhere between b and a. Note that b <= a since the positive values are on the left of the relative bound and the negative values are on the right of the relative bound. See plots for intuition.
            negbound = max(-B,b)
            posbound = min(B,a)
            if negbound == posbound:
                guess.append(negbound)
            else:
                guess.append(list(range(negbound,posbound+1)))
    return guess

def compare_guess(guess,secret):
    # Compare the guessed secret with the actual secret
    # If the guess is correct, return True, the number of determined entries and the keyspace
    known_exponents = 0
    keyspace = 1
    for i in range(n):
        if type(guess[i]) == type(secret[i]):
            if guess[i] == secret[i]:
                known_exponents += 1
            else:
                return False,-1,-1
        else:
            if secret[i] in guess[i]:
                keyspace *= len(guess[i])
            else:
                return False,-1,-1
    return True,known_exponents,keyspace


def experiment(views=[],fixed_key=False,known_sigs = optimal_sigs):

    # Generate secret exponent
    print("Secret exponent is:")
    if fixed_key: # For testing purposes
        secret = FIXED_KEY
    else:
        secret = gen_secret()
    print(secret)
    print("====================================")
    # Generate biased ephemeral exponents. This assumes that given t repetitions, 1/2 of them are 0 challenges (which are the biased cases.)
    biased_ephemerals = [gen_biased_ephemeral(secret) for _ in range(floor(known_sigs*t/2))]
    # If the function is called with views, we plot the distribution of good samples for the given vector entries
    if len(views) != 0:
        good_samples = graph_bias(views, biased_ephemerals)
        print("Plotting distribution of good samples from {} signatures in vector position {}, note that e_{} = {}. The positive elements should be on the left of the relative bound, and negative on the right".format(known_sigs,views,views,[secret[view] for view in views]))
        for i,view in enumerate(views):
            minusplot = list_plot([(res,freq) for (res,freq) in good_samples[i][0].items()],title="View of good distribution of v_{} = {}".format(view,secret[view]),legend_label="negative values")
            plusplot = list_plot([(res,freq) for (res,freq) in good_samples[i][1].items()],color='red',title="View of good distribution of v_{} = {}".format(view,secret[view]),legend_label="positive values")
            pole = line([(secret[view],0),(secret[view],known_sigs/200)],color='green',legend_label="relative bound")
            show(minusplot+plusplot+pole)
    guessed_secret = attack(biased_ephemerals)
    print("Guess for secret exponent is:")
    print(guessed_secret)
    print("====================================")
    is_contained,known_exponents,keyspace = compare_guess(guessed_secret,secret)
    if is_contained:
        print("Given {} signatures, the attack has recovered {} known exponents, overall keyspace reduced to size 2^{:.2f}".format(int(known_sigs), known_exponents,float(log(keyspace,2))))
        print("====================================")
    else:
        print("Possible key space does not contain the secret key!")
    return float(log(keyspace,2))



def graph_bias(views, biased_ephemerals):
    #Graph the distribution of 'good' samples
    good_samples = []
    for i in views:
        #We are interested in samples which are 'good' in the sense that they leak information about the secret key. This is when the samples lie in the bounds below
        i_plus_entries = []
        i_minus_entries = []
        for v in biased_ephemerals:
            if v[i] >= (delta-1)*B:
                #Normalise the values so they can be graphed nicely
                i_plus_entries.append(v[i]-delta*B)
            elif v[i] <= -(delta-1)*B:
                i_minus_entries.append(v[i]+delta*B)
        #Cumulative sums for the histogram
        good_samples.append([Counter(i_minus_entries),Counter(i_plus_entries)])
    return good_samples

def testing_harness(num_iters, bound=4,step=2):
    #Repeats the experiment for different values of known signatures and computes keyspaces.
    data = []
    with open('data.csv','w') as csvfile:
        csvwriter = csv.writer(csvfile)
        for num_sigs in range(1,bound*optimal_sigs,step):
            i_vals = []
            for _ in range(num_iters):
                i_vals.append(experiment(known_sigs=num_sigs))
            row = (num_sigs,numpy.mean(i_vals))
            data.append(row)
            csvwriter.writerow(row)
    return data

def plot_keyspace(filename):
    with open(filename, 'r') as file:
        data = list(csv.reader(file))
    line = list_plot(data,color='red',linestyle="--",plotjoined=True)
    boxes = list_plot(data,plotjoined=False,marker='s',color='red')
    (line+boxes).show(dpi=1000)

#Views show the distribution of good samples for the given vector entries
experiment()