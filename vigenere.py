import sys
import itertools
import argparse
import re
import numpy as np
from nltk import ngrams
import operator
from collections import Counter
import matplotlib.pyplot as plt

fn_factors = "vigenere_factors.pdf"
fn_letters = "Letter_distribution.pdf"
fn_rotation = "Rotation_detection.pdf"

alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Monograms taken from
# http://practicalcryptography.com/media/cryptanalysis/files/english_monograms.txt

monograms_off = b"""E 529117365
T 390965105
A 374061888
O 326627740
I 320410057
N 313720540
S 294300210
R 277000841
H 216768975
L 183996130
D 169330528
C 138416451
U 117295780
M 110504544
F 95422055
G 91258980
P 90376747
W 79843664
Y 75294515
B 70195826
V 46337161
K 35373464
J 9613410
X 8369915
Z 4975847
Q 4550166"""

def eta_inv(number):
    '''
    Takes int, returns ASCII
    '''
    number = number % 26
    return chr(int(number)+65)

def eta(letter):
    '''
    Takes ASCII, returns int
    '''
    return ord(letter)-65

def rotate_alphabet(letter,length=0):
    '''
    Takes a letter and shifts it by length along the alphabet
    '''
    number = ord(letter)-65
    number = (number + length) % 26
    return chr(int(number)+65)


def factors(n,limit=False,skipone=True):
    '''
    Generate factors of a given integer n up to limit
    E.g.: factors(20,8) == [2,4,5]
    '''
    if not limit: limit = n

    if skipone == True: start = 2
    else: start = 1

    result = []

    for i in range(start, limit + 1):
        if n % i == 0:
            result.append(i)

    return result

parser = argparse.ArgumentParser(description='Vigener En-/Decryption')
group = parser.add_mutually_exclusive_group()
group.add_argument('-d', '--decrypt', help="Decrypt, Reading from stdin", default=False, action="store_true")
group.add_argument('-e', '--encrypt', help="Encrypt, Reading from stdin", default=False, action="store_true")
parser.add_argument('-k', '--keylength', help="Determine keylength. If 0, keylength will be guessed.; default:0", default=0, type=int)
parser.add_argument('-K', '--key', help="Define de-/encryption key; default:0", default="", type=str)
parser.add_argument('-l', '--keylengthlimit', help="Maximum keylength taken into consideration; default:20", default=20, type=int)
parser.add_argument('-p', '--plot', help="Generate plots", default=False, action="store_true")
parser.add_argument('--threshholdkeylength', help="Tolerance/threshhold for determing min/max; default:0.15", default=0.15, type=float)
parser.add_argument('--threshholdrotation', help="Tolerance/threshhold for determing min/max; default:0.15", default=0.15, type=float)
parser.add_argument('-v', '--verbose', help="More talkative", default=False, action="store_true")
args = parser.parse_args()

if args.decrypt:

    # Reading encrypted material from stdin
    crypt = sys.stdin
    crypt = str(*crypt).strip()

    verify = [ord(l) for l in crypt if ord(l) <= 90 and ord(l) >= 65]
    if len(verify) == 0:
        raise RuntimeError('No input given!')
    elif not len(verify) == len(crypt):
        raise RuntimeError('Only acception capital letters. Non-capital letters found! Use e.g. tr [a-z] [A-Z].')

    if args.key == "":

        # Read monograms_off into dict, sort it by letter and slice it into
        # two lists with keys, values
        monograms_off = np.genfromtxt(monograms_off.splitlines(),dtype=None)
        monograms_off = dict(monograms_off)
        monograms_off = sorted(monograms_off.items(), key=operator.itemgetter(0), reverse=False)
        monox,monoy = zip(*monograms_off)
        temp_sum = np.sum(monoy)
        # Relative frequencies
        monoy = [m/temp_sum for m in monoy]

        if not args.keylength:
            # Generating list of tri- and quadgrams from crypt
            trigrams = [''.join(list(g)) for g in ngrams(crypt, 3)]
            quadgrams = [''.join(list(g)) for g in ngrams(crypt, 4)]
            
            # Count occurences of ngrams
            #monograms = Counter(crypt)
            Ntrigrams = len(trigrams)
            trigrams = Counter(trigrams)
            Nquadgrams = len(quadgrams)
            quadgrams = Counter(quadgrams)
            
            # Sort monograms by alphabet
            #monograms = sorted(monograms.items(), key=operator.itemgetter(0))
            #letters, frequencies = zip(*monograms)
            #temp_sum = np.sum(frequencies)
            #frequencies = [m/temp_sum for m in frequencies]
            #monograms = {"keys": letters, "values": frequencies}

            spectrum = {}
            
            tmp = sorted(trigrams.items(), key=operator.itemgetter(1), reverse=True)
            letters, frequencies = zip(*tmp)
            frequencies = [m/Ntrigrams for m in frequencies]
            spectrum["trigrams"] = {"keys": letters, "values": frequencies}

            tmp = sorted(quadgrams.items(), key=operator.itemgetter(1), reverse=True)
            letters, frequencies = zip(*tmp)
            frequencies = [m/Nquadgrams for m in frequencies]
            spectrum["quadgrams"] = {"keys": letters, "values": frequencies}

            keylength = {}

            # For every set of ngrams (tri-, quad-, ...)
            for key in spectrum:
                # For every ngram ...
                for ngram in spectrum[key]["keys"]:
                    # ... find all occurences in crypt
                    occ = [m.start() for m in re.finditer(ngram, crypt)]
                    for i in reversed(range(len(occ)-1)):
                        # Take the pairwise distance in letters and generate
                        # factors until --keylenghtlimit
                        for factor in factors(occ[i+1]-occ[i],args.keylengthlimit):
                            # Sum the occurences of the factors up in dict keylength
                            if not factor in keylength: keylength[factor] = 1
                            else: keylength[factor] += 1
            

            if len(keylength.items()) == 0:
                raise RuntimeError("Too few input to determine key!")

            # Sort dict keylength by keys
            keylength = sorted(keylength.items(), key=operator.itemgetter(0), reverse=False)
            plotx,ploty = zip(*keylength)
            #temp_sum = np.sum(ploty)
            #ploty = [m/temp_sum for m in ploty]

            # Take the topmost keylength acc to --threshhold
            factormax = np.max(ploty)
            factor_threshh = (1-args.threshholdkeylength) * factormax 
            factor_candidates = [plotx[ploty.index(f)] for f in ploty if f > factor_threshh]

            if args.plot:
                # Plot keylength spectrum
                fig = plt.figure()
                ax = fig.add_subplot(111)
                pos = range(len(plotx))
                ax.set_xticks(pos)
                ax.set_xticklabels(plotx)
                ax.set_title('Vigenere keylength analysis')

                # Actual plotting
                ax.bar(pos, ploty, color='r')
                label, = ax.plot([pos[0], pos[-1]], [factor_threshh, factor_threshh], "k--",label='Threshhold at {} %'.format(int(args.threshholdkeylength*100)))

                ax.legend(handles=[label],loc='best')

                fig.savefig(fn_factors)
                if args.verbose:print("Vigenere keylength analysis saved to {}".format(fn_factors))

            if len(factor_candidates) == 1:
                print("Keylength probably {}".format(factor_candidates[0]))
                keylength = factor_candidates[0]
            else:
                print("Keylength probably one of {}".format(factor_candidates))
                keylength = np.max(factor_candidates)
                print("Taking {}. Otherwise specify via --keylength".format(keylength))

        else:keylength = args.keylength

        if args.verbose:print("Factor analysis saved in {}".format(fn_factors))

        # positions: Dict with every keylength letter from crypt
        # E.g.: keylength = 7: crypt[0], crypt[7], ...
        positions = {}

        for j in range(keylength):
            positions[j] = [crypt[i] for i in np.arange(j, len(crypt), keylength)]

        monograms = {}
        # Generate letter-spectrum
        for j in range(keylength):
            monograms[j] = Counter(positions[j])

        # Check if complete alphabet present
        for letter in alphabet:
            for j in range(keylength):
                # If letter missing, set occ to zero
                if letter not in monograms[j]:monograms[j][letter]=0


        passphrase = {}

        for letter in range(keylength):
            # Sort spectrum by letters
            monograms[letter] = sorted(monograms[letter].items(), key=operator.itemgetter(0))
            letters, frequencies = zip(*monograms[letter])
            temp_sum = np.sum(frequencies)
            frequencies = [m/temp_sum for m in frequencies]
            monograms[letter] = {"keys": letters, "values": frequencies}

            if args.plot:
                # Generate figure for plotting letter distribution
                fig = plt.figure()
                ax = fig.add_subplot(111)
                pos = range(len(letters))
                ax.set_xticks(pos)
                ax.set_xticklabels(letters)
                ax.set_title('Letter distribution of letter {}'.format(letter))

                ax.bar(pos, frequencies, color='r')
                fig.savefig("".join(["letter{}".format(letter),fn_letters]))
                if args.verbose:print("Letter distribution for letter {} saved to {}".format(letter,"".join(["letter {}".format(letter),fn_letters])))

            # Rotation detection: Which Caesar rotation has the minimum
            # deviation from monograms_off/the english language
            rot_detect = {}
            # Go though all 26 shifts
            for rot in range(26):
                iii=0
                for value in monograms[letter]["values"]:
                    # Take sum the absolut values of the difference of the
                    # current distribution and the distr. of the english lang.
                    if not rot in rot_detect:
                        rot_detect[rot] = np.abs(value-monoy[ (iii+rot) % len(monoy) ])
                    else:
                        rot_detect[rot] += np.abs(value-monoy[ (iii+rot) % len(monoy) ])
                    iii += 1

            length, values = zip(*rot_detect.items())

            # Take the lowest deviation acc to --threshhold
            deviation_min = np.min(values)
            deviation_threshh = (1+args.threshholdrotation) * deviation_min
            deviation_candidates = [length[values.index(f)] for f in values if f < deviation_threshh]

            if args.plot:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                pos = range(len(length))
                ax.set_xticks(pos)
                ax.set_xticklabels(length)
                ax.set_title('Rotation detection for letter {}:{}'.format(letter,[rotate_alphabet('A',-l) for l in deviation_candidates]))

                ax.bar(pos, values, color='r')
                label, = ax.plot([pos[0], pos[-1]], [deviation_threshh, deviation_threshh], "k--",label='Threshhold at {} %'.format(int(args.threshholdrotation*100)))
                ax.legend(handles=[label],loc='best')
                fig.savefig("".join(["letter{}".format(letter),fn_rotation]))
                if args.verbose:print("Rotation detection for letter {} saved to {}".format(letter,"".join(["letter {}".format(letter),fn_letters])))

            passphrase[letter] = deviation_candidates

        passwords = list(itertools.product(*passphrase.values()))

        print("Possible passwords:")
        for pw in passwords:
            print("".join([rotate_alphabet('A',-l) for l in pw]))

    else:

        if args.key == "":
            raise ValueError("No encryption key given!")
        
        PW = args.key
        
        PW=(int(len(crypt)/len(PW))+1)*PW
        PW=PW[:len(crypt)]

        decrypt = []
        for i,j in zip(crypt,PW):
            decrypt.append(eta_inv(eta(i)-eta(j)))

        decrypt = "".join(decrypt)
        print(decrypt)
        
elif args.encrypt:

    if args.key == "":
        raise ValueError("No encryption key given!")

    # Reading material from stdin
    MS = sys.stdin
    MS = str(*MS).strip()

    verify = [ord(l) for l in MS if ord(l) <= 90 and ord(l) >= 65]
    if len(verify) == 0:
        raise RuntimeError('No input given!')
    elif not len(verify) == len(MS):
        raise RuntimeError('Only acception capital letters. Non-capital letters found! Use e.g. tr [a-z] [A-Z].')


    PW = args.key
    
    PW=(int(len(MS)/len(PW))+1)*PW
    PW=PW[:len(MS)]

    crypt = []
    for i,j in zip(PW,MS):
        crypt.append(eta_inv(eta(i)+eta(j)))

    crypt = "".join(crypt)

    print(crypt)

else:
    sys.exit(1)
