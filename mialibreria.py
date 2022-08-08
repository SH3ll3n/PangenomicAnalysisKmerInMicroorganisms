import math
import matplotlib.pyplot as plt
import os
from pprint import pprint
import random

########################################
def estrai_kmeri(s,k):
    """
    Data una stringa s estrae le sottostringhe di lunghezza k (k-meri).
    ---------
    Parametri:
        s (str) : la stringa in input
        k (int) : la lunghezza delle sottostringhe
    ---------
    Ritorna:
        (set) l'insieme dei k-meri di s 
    """
    D = set()
    
    for i in range(0, len(s)-k+1):
        D.add( s[i:i+k] )
            
    return D
    
#################################


def leggi_sequenze_fasta(miofile):
    """
    Legge una o più sequenze fasta da un file e le ritorna come lista di stringhe.
    ----------
    Parametri:
        miofile (str) : percorso assoluto o relativo al file da leggere
    ----------
    Ritorna:
        list( (str) ) : stringhe lette nel file nello stesso ordine in cui sono scritte
    """
    ifile =  open(miofile,'r')

    stringhe = list()

    s = ''
    for l in ifile:
        #print(len(l))
        l = l.strip()
        if (len(l)>0):
            if l[0] == '>':
                if len(s)>0:
                    stringhe.append( s )
                    s = ''
            else:
                s += l
    if len(s)>0:
        stringhe.append( s )
        s = ''

    ifile.close()
    #print('@', len(stringhe))
    return stringhe



#################################

#stringhe = leggi_sequenze_fasta('Babesia')

#genoma = stringhe[0]
#if 'N' in genoma:
 #   print('ATTENZIONE carattere N in genoma')
    
#genoma[0:100]

#######################
random.seed(0)

def genera_random(l):
    """
    Genera una stringa random di lunghezza l sull'alfabeto A,C,G,T
    """
    random.seed(0)
    s = ''

    for i in range(l):
        x = random.randint(0,3)
        if x == 0:
            s += 'A'
        elif x == 1:
            s += 'C'
        elif x == 2:
            s += 'G'
        else:
            s += 'T'
    return s



########################################
def genera_random_noseed(l):
    """
    Genera una stringa random di lunghezza l sull'alfabeto A,C,G,T
    """
    random.seed(None)
    s = ''

    for i in range(l):
        x = random.randint(0,3)
        if x == 0:
            s += 'A'
        elif x == 1:
            s += 'C'
        elif x == 2:
            s += 'G'
        else:
            s += 'T'
    return s


########################################
def genera_random_seed(l, seed):
    """
    Genera una stringa random di lunghezza l sull'alfabeto A,C,G,T
    """
    random.seed(seed)
    s = ''

    for i in range(l):
        x = random.randint(0,3)
        if x == 0:
            s += 'A'
        elif x == 1:
            s += 'C'
        elif x == 2:
            s += 'G'
        else:
            s += 'T'
    return s



########################################
def mfl(s):
    """
    Calcola il valore di mfl (minimal forbidden length) della stringa s
    , ovvero il minimo valore di k per cui esiste un k-mero teorico 
    che non è presente nella stringa
    """
    k = 1
    while k < len(s):
        kmeri = estrai_kmeri(s,k)
        #print(k, len(kmeri), 4**k)
        if len(kmeri) < 4**k:
            return k
        k += 1
    return -1


########################################
def mcl(s):
    """
    Calcola il maximum complete lenght della stringa s.
    """
    return mfl(s)-1






########################################
def molteplicita(s,w):
    """
    Calcola la molteplicita' di w in s.
    ----------
    Parametri:
        s (str)
        w (str)
    ----------
    Ritorna:
        (int) la molteplictia' di w in s
    """
    molt = 0
    for i in range(len(s) - len(w) + 1):
        if s[i:i+len(w)] == w:
            molt += 1
    return molt


########################################
def get_mults(s,k):
    """
    Ritorna la lista dei k-meri della stringa S arrichita delle loro molteplicità.
    ----------
    Parametri:
        s (str) : la stringa data
        k (int) : la lunghezza dei k-meri
    ----------
    Ritorna:
        (dict) : un dizionario in cui le chiavi sono i k-mer presenti nella stringa e i valori sono le loro molteplicità
    """
    mults = dict()
    # chiavi  = k-meri
    # valori = molteplicità
    for i in range(len(s)-k+1):
        kmer = s[i:i+k]
        mults[ kmer ] = mults.get(kmer, 0) + 1
    return mults



########################################
def mhl(s):
    """
    Calcola il minimal hapax length della stringa s, 
    ovvero la lunghezza dell'hapax più corto.
    """
    k = 1
    while k < len(s):
        molts = get_mults(s,k)
        for kmero in molts:
            if molts[kmero] == 1:
                return k
        k += 1
    return -1


########################################
def mhl_info(s):
    """
    Calcola il minimal hapax length della stringa s, 
    ovvero la lunghezza dell'hapax più corto.
    ----------
    Parametri:
        s (str) : la stringa su cui calcolare mhl
    ----------
    Ritorna:
        (int) : mhl(s)
        (str) : il primo k-mero hapax
    """
    k = 1
    while k < len(s):
        molts = get_mults(s,k)
        for kmero in molts:
            if molts[kmero] == 1:
                return k,kmero
        k += 1
    return -1,''


########################################
def mrl(s):
    """
    Calcola la lunghezza del repeat più lungo in s
    ----------
    Parametri:
        s (str) : la stringa su cui calcolare mrl
    ----------
    Ritorna:
        (int) : mrl(s)
    """
    k = 1
    trovato_repeat = False
    while k < len(s):
        trovato_repeat = False
        molts = get_mults(s,k)
        
        for kmero in molts:
            if molts[kmero] > 1:
                trovato_repeat = True
                break
                
        if trovato_repeat == False:
            return k-1
        k += 1
    return -1



########################################
def entropia(s,k):
    """
    Calcola l'neoptria empitica k-esima dei k-meri in s.
    """
    molts = get_mults(s,k)
    #tot = 0
    #for kmer,m in molts.items():
    #    tot += m
    tot = sum(molts.values())
    en = 0.0
    for kmero,m in molts.items():
        en += (m/tot) * math.log( (m/tot), 2 )
    en *= -1
    return en


########################################
def leggi_sequenze_fasta(miofile):
    """
    Legge una o più sequenze fasta da un file e le ritorna come lista di stringhe.
    ----------
    Parametri:
        miofile (str) : percorso assoluto o relativo al file da leggere
    ----------
    Ritorna:
        list( (str) ) : stringhe lette nel file nello stesso ordine in cui sono scritte
    """
    ifile =  open(miofile,'r')

    stringhe = list()

    s = ''
    for l in ifile:
        #print(len(l))
        l = l.strip()
        if (len(l)>0):
            if l[0] == '>':
                if len(s)>0:
                    stringhe.append( s )
                    s = ''
            else:
                s += l
    if len(s)>0:
        stringhe.append( s )
        s = ''

    ifile.close()
    #print('@', len(stringhe))
    return stringhe


########################################
def visualizza_distr(x, y, xlabel=None, ylabel=None, xticks=None, title=None, specie=None):
    plt.bar(x , y)
    if xticks:
        plt.xticks(x, xticks)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)
        plt.savefig(f"./immagini/{specie}/{title}.png")
    plt.show()
