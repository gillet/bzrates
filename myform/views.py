#-*-coding:utf8-*-
"""
==================================================================================
bz-rates: a web-tool to accurately estimate mutation rates from fluctuation assays
==================================================================================
Copyright (C) 2015 Alexandre Gillet-Markowska

This file is part of bz-rates.

bz-rates is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

bz-rates is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with bz-rates.  If not, see <http://www.gnu.org/licenses/>.
==================================================================================
Contact: Alexandre Gillet-Markowska - alexandre(dot)gillet(at)yahoo(dot)fr
==================================================================================
"""
#-------------------------- Dependencies ------------------------
from django.forms import ModelForm
from django.shortcuts import render_to_response, HttpResponseRedirect, HttpResponse
from django import forms
from django.template import RequestContext
import re
import numpy as np
from scipy import integrate
from scipy.stats import norm


class NMutField(forms.CharField):
    def to_python(self, value):
        """Normalize data to a list of strings."""
        Split = [filter(None, re.split(r" |\t", elt)) for elt in re.split("\r\n?|\n", value.strip())]
        Flat = [item for sublist in Split for item in sublist]

        if not set([len(elt) for elt in Split]) == set([2]):
            raise forms.ValidationError("2-columns only: Nmutant | Nplated")

        if not all(map(is_integer, Flat)):
            raise forms.ValidationError("Only integers are allowed")

        if not all([int(float(x[0]))<=int(float(x[1])) for x in Split]):
            raise forms.ValidationError("The number of mutants 'Nmutant' cannot be larger than the number of plated cells 'Nplated'")

        return value



class MutForm(forms.Form):
    N0 = forms.IntegerField(help_text="Initial number of cells in the culture", initial=1, label="N0", widget=forms.NumberInput(attrs={'class': 'narrow-select', 'required': "", 'min': 1}))
    b = forms.DecimalField(help_text="Mutant fitness relatively to wild type", label="b", widget=forms.NumberInput(attrs={'class': 'narrow-select', 'min': 10e-100}), required=False)
    z = forms.DecimalField(help_text="Plating efficiency", initial=1, label="z", widget=forms.NumberInput(attrs={'class': 'narrow-select', 'required': "", 'min': 10e-100, 'max': 1}))    
    #fluctuation = NMutField(help_text="Nmutant | Nplated cells", widget=forms.Textarea(attrs={'class': 'dataFluc', 'required': ""}), label="")
    fluctuation = forms.CharField(help_text="Nmutant | Nplated cells", widget=forms.Textarea(attrs={'class': 'dataFluc', 'required': ""}), label="")



#______________________________________________________________________________
def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
#______________________________________________________________________________
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#______________________________________________________________________________
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
#______________________________________________________________________________
def is_integer(s):
    """Looks for numbers with no decimal part"""
    try:
        c = float(s)
        return int(c) == c
    except ValueError:
        return False

#______________________________________________________________________________
# Estimateur de Lea et Coulson (Méthode de la médiane)
def LeaCoulson(Nmut , prec=6, a = 0., b = 30.):
    """ Calcul de m selon la méthode de la médiane (Lea et Coulson)
    Input : Nmut = mutants par culture
           prec = nombre de décimales pour tester l'égalité
           a = borne inférieure de l'intervalle à tester
           b = borne supérieure de l'intervalle à tester
           """
    r0 = np.median(Nmut)
    try:
        assert r0 != 0
    except AssertionError:
        #return None
        return 0
    #print 'La médiane r0 vaut 0'
    else:
        m = float(a + b)/2
        expr = round(r0/m - np.log(m) - 1.24, prec) #(fonction décroissante de m)
        while expr != 0:
            if expr > 0:
                a = m
            elif expr < 0:
                b = m
            m = (a + b)/2
            expr = round(r0/m - np.log(m) - 1.24, prec)
        return round(m , 5)

#______________________________________________________________________________
# Fonction de vraisemblance (son max est atteint pour le lot de paramètres (ici : m)
# qui explique le mieux la distribution).
def likelihood(p,cdistrib,robs):
    """ Donne la fonction de vraisemblance associée à la distribution des pr et des cr"""
    L = 1
    for r in robs:
        L *= p[r]**cdistrib[r]
    return L
#______________________________________________________________________________
# Log de la fonction de vraisemblance. Équivalent à l'autre (sauf vitesse de calcul ??)
def loglikelihood(p,cdistrib,robs):
    """Fonction de vraisemblance passée au logarithme"""
    L = 0
    for r in robs:
        L += np.log(p[r])*cdistrib[r]
    return L
#______________________________________________________________________________
# Calcul de la liste des probabilités selon la méthode de Ma-Sandri-Sarkar
def probaLCFDiff(m,rmax, b, N0, Nt):
    """Équation (4) de Jaeger & Sarkar 1995 Genetica : formule de récurrence de pr
    avec diff growth des mutants (b = fitness des mutants 0 < b < Inf).
    Renvoit une liste contenant les pr
    pr = probabilité qu'une culture ait r mutants"""
    a = (N0/Nt)**b
    p = [np.exp((m/b)*(a-1))] # initialisation
    for r in range(1,rmax+1):
        termesom = [(m/(r*b))*p[i]*(r-i)* \
                    ((((1-a)**(r-i))/(r-i)) - (((1-a)**(r-i+1))/(r-i+1))) \
                    for i in range(r)] #fixed (5) from Jaeger & Sarkar 1995

        p += [sum(termesom)]
    return p


#______________________________________________________________________________
# Ma-Sandri-Sarkar Maximum Likelihood Estimator
def LCFDiff(N0, Nt, b, m0, Nmut, prec = 5):
    """ Calcule m : nombre de mutations attendues dans chaque culture,
    selon la méthode du maximum de vraisemblance de LCF with Finite Number of cell Divisions with
    multiple seed and differential fitness
    Input :	m0 = m préliminaire
            b = fitness des mutants, 0 < b < Inf
            N0 = seed (nombre initial de cellules dans la culture)
            Nt = nombre final de cellulles dans la culture etallee
            Nmut = distribution du nombre de mutants par culture.
            prec = nombre de chiffres significatifs dans le résultat.
    """
    deft = {}
    deft['rmax'] = 500
    rmax = min(deft['rmax'] , max(Nmut)) # r : nb de mutants par culture
    if rmax == 0:
        m = 0
    else:
        robs = list(set(Nmut)) # Uniquement les valeurs de r observées, dans l'ordre
        robs.sort()
        for i,r in enumerate(robs): # On remplace les r observés supérieurs à rmax par rmax
            if r > rmax:
                robs[i] = rmax
        cdistrib = [Nmut.count(r) for r in range(rmax+1)] # Nombre de cultures avec r mutants
        m = round(m0,5)
        prectest = -2 # précision du résultat (nb chiffres signif.) en cours de construction
        while prectest < prec:
            Lcourbe = [] # Contiendra les valeurs de la fonction de vraisemblance correspondant à chaque m testé.
            mtest = list(np.arange( max(0 , m - 10**(-prectest)) , m + 10**(-prectest) , 10**(-prectest-1))) #attention, invalide si l'ecart initial entre 'm' et m0 >10
            for m in mtest:
                #p = proba(m,rmax)
                p = probaLCFDiff(m, rmax, b, N0, Nt)
                Lcourbe.append(loglikelihood(p,cdistrib,robs))
            m = mtest[Lcourbe.index(max(Lcourbe))]
            prectest += 1

    #print "MSS-MLE :	m0 = %s	m = %s" % (round(m0 , 4) , m)
    return m


#______________________________________________________________________________
# fonction de correction de la fraction étalée et/ou de la plating efficiency
def corr(m0 , z):
    """Applique une correction de la fraction étalée sur m0 (estimateur p0)"""
    if z == 1:
        mcor = m0
    else:
        mcor = m0 * (z - 1) / (z * np.log(z))
    return mcor

#______________________________________________________________________________
# écart-type de ln(m)
def sd_ln(m , c):
    """Input : m , c (nb de cultures)
    output : écart-type de ln(m)
    Formule : Foster 2006"""
    return 1.225 * m**(-0.315) / np.sqrt(c)
#______________________________________________________________________________
# Intervalle de confiance. 	Voir Foster 2006 (équations 29 et 30)
#						et Stewart 1994 ( équation 4)
def CL(M , m , c ):
    """ Intervalle de confiance de M """
    if M == 0:
        #CL = None,None
        CL = [0,0]
    else:
        sd = sd_ln(m , c)
        ecart = [-1.96 * sd * np.exp(1.96*sd*0.315) , 1.96 * sd * np.exp(1.96*sd*0.315)] # Foster 2006 (29-30)
        CL = [ M * np.exp(x) for x in ecart] # Stewart 1994 (4)
        CL = [round(x , 13) for x in CL]
    return tuple(CL)

#______________________________________________________________________________
#
def CalculateAllLCFDiff(Nmut , N0, Np, b, z):
    """Calcule m et M selon les estimateurs p0, L&C, et MSS,
    avec les intervalles de confiance.
    Nt peut être une liste ou un nombre."""
    c = len(Nmut)
    p0 = float(Nmut.count(0)) / c #proportion de cultures sans mutants
    if p0 ==0: #m0 cannot be calculated if 0 cultures with no mutants
        m0, m0cor = 0, 0
    else:
        m0 = - np.log(p0) #estime m0 avec le terme zero de la loi de poisson
    mLC = LeaCoulson(Nmut) #estime m0 avec LC
    m = mLC or m0 # 1è estimation pour commencer le calcul

    Npsd = np.std(Np)
    Npmoy = np.mean(Np)
    Ntmoy = Npmoy/z

    m = LCFDiff(N0, Ntmoy, b, m, Nmut)
    M = m / Npmoy

    mcorr = corr(m , z)
    Mcorr = mcorr / Ntmoy

    CLinf , CLsup = CL(Mcorr , corr(m , z) , c)

    return ["Maximum Likelihood", \
            float('{0:.2e}'.format(m)), round(M,13), \
            float('{0:.2e}'.format(mcorr)), round(Mcorr,13), \
            '{0:.3e}'.format(round(CLinf,13)) , '{0:.3e}'.format(round(CLsup,13)), \
            '{0:.3e}'.format(int(round(Npmoy))) , '{0:.3e}'.format(int(Npsd))]

def GF(Nmut, N0, Np, z):
    """Compute mutation rates with GF"""

    Npsd = np.std(Np)
    Npmoy = np.mean(Np)
    Ntmoy = np.mean(Np)/z

    m, rho = GFEst(Nmut)
    M = m / Npmoy

    mcorr = corr(m, z)
    Mcorr = corr(m , z) / Ntmoy

    CL = GFConfint(Nmut)
    #CLinf, CLsup = min(CL['alpha'].values()), max(CL['alpha'].values())
    mCLinf, mCLsup = min(CL['alpha'].values()), max(CL['alpha'].values())
    CLinf, CLsup = (mcorr-mCLinf)/Ntmoy,(mcorr+mCLsup)/Ntmoy

    rhoinf, rhosup = min(CL['rho'].values()), max(CL['rho'].values())

    return ["Generating Function", \
            float('{0:.2e}'.format(m)), round(M,13), \
            float('{0:.2e}'.format(mcorr)), round(Mcorr,13), \
            '{0:.3e}'.format(round(CLinf,13)) , '{0:.3e}'.format(round(CLsup,13)), \
            round(rho,3), round(rhoinf,3), round(rhosup,3), \
            '{0:.3e}'.format(int(round(Npmoy))) , '{0:.3e}'.format(int(Npsd))]

#______________________________________________________________________________
# """-----------------------------------------------------------------
#
# Generating Function
#
# Original R version 1.0 03/15/2012 by
# A. Hamon and B. Ycart
# Universite de Grenoble and CNRS, France
# http://ljk.imag.fr/membres/Bernard.Ycart/LD/
#
# Translated to Python 2.7 04/01/2015 by
# A. Gillet-Markowska
# Universite Pierre et Marie-Curie
# www.lcqb.upmc.fr/masstor
#
# -----------------------------------------------------------------
# """

#-------------------------- Data sets ----------------------------

# From p.504 of D. E. Luria and M. Delbruck, Mutations  of bacteria from
# virus sensitivity to virus resistance, Genetics, 28, (1943), 491-511
LD43a = (10,18,125,10,14,27,3,17,17,29,41,17,20,31,30,7,17,30,10,40,45,\
          183,12,173,23,57,51,6,5,10,8,24,13,165,15,6,10,38,28,35,107,13)
LD43b = (1,0,3,0,0,5,0,5,0,6,107,0,0,0,1,0,0,64,0,35,\
          1,0,0,7,0,303,0,0,3,48,1,4)
LD43c = (0,0,0,0,8,1,0,1,0,15,0,0,19,0,0,17,11,0,0)
# From L. Boe, T. Tolker-Nielsen, K. M. Eegholm, H. Spliid and A. Vrang,
# Fluctuation analysis of mutations to nalidixic acid resistance in
# Escherichia Coli, J. Bacteriol., 176(10) (1994) 2781-2787
B94 = list(np.hstack([[0]*543,[1]*169,[2]*92,[3]*57,[4]*42,[5]*25,\
        [6]*23,[7]*13,[8]*11,[9]*5,[10]*8,[11]*8,\
        [12]*9,[13]*5,[14]*5,15,[16]*8,[17]*2,[18]*5,\
        [19]*6,[20]*2,[21]*5,23,[24]*2,25,[26]*5,[27]*3,\
        28,[29]*3,[30]*3,31,32,34,35,36,37,[39]*2,40,41,41,49,52,\
        57,59,66,68,69,[73]*2,74,105,107,116,132,137,140,146,151,152,\
        192,258,265,320,482,[513]*4]))
# From p.22 of W. A. Rosche and P. L. Foster, Determining mutation rates
# in bacterial populations, Methods, 20(1), (2000), 1-17
RF00 = list(np.hstack([[0]*11,[1]*17,[2]*12,[3]*3,[4]*4,5,6,[7]*2,9]))
# From p.239 of Q. Zheng, Statistical and algorithmic methods for
# fluctuation analysis with SALVADOR as an implementation,
# Mathematical Biosciences, 176, (2002), 237-252
Z02 = [33,18,839,47,13,126,48,80,9,71,196,66,28,17,27,37,\
        126,33,12,44,28,67,730,168,44,50,583,23,17,24]


#-------------------- Generating Functions -----------------------
def hr(rho, z):
    """returns the evaluation at z in [0,1], of the generating
    function of the Yule distribution with parameter rho
    """
    eps = 10**(-8)              # precision
    if abs(z) < eps:
        h = 0                   # hr(rho,0)=0
    else:
        if abs(1-z) < eps:
            h = 1               # hr(rho,1)=1
        else:
            I = integrate.quad(lambda v: (v**rho)/(1-z+z*v), 0, 1)
            h = I[0]*rho*z      # multiply the integral
    return h



def h1r(rho,z):
    """returns the evaluation at z in [0,1], of the derivative
    in rho of the generating function of the Yule
    distribution with parameter rho
    """
    eps = 10**(-8)              # precision
    if abs(z) < eps:
        h = 0                   # h1r(rho,0)=0
    else:
        if abs(1-z) < eps:
            h = 0               # h1r(rho,1)=0
        else:
            I = integrate.quad(lambda v: ((v**rho)*(1+rho*np.log(v)))/(1-z+z*v), 0, 1) # integrate from 0 to 1
            h = I[0]*z                               # multiply the integral
    return h



def hrSolve(y, z1, z2):
    """returns the value of rho such that
    (1-hr(rho,z1))/(1-hr(rho,z2))=y
    the equation is solved by Newton's method
    """
    def f(r):
        return ((1-hr(r,z1))/(1-hr(r,z2)))-y # equation to solve: f(r)=0
    def df(r):
        return ((1-hr(r,z1))*h1r(r,z2)-(1-hr(r,z2))*h1r(r,z1))/(1-hr(r,z2))**2 # derivative of f
    R = list(np.arange(0.3,3.01,0.01))   # vector of rho values
    nr = len(R)                          # number of rho values
    rb = 1                               # initialization
    best = f(rb)                         # initialization
    for i in range(0,nr):                # grid search
        r = R[i]                         # value in the grid
        va = abs(f(r))                   # its score
        if va < best:                    # if better than previous
            best = va                    # update best score
            rb = r                       # update best value
    eps = 10**(-5)                       # precision
    ro = rb                              # initialize by grid search result
    it = 1                               # number of iterations
    crit = True                          # stopping criterion
    while crit:
        rn = ro-f(ro)/df(ro)             # next value
        crit = (abs(rn-ro) > eps)        # check precision
        it = it+1                        # count iteratio
        if it > 10:
            crit = False                 # convergence failure
            rn = rb                      # default grid search result
        else:
            ro = rn
    return rn


def gar(alpha,rho,z):
    """
    returns the evaluation at z of the generating function
    of the LD(alpha,rho)
    """
    g = np.e**(alpha*(hr(rho,z)-1))  # from the gf of the Yule distribution
    return g


def covg(alpha,rho,z1,z2):
    """#   returns the covariance of (z1^X,z2^X), where X
    #   follows the LD(alpha,rho)
    #"""
    v = gar(alpha,rho,z1*z2) - gar(alpha,rho,z1)*gar(alpha,rho,z2)
    return v


def GFCovar(alpha,rho,z1,z2,z3):
    """#   returns the covariance matrix of the GF estimator of
    #   (alpha,rho) for the LD(alpha,rho)
    #"""
    c11 = covg(alpha,rho,z1,z1)           # variance of z1^X
    c22 = covg(alpha,rho,z2,z2)           # variance of z2^X
    c33 = covg(alpha,rho,z3,z3)           # variance of z3^X
    c12 = covg(alpha,rho,z1,z2)           # covariance (z1^X,z2^X)
    c13 = covg(alpha,rho,z1,z3)           # covariance (z1^X,z3^X)
    c23 = covg(alpha,rho,z2,z3)           # covariance (z2^X,z3^X)
    M1 = np.vstack((np.array((c11,c12,c13)), np.array((c12,c22,c23)), np.array((c13,c23,c33))))          # covariance matrix
    D = alpha*((hr(rho,z2)-1)*h1r(rho,z1) - (hr(rho,z1)-1)*h1r(rho,z2))    # common denominator
    R1 = (hr(rho,z2)-1)/(D*gar(alpha,rho,z1))         # coefficient of first coordinate
    R2 = -(hr(rho,z1)-1)/(D*gar(alpha,rho,z2))        # coefficient of second coordinate
    R3 = 0                                            # coefficient of third coordinate
    A1 = (alpha*h1r(rho,z3)*R1)/(1-hr(rho,z3))  # coefficient of first coordinate
    A2 = (alpha*h1r(rho,z3)*R2)/(1-hr(rho,z3))  # coefficient of second coordinate
    A3 = 1/(gar(alpha,rho,z3)*(hr(rho,z3)-1))   # coefficient of third coordinate
    CO = np.vstack((np.array((A1,R1)),np.array((A2,R2)),np.array((A3,R3)))) # matrix of coefficients
    V = np.dot(np.dot(np.transpose(CO),M1), CO)    # covariance matrix
    return V


#--------------- Generating Function Estimates -------------------

def tuning():
    """ returns the tuning constants for the GF estimators """
    tu = {'z1':0.1,'z2':0.9,'z3':0.8,'q':0.1}
    return tu




def GFEstAlpha(S,rho,b):
    """#   returns the GF estimate of alpha for a sample S of the
    #   LD(alpha,rho), given rho and a scaling factor b
    #  """
    tu = tuning()                         # tuning constants
    z3 = tu['z3']                         # fixed value
    g = np.mean(z3**(S/b))                # empirical generating function
    a = np.log(g)/(hr(rho,z3**(1/b))-1)   # estimate of alpha
    return a


def GFEstRho(S,b):
    """#   returns the GF estimate of rho for a sample S of the
    #   LD(alpha,rho), given a scaling factor b
    #  """
    tu = tuning()                        # tuning constants
    z1 = tu['z1']                        # lower value
    z2 = tu['z2']                        # higher value
    g1 = np.mean(z1**(S/b))              # empirical generating function at z1
    g2 = np.mean(z2**(S/b))              # empirical generating function at z2
    y = np.log(g1)/np.log(g2)            # get ratio of logs
    r = hrSolve(y,z1**(1/b),z2**(1/b))   # find corresponding r
    return r


def GFEst(S):
    """#   returns the GF estimates of alpha and rho
    #   for a sample S of the LD(alpha,rho)
    #  """
    tu = tuning()                       # tuning constants
    q = tu['q']                         # quantile
    b = np.percentile(S,q*100)+1        # scaling factor
    r = GFEstRho(S,b)                   # first estimate of rho
    a = GFEstAlpha(S,r,b)               # first estimate of alpha
    return np.hstack((a,r))


def GFConfint(S,level=0.95):
    """#   returns confidence intervals on the GF estimates of
    #   alpha and rho for a sample S of the LD(alpha,rho)
    #"""
    tu = tuning()                          # tuning constants
    z1 = tu['z1']                          # fixed value
    z2 = tu['z2']                          # fixed value
    z3 = tu['z3']                          # fixed value
    q = tu['q']                            # quantile
    b = np.percentile(S,q*100)+1           # scaling factor
    ar = GFEst(S)                          # point estimates
    a, r = ar[0], ar[1]                    # estimates of alpha and rho
    z1 = z1**(1/b)                         # rescale z1
    z2 = z2**(1/b)                         # rescale z2
    z3 = z3**(1/b)                         # rescale z3
    M = GFCovar(a,r,z1,z2,z3)              # asymptotic covariance matrix
    sig = np.sqrt(np.diag(M))              # asymptotic standard deviations
    p = (1-level)/2                        # lower probability
    ma = norm.ppf(1-p)*sig/np.sqrt(len(S)) # quantile of standard normal *  margin of error
    M = np.transpose(np.vstack((ar-ma, ar+ma)))         # matrix of two intervals
    M = { 'alpha': {str(p*100)+"%": M[0][0], str((1-p)*100)+"%": M[0][1]}, \
          'rho' : {str(p*100)+"%": M[1][0], str((1-p)*100)+"%": M[1][1]}} #rho: lower, upper
    return M


def GFTestAlpha(S,a0):
    """#   returns the p-value of the one-sided test of alpha=a0,
    #   using the GF estimate of alpha for a sample S of the
    #   LD(alpha,rho), rho being unknown
    #  """
    tu = tuning()                         # tuning constants
    z1 = tu['z1']                         # fixed value
    z2 = tu['z2']                         # fixed value
    z3 = tu['z3']                         # fixed value
    q = tu['q']                           # quantile
    b = np.percentile(S,q*100)+1          # scaling factor
    estim = GFEst(S)                      # point estimates
    a, r = estim[0], estim[1] # estimates of alpha and rho
    z1 = z1**(1/b)                        # rescale z1
    z2 = z2**(1/b)                        # rescale z2
    z3 = z3**(1/b)                        # rescale z3
    M = GFCovar(a,r,z1,z2,z3)             # asymptotic covariance matrix
    V = M[0][0]                           # asymptotic variance of alpha
    T = a-a0                              # difference with null hypothesis
    T = T/np.sqrt(V)                      # standard deviation
    T = T*np.sqrt(len(S))                 # test statistic
    if T < 0:
        p = norm.cdf(T)                   # left tail
    else:
        p = 1 - norm.cdf(T)               # right tail
    return p



def GFTestRho(S,r0):
    """#   returns the p-value of the one-sided test of rho=r0,
    #   using the GF estimate of rho for a sample S of the
    #   LD(alpha,rho)
    #   """
    tu = tuning()                         # tuning constants
    z1 = tu['z1']                         # fixed value
    z2 = tu['z2']                         # fixed value
    z3 = tu['z3']                         # fixed value
    q = tu['q']                           # quantile
    b = np.percentile(S,q*100)+1          # scaling factor
    estim = GFEst(S)                      # point estimates
    a, r = estim[0], estim[1]             # estimates of alpha and rho
    z1 = z1**(1/b)                        # rescale z1
    z2 = z2**(1/b)                        # rescale z2
    z3 = z3**(1/b)                        # rescale z3
    M = GFCovar(a,r,z1,z2,z3)             # asymptotic covariance matrix
    V = M[1][1]                           # asymptotic variance of rho
    T = r-r0                              # difference with null hypothesis
    T = T/np.sqrt(V)                      # standard deviation
    T = T*np.sqrt(len(S))                 # test statistic
    if T < 0:
        p = norm.cdf(T)                   # left tail
    else:
        p = 1 - norm.cdf(T)               # right tail
    return p


#______________________________________________________________________________

def contact(request):
    if request.method == 'POST': # If the form has been submitted...
        mut_form = MutForm(request.POST) # A form bound to the POST data

        if mut_form.is_valid(): # All validation rules pass
            z = mut_form.cleaned_data['z']
            b = mut_form.cleaned_data['b']
            N0 = mut_form.cleaned_data['N0']

            fluctuation = mut_form.cleaned_data['fluctuation']
            #fluctuationList = list(chunks(re.findall(r"[\w']+", fluctuation),2))
            test = re.findall(r"\w\.?\w*", fluctuation)
            fluctuationList = list(chunks(re.findall(r"[\w']+\.?[\w']*", fluctuation),2))
            fluctuationList = list(chunks(re.findall(r"[\w]+\.?[\w]*", fluctuation),2))
            #f = [[int(x[0]), int(float(x[1]))] for x in fluctuationList]
            f = [[int(x[0]), int(x[1].split(".")[0])] for x in fluctuationList]
            z = float(z)
            N0 = int(N0)
            Nmut = [x[0] for x in f]
            Np = [x[1] for x in f]

            if not b:
                res = GF(Nmut, N0, Np, z)
            else:
                b = float(b)
                res = CalculateAllLCFDiff(Nmut , N0, Np, b, z)

            return render_to_response('contact.html', { 'mut_form': mut_form, 'res': res}, context_instance=RequestContext(request))


    else: #if it is a GET
        mut_form = MutForm() # An unbound form
        res = 'vide'
    return render_to_response('contact.html', { 'mut_form': mut_form, 'res': res}, context_instance=RequestContext(request))



