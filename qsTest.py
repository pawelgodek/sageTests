from sage.all import *
from itertools import chain, combinations
from time import time

initTime = int(round(time() * 1000))


def timeOfEvent():
    currTime = int(round(time() * 1000))
    return " ET: " + repr(currTime - initTime) + "ms"


def poproz(x):
    if x <= 0:
        x = (0 - x) + 1
    else:
        x = 0 - x
    return x


def isInteger(x):
    if x == x.__int__():
        return true
    else:
        return false


def izZeroz(v):
    t = len(v)
    i = 0

    zero = 0
    while i < t:
        zero = zero + (v[i]).__mod__(2)
        i = i + 1
    if zero == 0:
        return true
    else:
        return false


def powerset(v):
    s = list(v)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def getFactorBase(t, n): #wypelnianie bazy dzielnikow
    S = [-1]
    P = Primes()

    nOt = 0
    i = 0

    print("Wypelnianie bazy (p1-p" + repr(t) + ")", end='')

    while nOt < t-1:
        mP = P.unrank(i)
        if kronecker(n, mP) == 1: #symbol legendre'a
            S.append(mP)
            nOt = nOt + 1
        i = i + 1

    print("..OK" + timeOfEvent())
    return S


def solveRov(S, x, n):
    t = len(S)
    factorList = [0] * t

    m = sqrt(n).__int__() #m
    a = m + x #wyliczanie parametru a

    qx = (x+m)**2 - n #wyliczanie wartosci funkcji q w punkcie x

    if qx < 0: #czy wartosc q(x) jest ujemna
        factorList[0] = 1
        qx = (S[0]) * qx

    j = 1
    #rozklad q(x) na dzielniki z bazy
    tmp1 = qx

    while j < t:
        tmp2 = tmp1
        while true:
            tmp2 = tmp2 / S[j]
            if isInteger(tmp2):
                factorList[j] = (factorList[j] + 1)
                tmp1 = tmp2
            else:
                break

        j = j + 1


    if tmp1 == 1: #czy liczba posiada dzielniki z poza bazy
        smooth = true
    else:
        smooth = false

    return smooth, a, factorList


def rowoization(S, n):
    t = len(S)
    rowoiz = [int] * (t+1), [[int] * t] * (t+1)

    smooth, a, factorList = solveRov(S, 0, n)

    i = 1
    j = 0
    if smooth == true:
        rowoiz[0][j] = a
        rowoiz[1][j] = factorList
        rowoiz[1][j].append(j)
        j = j + 1

    while j < t + 1:
        smooth, a, factorList = solveRov(S, i, n) #skladanie uzyskanych wynikow w macierz

        if smooth == true:
            rowoiz[0][j] = a
            rowoiz[1][j] = factorList
            rowoiz[1][j].append(j) #numer wiersza
            j = j + 1
        i = poproz(i)

    print("Znaleziono dzielniki q(x) " + timeOfEvent() + "\n")
    return rowoiz


def solveRowoiz(rowoiz): #wyszukiwanie zbiorow wektorow ktorych suma daje wektor zerowy
    cRowo = list(powerset(rowoiz[1])) #generowanie wszystkich mozliwych podzbiorow bez powtorzen
    lRowo = len(cRowo) #dalej sprawdzanie czy suma wektorow w podzbiorze jest wektorem zerowym

    t = len(rowoiz[1][0])
    i = 1

    vector = []
    vectors = []

    testvec = [0] * t

    while i < lRowo:
        llRowo = len(cRowo[i])
        j = 0

        while j < llRowo:
            k = 0
            while k < (t - 1):
                testvec[k] = testvec[k] + cRowo[i][j][k]
                k = k + 1
            vector.append(cRowo[i][j][t - 1])
            j = j + 1

        if izZeroz(testvec):
            vectors.append(vector) #lista wektorow ktorych suma jest wektorem zerowym

        vector = []
        testvec = [0] * t
        i = i + 1

    return vectors


def comradesVerification(vectors, rowoiz, n, S):
    i = 0
    vLen = len(vectors) #ilosc zbiorow wektorow do zbadania
    lLen = len(rowoiz[1]) - 1 #rzad wektora

    factors = [] #dzielniki

    while i < vLen:
        j = 0
        cLen = len(vectors[i]) #ilosc wektorow do zbadania w zbiorze

        els = [0] * lLen

        x = 1 #wyliczanie parametru x
        while j < cLen:
            x = (x * rowoiz[0][vectors[i][j]]).__mod__(n)
            j = j + 1

        j = 0 #wyliczanie parametrow li
        while j < lLen:
            k = 0
            while k < cLen:
                els[j] = els[j] + rowoiz[1][vectors[i][k]][j]
                k = k + 1
            els[j] = int(els[j] / 2)
            j = j + 1

        j = 0
        y = 1 #wyliczanie parametru y
        while j < lLen:
            y = (y * S[j]**els[j]).__mod__(n)
            j = j + 1


        if not((( x + y ) / n).is_integer() or (( x - y ) / n).is_integer()): #czy x przystaje do (+-)y modulo n
            factor = gcd((x - y), n) #NWD - dzielnik n

            if not(factor in factors):
                factors.append(factor)

                cfactor = n / factor

                print("Znaleziono dzielnik: " + repr(factor) + timeOfEvent() + "\n" + repr(factor) + " * " + repr(cfactor) + " = " + repr(n))
        i = i + 1

def intro():
    intro = "##############################################\n" + \
            "################ Pawel Godek #################\n" + \
            "######## Algorytm Sita Kwadratowego ##########\n" + \
            "##############################################\n"

    print(intro)

def main():
    intro()
    n = 24961
    t = 6
    S = getFactorBase(t, n)
    rowoiz = rowoization(S, n)
    vectors = solveRowoiz(rowoiz)
    comradesVerification(vectors, rowoiz, n, S)


main()
