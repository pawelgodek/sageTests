from sage.all import *
from itertools import combinations
from time import time


initTime = int(round(time() * 1000))


def timeOfEvent(): #podaje czas od poczatku programy do momentu wywolania funkcji w ms
    currTime = int(round(time() * 1000))
    return " ET: " + repr(currTime - initTime) + "ms"


def getFactorBase(t, n): #wypelnianie bazy
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


def solveRov(S, n):
    print("Pary (a:b)")
    t = len(S)
    rowoiz = [int] * (t + 1), [[int] * t] * (t + 1)

    m = int(math.sqrt(n)) #m
    x = 0

    w = 0
    while w < t + 1:
        factorList = [0] * t
        a = m + x #wyliczanie parametru a

        qx = pow((x+m), 2) - n #wyliczanie wartosci funkcji q w punkcie x

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
                if tmp2 == tmp2.__int__():
                    factorList[j] = (factorList[j] + 1)
                    tmp1 = tmp2
                else:
                    break

            j = j + 1
        ###dopisywanie do zbierczej tablicy
        if tmp1 == 1: #czy liczba posiada dzielniki z poza bazy
            rowoiz[0][w] = a
            rowoiz[1][w] = factorList
            rowoiz[1][w].append(w)  # numer wiersza
            print("({}:{})".format(a, qx))
            w = w + 1

        if x < 0: #genrowanie x na podstawie poprzdniej wwartosci egz 0 1 -1 2 -2 etc
            x = (x - 1) * (-1)
        elif x > 0:
            x = (0 - x)
        else:
            x = 1

    print("Znaleziono t+1 dzielnikow q(x) " + timeOfEvent() + "\n")
    return rowoiz


def solveRowoiz(rowoiz, S, n): #wyszukiwanie zbiorow wektorow ktorych suma daje wektor zerowy
    vectors = rowoiz[1]

    numOfVectors = (len(vectors) - 1)
    lenOfVector = len(vectors[0])

    listOfZeroVectors = []

    iterVectorList = 0
    while iterVectorList <= numOfVectors:

        vectorComb = list((combinations(vectors, (iterVectorList + 1)))) #generowanie zbioru kombinacji bez powtorzen

        lenOfVectorComb = len(vectorComb)
        numOfElemInComb = len(vectorComb[0])

        iterVectorComb = 0


        while iterVectorComb < lenOfVectorComb: #iterowanie zbioru kombinacji
            iterElemInComb = 0
            testVec = [0] * (lenOfVector - 1)

            mayZeroVestors = []

            while iterElemInComb < numOfElemInComb:
                iterElemInVec = 0

                listVectorComb = list(vectorComb[iterVectorComb])

                while iterElemInVec < (lenOfVector-1): #sumowanie wektorow
                    testVec[iterElemInVec] = (testVec[iterElemInVec] + listVectorComb[iterElemInComb][iterElemInVec]).__mod__(2)
                    iterElemInVec = iterElemInVec + 1

                mayZeroVestors.append(listVectorComb[iterElemInComb][lenOfVector - 1])

                iterElemInComb = iterElemInComb + 1

            isNotVectorOfZeros = 0
            iterElemInVec = 0

            while iterElemInVec < (lenOfVector-1): #sumowanie elementow wektorow
                isNotVectorOfZeros = isNotVectorOfZeros + testVec[iterElemInVec]
                iterElemInVec = iterElemInVec + 1

            if isNotVectorOfZeros == 0: #czy suma wektorow jest rowna wektorowi zerowemu
                listOfZeroVectors.append(list(mayZeroVestors))


                factors = []

                j = 0
                cLen = len(mayZeroVestors)  # ilosc wektorow do zbadania w zbiorze
                lLen = len(rowoiz[1]) - 1

                els = [0] * lLen

                x = 1  # wyliczanie parametru x
                while j < cLen:
                    x = (x * rowoiz[0][mayZeroVestors[j]]).__mod__(n)
                    j = j + 1

                j = 0  # wyliczanie parametrow li
                while j < lLen:
                    k = 0
                    while k < cLen:
                        els[j] = els[j] + rowoiz[1][mayZeroVestors[k]][j]
                        k = k + 1
                    els[j] = int(els[j] / 2)
                    j = j + 1

                j = 0
                y = 1  # wyliczanie parametru y
                while j < lLen:
                    y = (y * S[j] ** els[j]).__mod__(n)
                    j = j + 1

                if not (((x + y) / n).is_integer() or ((x - y) / n).is_integer()):  # czy x przystaje do (+-)y modulo n
                    factor = gcd((x - y), n)  # NWD - dzielnik n

                    if not (factor in factors):
                        factors.append(factor)

                        cfactor = n / factor

                        print("Znaleziono dzielnik: " + repr(factor) + timeOfEvent() + "\n" + repr(factor) + " * " + repr(cfactor) + " = " + repr(n))
                        return 0



            iterVectorComb = iterVectorComb + 1
        iterVectorList = iterVectorList + 1


def intro():
    intro = "##############################################\n" + \
            "################ Pawel Godek #################\n" + \
            "######## Algorytm Sita Kwadratowego ##########\n" + \
            "##############################################\n"

    print(intro)

def main():
    intro()
    n = 840260382772330973
    t = 103

    print("n = {} t = {}".format(n, t))
    S = getFactorBase(t, n)
    print("Baza: ", S)
    rowoiz = solveRov(S, n)
    solveRowoiz(rowoiz, S, n)


main()
