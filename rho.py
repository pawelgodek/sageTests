from sage.all import *

def order(a, p):
    i = 0
    while(not(false)):
        i += 1
        if pow(a, i, p) == 1:
            break
    return i


def pollard_rho_d(B, a, p):
    oa = order(a, p)
    row = [1, 0, 0]
    tab = [row]

    i = 1
    while(not(false)):
        row = []
        if tab[i-1][0] % 3 == 1:
            row.append((B * tab[i-1][0]) % p)
            row.append(tab[i-1][1])
            row.append((tab[i-1][2] + 1) % oa)
            tab.append(row)
        elif tab[i-1][0] % 3 == 0:
            row.append((tab[i - 1][0] * tab[i - 1][0]) % p)
            row.append((tab[i - 1][1] * 2) % oa)
            row.append((tab[i - 1][2] * 2) % oa)
            tab.append(row)
        elif tab[i-1][0] % 3 == 2:
            row.append((tab[i - 1][0] * a) % p)
            row.append((tab[i - 1][1] + 1) % oa)
            row.append(tab[i - 1][2])
            tab.append(row)

        if i % 2 == 0:
            if tab[i][0] == tab[int(i/2)][0]:
                r = (tab[int(i/2)][2] - tab[i][2]) % oa
                if r != 0:
                    x = (pow(r, -1, oa) * (tab[i][1] - tab[int(i/2)][1])) % oa
                    return x
                else:
                    return -1
        i += 1


#   Generowanie liczby pierwszej
def gen_dlp_prime(nbits):
    p_0 = random_prime(pow(2, nbits), False, pow(2, nbits-1))
    a = 2
    p = a*p_0 + 1
    while p not in Primes():
        a += 1
        p = a*p_0 + 1
    return (p, p_0, a)

#   Generowanie problemu DLP
#   Funkcja zwraca czworke (p, g, x, a), gdzie p jest liczba pierwsza, g jest generatorem F_p*,
#   x jest losowa liczba z przedzialu od 2 do p-1, taka ze a = g^x mod p
def rnd_dlp(nbits):
    p, p_0, a = gen_dlp_prime(nbits)
    g = 2
    while pow(g, a, p) == 1:
        g+=1
    g = pow(g, a, p)
    x = randrange(2, p)
    return (p, g, x, pow(g, x, p))


def dlp_set_generator(nbits):
    po = 0
    bl = 0
    print("----- a = g^x (mod p) -----")
    p, g, x, a = rnd_dlp(nbits)

    cx = pollard_rho_d(int(a), int(g), int(p))

    print("Nbits = {}\n\tp = {}\n\tg = {}\n\ta = {}\n\tx = {}".format(nbits, p,g,a,x))
    if cx != -1:
        if pow(g, x, p) == pow(g, cx, p):
            print("Zgodne cx = ", cx)
            print("\tg^{} mod p = ".format(x), pow(g, x, p))
            print("\tg^{} mod p = ".format(cx), pow(g, cx, p))
            po += 1
        else:
            print("Blad obliczen cx = ", cx)
            print("\tg^{} mod p = ".format(x), pow(g, x, p))
            print("\tg^{} mod p = ".format(cx), pow(g, cx, p))
            bl += 1
    else:
        print("Algorytm zakonczyl sie niepowodzeniem")

    if bl == 0:
        return false
    else:
        return true


def dlp_tests(nbits, t):
    bl = 0
    for i in range(t):
        bl += dlp_set_generator(nbits)
    print("Powodzenia: {}\nNiepowodzenia: {}\nRazem: {}\n".format(t-bl, bl, t))
    if bl == 0:
        return false
    else:
        return true

fal = 0
print("\n\n-------TEST1")
fal += not(dlp_tests(10, 3))
print("\n\n-------TEST2")
fal += not(dlp_tests(15, 2))
print("\n\n-------TEST3")
fal += not(dlp_tests(20, 1))

if fal == 0:
    print("Nie wszystkie testy zakonczyly sie powodzeniem")
else:
    print("Wszystkie testy zakonczone powodzeniem")