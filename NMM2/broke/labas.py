import cmath

a = 2.0
kapa1 = 1.0
kapa2 = 1.0
gama1 = 0.0
gama2 = 0.0
N = 100
h = 0.01
tau = 0.01
T = 1.0
i=1j
a=1
d=1
c=1


# Pradzia Pirmam testui
def u_tikslus(x, t):
    return (2 - i*t**2)*(-x + x**3/3)

def f(x, t):
    u = u_tikslus(x, t)
    udt = -2/3*i*t*x*(-3 + x**2)
    uddxx = 2*(2 - i*t**2)*x
    f = udt - (a**2+i)*uddxx - i*c*u - i*d*(abs(u)**2)*u

    return f


def netiktis1(h, tau, x, t):
    uj = u_tikslus(x, t)
    uj_plius_1 = u_tikslus(x + h, t)
    uj_minus_1 = u_tikslus(x - h, t)

    uj_stog = u_tikslus(x, t + tau)
    uj_stog_plius_1 = u_tikslus(x + h, t + tau)
    uj_stog_minus_1 = u_tikslus(x - h, t + tau)

    fj = f(x, t)
    fj_stog = f(x, t + tau)

    A = (uj_stog_plius_1 - (2 * uj_stog) + uj_stog_minus_1 + uj_plius_1 - (2 * uj) + uj_minus_1) / (h ** 2)
    B = ((uj_stog+uj)/2.0)
    C = (abs((uj_stog+uj)/2.0)**2) * ((uj_stog+uj)/2.0)
    D = (fj_stog + fj) / 2

    kaire = (uj_stog - uj) / tau
    desine = ((a**2+i)*  0.5 * A) + (i*c*B) + (i*d*C) + D

    netiktis = abs(kaire - desine)

    return netiktis


def testas_1():
    x = 0.79
    t = 1.89
    h = 0.01
    tau = 0.01

    return netiktis1(h, tau, x, t) / netiktis1(h/10.0, tau / 10.0, x, t)


# Pabaiga Pirmam testui

# Pradzia Antras testas
def C(h, tau):
    Const = 2.0 + 2.0*h**2/((a**2+i)*tau)
    return Const


def F(uj, uj_stog, uj_plius_1, uj_minus_1, uj_stog_plius_1, uj_stog_minus_1, fj, fj_stog, h, tau):
    A = (uj_plius_1 - (2 * uj) + uj_minus_1) / (h ** 2)
    B = i*c*((uj_stog+uj)/2.0)
    C = i*d*(abs((uj_stog+uj)/2.0)**2) * ((uj_stog+uj)/2.0)
    D = (fj_stog + fj) / 2
    daugiklis = 2.0*h**2/(a**2+i)

    Fj = -1*daugiklis/tau*uj - A - daugiklis*B - daugiklis*i*d*C - daugiklis*D
    return Fj


def netiktis2(h, tau, x, t):
    uj = u_tikslus(x, t)
    uj_plius_1 = u_tikslus(x + h, t)
    uj_minus_1 = u_tikslus(x - h, t)

    uj_stog = u_tikslus(x, t + tau)
    uj_stog_plius_1 = u_tikslus(x + h, t + tau)
    uj_stog_minus_1 = u_tikslus(x - h, t + tau)

    fj = f(x, t)
    fj_stog = f(x, t + tau)

    netiktis_ = abs(uj_stog_plius_1 - C(h, tau) * uj_stog + uj_stog_minus_1 + F(uj, uj_stog, uj_plius_1, uj_minus_1,
                                                                                uj_stog_plius_1, uj_stog_minus_1, fj,
                                                                                fj_stog, h, tau))

    return netiktis_


def testas_2():
    x = 0.79
    t = 1.89
    h = 0.01
    tau = 0.01

    return netiktis2(h, tau, x, t) / netiktis2(h / 10.0, tau / 10.0, x, t)


# Pabaiga Antras testas

# Pradzia Treciam testui
def alpha_plius_1(alpha, h, tau):
    return (1.0 / (C(h, tau) - alpha))


def beta_plius_1(alpha_plius, beta, F_):
    return (alpha_plius * (beta + F_))


def Thomas(alpha, uj_stog_new, beta):
    uj_stog_new_minus_1 = alpha * uj_stog_new + beta
    return uj_stog_new_minus_1


def testas_3():
    alpha = [0] * (N + 1)
    beta = [0] * (N + 1)
    p = [0] * (N + 1)
    y = [0] * (N + 1)
    Fj = [0] * (N + 1)

    alpha[1] = kapa1
    beta[1] = gama1

    for i in range(1, N, 1):
        p[i] = cmath.sin(i ** 2 - 7) + cmath.cos(i * 8) * 5 * 1j

    p[0] = kapa1 * p[1] + gama1
    p[N] = kapa2 * p[N - 1] + gama2

    for i in range(1, N, 1):
        Fj[i] = C(h, tau) * p[i] - p[i + 1] - p[i - 1]
        alpha[i + 1] = alpha_plius_1(alpha[i], h, tau)
        beta[i + 1] = beta_plius_1(alpha[i + 1], beta[i], Fj[i])

    y[N] = ((kapa2 * beta[N]) + gama2) / (1.0 - (kapa2 * alpha[N]))

    for i in range(N, 0, -1):
        y[i - 1] = Thomas(alpha[i], y[i], beta[i])

    skirtumas = difference(y, p)

    # if(skirtumas < 10**(-12)):
    #  return True

    return skirtumas


# Pabaiga Treciam testui

def difference(firstArgument, secondArgument):
    maxDiff = 0.0

    for i in range(-1, N, 1):
        temp = abs(firstArgument[i] - secondArgument[i])
        if (temp > maxDiff):
            maxDiff = temp

    return maxDiff


print("1 testo rezultatas: ", testas_1())
print("2 testo rezultatas: ", testas_2())
print("3 testo rezultatas: ", testas_3())


# Pradzia Thomas algoritmo implementacija
def ThomasImpl(tau, h, N):
    t = 0.0
    nuokrypis = 0.0

    alpha = [0] * (N + 1)
    beta = [0] * (N + 1)
    uj_tikslus = [0] * (N + 1)
    uj = [0] * (N + 1)
    uj_stog_old = [0] * (N + 1)
    uj_stog_new = [0] * (N + 1)
    fj = [0] * (N + 1)
    fj_stog = [0] * (N + 1)
    Fj = [0] * (N + 1)

    alpha[1] = kapa1
    beta[1] = gama1

    for i in range(0, N + 1, 1):
        uj[i] = u_tikslus(i * h, t)
        uj_stog_old[i] = uj[i]

    while (t < T):
        progresas = 1.0
        for i in range(0, N + 1, 1):
            fj[i] = f(i * h, t)
            fj_stog[i] = f(i * h, t + tau)

        while (progresas > 10 ** (-8)):
            for i in range(1, N, 1):
                Fj[i] = F(uj[i], uj_stog_old[i], uj[i + 1], uj[i - 1], uj_stog_old[i + 1], uj_stog_old[i - 1], fj[i],
                          fj_stog[i], h, tau)
                alpha[i + 1] = alpha_plius_1(alpha[i], h, tau)
                beta[i + 1] = beta_plius_1(alpha[i + 1], beta[i], Fj[i])

            uj_stog_new[N] = ((kapa2 * beta[N]) + gama2) / (1.0 - (kapa2 * alpha[N]))

            for i in range(N, 0, -1):
                uj_stog_new[i - 1] = Thomas(alpha[i], uj_stog_new[i], beta[i])

            progresas = difference(uj_stog_new, uj_stog_old)

            for i in range(0, N + 1, 1):
                uj_stog_old[i] = uj_stog_new[i]

        for i in range(0, N + 1, 1):
            uj[i] = uj_stog_old[i]

        t = t + tau

        for i in range(0, N + 1, 1):
            uj_tikslus[i] = u_tikslus(i * h, t)

        temp = difference(uj_stog_new, uj_tikslus)
        if (temp > nuokrypis):
            nuokrypis = temp

    return nuokrypis


def ThomasImpl_tikrinimas():
    return ThomasImpl(tau, h, N) / ThomasImpl(tau / 10.0, h / 10.0, N * 10)


# Pabaiga Thomas algoritmo implementacija

print("Thomas implementacijos patikrinimas: ", ThomasImpl_tikrinimas())

