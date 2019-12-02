# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

INFO = dict(
        CFL=0.9,
        u0=[1],
        dx=0.1,
        a=1,
        i_max=100,
        n_a=6,
        n_b=9,
        dom_x=1,
        )


def init(info, u, x):
    info["dt"] = info.get("dt") or info.get("CFL") * info["dx"] / info["a"]
    i_max = info.get("i_max")
    dx = info.get("dx")
    n_a = info.get("n_a")    
    n_b = info.get("n_b")

    for i in range(-i_max, i_max):
        x[i] = i * dx
#        x_dx = x[i]/dx
#        u[i] = np.sin(2 * np.pi * x_dx / n_a) * np.exp(-np.log(2) * (x_dx / n_b)**2 )
        u[i] = np.sin(2 * np.pi * i / n_a) * np.exp(-np.log(2) * (i / n_b)**2)

    u = condition_limites(u, i_max)

    plot_solution(u, x)
    return x, u


def condition_limites(u, i_max):
    u[0] = u[i_max - 1]
    u[i_max] = u[1]
    return u


def flux_num_roe(u, info):
    i_max = info["i_max"]
    a = info["a"]
    dx = info["dx"]    
    dt = info["dt"]
    
    f_num = np.zeros(i_max)
    a_plus, a_moins = get_a(a)
        
    for j in range(i_max):
        f_num[j] = a_moins * u[j+1] + a_plus * u[j]
        
    # for j in range(i_max):
    #     u[j] = u[j] - dt/dx * (f_num[j] - f_num[j-1])
        
    return u


def calc_u(u, f_num_roe, f_num_3, info):
    i_max = info["i_max"]
    dx = info["dx"]
    dt = info["dt"]

    f_num = np.zeros(i_max)

    for j in range(i_max):
        f_num[j] = f_num_roe[j] + f_num_3[j]
        u[j] = u[j] - dt/dx * (f_num[j] - f_num[j-1])


def flux_num_3(u, psy_2, psy_3, info):
    i_max = info["i_max"]
    a = info["a"]
    dx = info["dx"]
    dt = info["dt"]

    f_num = np.zeros(i_max)
    for j in range(i_max):
        f_num[j] = 0.5 * (psy_2[j] + psy_3[j]) if a >= 0 else 0.5 * (psy_2[j] - psy_3[j+1])


def plot_solution(u, x):
    fig = plt.figure()
    ax = plt.axes()
#    plt.xlim(-1,1)
    ax.plot(x, u)
    fig.show()
    

def get_a(a):
    a_plus = 0.5 * (a + abs(a))
    a_moins = 0.5 * (a - abs(a))
    return a_plus, a_moins
    

# def main(info):
#     i_max = info.get("i_max")
#     u = np.zeros(i_max)
#     x = np.zeros(i_max)
#     init(info, u, x)
    
    # for t in range(n_max):
    #     continue


def calc_psy(psy, wc, info):
    i_max = info.get("i_max")
    dx = info.get("dx")
    dt = info.get("dt")
    a = info.get("a")
    nu = abs(a) * dt/dx

    psy_2 = np.zeros(i_max)
    psy_3 = np.zeros(i_max)

    c2 = abs(a) * (1-nu)
    c3 = (1 + nu) * c2/3
    # c4 = (-2 - nu) * c3/4
    # c5 = (2 + nu) * c4/5

    for i in range(i_max - 1):
        psy_2[i] = c2 * (wc[i+1] - wc[i])
        psy_3[i] = c3 * (wc[i] - wc[i-1]) - c3 * (wc[i+1] - wc[i])


if __name__ == "__main__":
    u = np.zeros(INFO["i_max"]+1)
    x = np.zeros(INFO["i_max"]+1)
    x, u = init(INFO, u, x)
    u = flux_num_roe(u, INFO)
    plot_solution(u, x)
