import networkx as nx
from math import fabs, ceil
import random as rd
from scipy.stats import beta as beta_dist


def CNC(H: nx.Graph, u: int, v: int) -> int:
    return len(set(H[u]) & set(H[v]))


def PIN(G: nx.Graph) -> nx.Graph:
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    for u, v in G.edges():
        if rd.random() <= G[u][v]['p']:
            H.add_edge(u, v)
    return H


def _change_supdis(G, H, Hsupdis, u, v, s):
    for w in set(G[u]) & set(G[v]):
        if (v, w) in H.edges():
            Hsupdis[(min(u, w), max(u, w))] += s
        if (u, w) in H.edges():
            Hsupdis[(min(v, w), max(v, w))] += s


def P_RSR(G: nx.Graph, G0: nx.Graph, k: int, cnexp: dict) -> tuple:
    H = nx.Graph(G0)
    for u, v in G.edges():
        G[u][v]['supexp'] = cnexp[(min(u, v), max(u, v))]
    Hsupdis = {(min(u, v), max(u, v)): CNC(H, u, v) - G[u][v]['supexp']
               for u, v in G.edges()}
    E_lst = list(G.edges())
    for _ in range(k):
        u, v = rd.sample(E_lst, 1)[0]
        s = -1 if (u, v) in H.edges() else 1
        delta = 0.0
        for w in set(G[u]) & set(G[v]):
            if (u, w) in H.edges():
                delta += fabs(Hsupdis[(min(v, w), max(v, w))] + s) - fabs(Hsupdis[(min(v, w), max(v, w))])
            if (v, w) in H.edges():
                delta += fabs(Hsupdis[(min(u, w), max(u, w))] + s) - fabs(Hsupdis[(min(u, w), max(u, w))])
        if delta < 0:
            _change_supdis(G, H, Hsupdis, u, v, s)
            if s == -1:
                H.remove_edge(u, v)
            else:
                H.add_edge(u, v)
    return H, Hsupdis


def _init_ir(G, H, cnexp):
    for u, v in G.edges():
        G[u][v]['supexp']   = cnexp[(min(u, v), max(u, v))]
        G[u][v]['roundexp'] = round(G[u][v]['supexp'])
    nbr = {u: set(H[u]) for u in H.nodes()}
    return {(min(u, v), max(u, v)): len(nbr[u] & nbr[v]) - G[u][v]['roundexp']
            for u, v in G.edges()}


def P_IR(G: nx.Graph, G0: nx.Graph, k: int, cnexp: dict) -> tuple:
    H = nx.Graph(G0)
    Hsupdis = _init_ir(G, H, cnexp)

    E_lst = list(G.edges())
    esize = len(E_lst)
    H_nbr = {u: set(H[u]) for u in H.nodes()}

    for _ in range(k):
        u, v = E_lst[rd.randint(0, esize - 1)]
        cn_uv = set(G[u]) & set(G[v])
        gt0 = lt0 = eq0 = 0
        for w in cn_uv:
            key_vw = (min(v, w), max(v, w))
            key_uw = (min(u, w), max(u, w))
            if w in H_nbr[v]:
                d = Hsupdis[key_uw]
                if d > 0: gt0 += 1
                elif d < 0: lt0 += 1
                else: eq0 += 1
            if w in H_nbr[u]:
                d = Hsupdis[key_vw]
                if d > 0: gt0 += 1
                elif d < 0: lt0 += 1
                else: eq0 += 1

        uv_in_H = v in H_nbr[u]
        if uv_in_H and gt0 > lt0 + eq0:
            for w in cn_uv:
                key_vw = (min(v, w), max(v, w))
                key_uw = (min(u, w), max(u, w))
                if w in H_nbr[v]: Hsupdis[key_uw] -= 1
                if w in H_nbr[u]: Hsupdis[key_vw] -= 1
            H.remove_edge(u, v)
            H_nbr[u].discard(v); H_nbr[v].discard(u)
        elif not uv_in_H and lt0 > gt0 + eq0:
            for w in cn_uv:
                key_vw = (min(v, w), max(v, w))
                key_uw = (min(u, w), max(u, w))
                if w in H_nbr[v]: Hsupdis[key_uw] += 1
                if w in H_nbr[u]: Hsupdis[key_vw] += 1
            H.add_edge(u, v)
            H_nbr[u].add(v); H_nbr[v].add(u)

    return H, Hsupdis


def P_BIR(G: nx.Graph, G0: nx.Graph, th_q: float, th_c: float,
          cnexp: dict) -> tuple:
    H = nx.Graph(G0)
    Hsupdis = _init_ir(G, H, cnexp)

    E_lst = list(G.edges())
    esize = len(E_lst)
    H_nbr = {u: set(H[u]) for u in H.nodes()}

    n0 = ceil(1.0 / th_q)
    n = n0
    while beta_dist.cdf(th_q, 2, n) < th_c:
        n <<= 1

    i = 1; m = 0
    while True:
        u, v = E_lst[rd.randint(0, esize - 1)]
        cn_uv = set(G[u]) & set(G[v])
        gt0 = lt0 = eq0 = 0
        for w in cn_uv:
            key_vw = (min(v, w), max(v, w))
            key_uw = (min(u, w), max(u, w))
            if w in H_nbr[v]:
                d = Hsupdis[key_uw]
                if d > 0: gt0 += 1
                elif d < 0: lt0 += 1
                else: eq0 += 1
            if w in H_nbr[u]:
                d = Hsupdis[key_vw]
                if d > 0: gt0 += 1
                elif d < 0: lt0 += 1
                else: eq0 += 1

        flipped = False
        uv_in_H = v in H_nbr[u]
        if uv_in_H and gt0 > lt0 + eq0:
            for w in cn_uv:
                key_vw = (min(v, w), max(v, w))
                key_uw = (min(u, w), max(u, w))
                if w in H_nbr[v]: Hsupdis[key_uw] -= 1
                if w in H_nbr[u]: Hsupdis[key_vw] -= 1
            H.remove_edge(u, v)
            H_nbr[u].discard(v); H_nbr[v].discard(u)
            flipped = True
        elif not uv_in_H and lt0 > gt0 + eq0:
            for w in cn_uv:
                key_vw = (min(v, w), max(v, w))
                key_uw = (min(u, w), max(u, w))
                if w in H_nbr[v]: Hsupdis[key_uw] += 1
                if w in H_nbr[u]: Hsupdis[key_vw] += 1
            H.add_edge(u, v)
            H_nbr[u].add(v); H_nbr[v].add(u)
            flipped = True

        if flipped: m += 1
        if i % n == 0:
            if beta_dist.cdf(th_q, m + 1, n - m + 1) >= th_c:
                break
            m = 0
        i += 1

    return H, Hsupdis
