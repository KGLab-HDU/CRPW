import networkx as nx
from math import fabs
import random as rd


def _swap_remove(lst, idx_map, key):
    i = idx_map.pop(key)
    last = lst[-1]
    lst[i] = last
    lst.pop()
    if i < len(lst):
        idx_map[last] = i


def _swap_add(lst, idx_map, key):
    idx_map[key] = len(lst)
    lst.append(key)


def ADR(G: nx.Graph, steps: int) -> nx.Graph:
    exp_deg = {u: sum(G[u][v]['p'] for v in G[u]) for u in G.nodes()}

    E_lst = []
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    P = 0.0
    for u, v in G.edges():
        E_lst.append((u, v, G[u][v]['p'], rd.random(), False))
        P += G[u][v]['p']
    E_lst = sorted(E_lst, key=lambda o: o[2], reverse=True)
    P = round(P)
    n = 0
    while P != n:
        for i in range(len(E_lst)):
            u, v, p, rn, sel = E_lst[i]
            if rn <= p and not sel:
                n += 1
                E_lst[i] = (u, v, p, rn, True)
                if n == P: break
            elif rn > p and not sel:
                E_lst[i] = (u, v, p, rd.random(), False)

    E0 = [(u, v) for u, v, p, rn, sel in E_lst if sel]
    E1 = [(u, v) for u, v, p, rn, sel in E_lst if not sel]
    H.add_edges_from(E0)

    adj = {u: set(H[u]) for u in H.nodes()}
    E0_idx = {e: i for i, e in enumerate(E0)}
    E1_idx = {e: i for i, e in enumerate(E1)}
    cur_deg = {u: len(adj[u]) for u in adj}

    for k in range(steps):
        i0 = rd.randint(0, len(E0) - 1)
        i1 = rd.randint(0, len(E1) - 1)
        u, v = E0[i0]; x, y = E1[i1]

        du = cur_deg[u] - exp_deg[u]; dv = cur_deg[v] - exp_deg[v]
        dx = cur_deg[x] - exp_deg[x]; dy = cur_deg[y] - exp_deg[y]
        d1 = fabs(du-1)+fabs(dv-1)-fabs(du)-fabs(dv)
        d2 = fabs(dx+1)+fabs(dy+1)-fabs(dx)-fabs(dy)

        if d1 + d2 < 0:
            adj[u].discard(v); adj[v].discard(u)
            adj[x].add(y);     adj[y].add(x)
            cur_deg[u] -= 1; cur_deg[v] -= 1
            cur_deg[x] += 1; cur_deg[y] += 1
            _swap_remove(E0, E0_idx, (u, v)); _swap_add(E0, E0_idx, (x, y))
            _swap_remove(E1, E1_idx, (x, y)); _swap_add(E1, E1_idx, (u, v))

    H2 = nx.Graph()
    H2.add_nodes_from(H.nodes())
    for u in adj:
        for v in adj[u]:
            if u < v: H2.add_edge(u, v)
    return H2


def TRPW(G: nx.Graph, steps: int) -> nx.Graph:
    exp_deg = {u: sum(G[u][v]['p'] for v in G[u]) for u in G.nodes()}

    def exp_tri_deg(u):
        s = 0.0
        nbr_u = set(G[u])
        for v in G[u]:
            for w in nbr_u & set(G[v]):
                s += G[u][v]['p'] * G[u][w]['p'] * G[v][w]['p']
        return s / 2

    E_lst = list(G.edges())
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    P = round(sum(G[u][v]['p'] for u, v in G.edges()))
    E0 = rd.sample(E_lst, P)
    E1 = list(set(map(tuple, E_lst)) - set(map(tuple, E0)))
    H.add_edges_from(E0)

    adj = {u: set(H[u]) for u in H.nodes()}
    t_dis = {}
    for u in G.nodes():
        tri = sum(len(adj[u] & adj[v]) for v in adj[u]) / 2
        t_dis[u] = tri - exp_tri_deg(u)

    E0_idx = {e: i for i, e in enumerate(E0)}
    E1_idx = {e: i for i, e in enumerate(E1)}
    cur_deg = {u: len(adj[u]) for u in adj}

    for k in range(steps):
        i0 = rd.randint(0, len(E0) - 1)
        i1 = rd.randint(0, len(E1) - 1)
        u, v = E0[i0]; x, y = E1[i1]

        du = cur_deg[u]-exp_deg[u]; dv = cur_deg[v]-exp_deg[v]
        dx = cur_deg[x]-exp_deg[x]; dy = cur_deg[y]-exp_deg[y]
        d1 = fabs(du-1)+fabs(dv-1)-fabs(du)-fabs(dv)
        d2 = fabs(dx+1)+fabs(dy+1)-fabs(dx)-fabs(dy)

        cn_uv = adj[u] & adj[v]; cn_xy = adj[x] & adj[y]
        c1 = len(cn_uv); c2 = len(cn_xy)
        tdisu, tdisv, tdisx, tdisy = t_dis[u], t_dis[v], t_dis[x], t_dis[y]
        t1 = fabs(tdisu-c1)-fabs(tdisu) + fabs(tdisv-c1)-fabs(tdisv)
        t2 = fabs(tdisx+c2)-fabs(tdisx) + fabs(tdisy+c2)-fabs(tdisy)
        for w in cn_uv: t1 += fabs(t_dis[w]-1)-fabs(t_dis[w])
        for w in cn_xy: t2 += fabs(t_dis[w]+1)-fabs(t_dis[w])

        if d1 + d2 + t1 + t2 < 0:
            adj[u].discard(v); adj[v].discard(u)
            adj[x].add(y);     adj[y].add(x)
            cur_deg[u] -= 1; cur_deg[v] -= 1
            cur_deg[x] += 1; cur_deg[y] += 1
            _swap_remove(E0, E0_idx, (u, v)); _swap_add(E0, E0_idx, (x, y))
            _swap_remove(E1, E1_idx, (x, y)); _swap_add(E1, E1_idx, (u, v))
            t_dis[u] = tdisu-c1; t_dis[v] = tdisv-c1
            t_dis[x] = tdisx+c2; t_dis[y] = tdisy+c2
            for w in adj[u] & adj[v]: t_dis[w] -= 1
            for w in adj[x] & adj[y]: t_dis[w] += 1

        if len(E1) == 0:
            continue
        i3 = rd.randint(0, len(E1) - 1)
        u3, v3 = E1[i3]
        cn_e3 = adj[u3] & adj[v3]
        c3 = len(cn_e3)
        tdi_u3, tdi_v3 = t_dis[u3], t_dis[v3]
        d3 = (fabs(tdi_u3+c3)-fabs(tdi_u3) + fabs(tdi_v3+c3)-fabs(tdi_v3)
              + sum(fabs(t_dis[w]+1)-fabs(t_dis[w]) for w in cn_e3))
        if d3 < 0:
            adj[u3].add(v3); adj[v3].add(u3)
            cur_deg[u3] += 1; cur_deg[v3] += 1
            _swap_remove(E1, E1_idx, (u3, v3))
            _swap_add(E0, E0_idx, (u3, v3))
            t_dis[u3] = tdi_u3+c3; t_dis[v3] = tdi_v3+c3
            for w in adj[u3] & adj[v3]: t_dis[w] += 1

    H2 = nx.Graph()
    H2.add_nodes_from(H.nodes())
    for u in adj:
        for v in adj[u]:
            if u < v: H2.add_edge(u, v)
    return H2


def _ABM_bisearch(Q, x, y, wxy, weight):
    l, h = 0, len(Q) - 1
    while l <= h:
        mid = (l + h) // 2
        if   weight[(x, y)] > wxy: h = mid - 1
        elif weight[(x, y)] < wxy: l = mid + 1
        else:                      l = mid; break
    Q[l:l] = [(x, y)]


def ABM(G: nx.Graph) -> nx.Graph:
    real_expdeg = {u: sum(G[u][v]['p'] for v in G[u]) for u in G.nodes()}
    int_expdeg  = {u: round(real_expdeg[u]) for u in G.nodes()}
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    for u, v in G.edges():
        if len(H[u]) < int_expdeg[u] and len(H[v]) < int_expdeg[v]:
            H.add_edge(u, v)
    A = {u for u in G.nodes() if len(H[u]) - real_expdeg[u] <= -0.5}
    B = {u for u in G.nodes() if -0.5 < len(H[u]) - real_expdeg[u] < 0}
    weight = {}
    G_ = nx.Graph()
    Q = []
    for u, v in G.edges():
        if (u, v) in H.edges(): continue
        w = (fabs(len(H[u]) - real_expdeg[u]) + 2 * fabs(len(H[v]) - real_expdeg[v])
             - fabs(1 + len(H[u]) - real_expdeg[u]) - 1)
        if u in A and v in B and w > 0:
            weight[(u, v)] = -w
            G_.add_edge(u, v)
            Q.append((u, v))
    Q = sorted(Q, key=lambda o: weight[o])
    while Q:
        u, v = Q[0]
        Q.remove(Q[0])
        H.add_edge(u, v)
        if v in A: u, v = v, u
        for x in list(G_[v]):
            G_.remove_edge(v, x)
            if (x, v) in Q: Q.remove((x, v))
        if -1 <= len(H[u]) - real_expdeg[u] <= -0.5:
            for x in list(G_[u]):
                Q.remove((u, x))
                wux = (fabs(len(H[u]) - real_expdeg[u]) + 2 * fabs(len(H[x]) - real_expdeg[x])
                       - fabs(1 + len(H[u]) - real_expdeg[u]) - 1)
                wux = -wux
                weight[(u, x)] = wux
                _ABM_bisearch(Q, u, x, wux, weight)
        elif len(H[u]) - real_expdeg[u] > -0.5:
            for x in list(G_[u]):
                Q.remove((u, x))
    return H
