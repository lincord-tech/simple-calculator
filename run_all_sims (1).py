"""
Sigma Model — все симуляции для проверки закона d = 2j + 1
===========================================================

Запуск: python run_all_sims.py

Результаты печатаются в конце в виде таблицы — скопируй и отправь.

Что проверяем:
  Sim 1: j=1/2 на 3D решётке без агрессивной заморозки (Rule B)
          Ожидание: d ≈ 4 (если закон d_lat+1 работает)

  Sim 2: j=3/2 на 2D решётке (повтор с лучшей динамикой)
          Ожидание: d ≈ 3-4

  Sim 3: j=3/2 на 3D решётке
          Ожидание: d ≈ 4-5

  Sim 4: j=2 на 2D решётке
          Ожидание: d ≈ 5 (если закон d=2j+1)

  Sim 5: j=1 на 3D решётке
          Ожидание: d ≈ 4

Каждая симуляция независима. Можно запускать по одной если нужно.
Просто закомментируй ненужные вызовы в конце файла.
"""

"""

# ================================================================
# FAST CAUSAL GRAPH BUILDER (same logic, much faster)
# Replace ONLY your old build_causal_graph(events) with this version
# ================================================================

from collections import defaultdict

def build_causal_graph(events):

    Same exact logic as original:
    connect events only between t and t+1
    if Manhattan distance <= 1

    But instead of O(N^2), uses coordinate lookup.
    

    graph = defaultdict(set)
    by_t = defaultdict(list)

    # group by time
    for idx, ev in enumerate(events):
        coords = ev[:-1]
        t = ev[-1]
        by_t[t].append((idx, coords))

    times = sorted(by_t.keys())

    for t in times:
        if (t + 1) not in by_t:
            continue

        curr = by_t[t]
        nxt = by_t[t + 1]

        # hash all next-time events by coordinates
        pos_map = defaultdict(list)
        for idx2, coords2 in nxt:
            pos_map[coords2].append(idx2)

        ndim = len(curr[0][1])

        for idx1, coords1 in curr:

            # generate all points with Manhattan distance <=1
            candidates = [coords1]   # same position

            for d in range(ndim):
                c1 = list(coords1)
                c1[d] += 1
                candidates.append(tuple(c1))

                c2 = list(coords1)
                c2[d] -= 1
                candidates.append(tuple(c2))

            for cand in candidates:
                if cand in pos_map:
                    for idx2 in pos_map[cand]:
                        graph[idx1].add(idx2)
                        graph[idx2].add(idx1)

    return graph



"""




import numpy as np
import math
from collections import defaultdict, deque

# ================================================================
# ОБЩИЕ УТИЛИТЫ
# ================================================================







from collections import defaultdict

def build_causal_graph(events):
    graph = defaultdict(set)
    by_t = defaultdict(list)

    # group by time
    for idx, ev in enumerate(events):
        coords = ev[:-1]
        t = ev[-1]
        by_t[t].append((idx, coords))

    times = sorted(by_t.keys())

    for t in times:
        if (t + 1) not in by_t:
            continue

        curr = by_t[t]
        nxt = by_t[t + 1]

        # hash all next-time events by coordinates
        pos_map = defaultdict(list)
        for idx2, coords2 in nxt:
            pos_map[coords2].append(idx2)

        ndim = len(curr[0][1])

        for idx1, coords1 in curr:

            # generate all points with Manhattan distance <=1
            candidates = [coords1]   # same position

            for d in range(ndim):
                c1 = list(coords1)
                c1[d] += 1
                candidates.append(tuple(c1))

                c2 = list(coords1)
                c2[d] -= 1
                candidates.append(tuple(c2))

            for cand in candidates:
                if cand in pos_map:
                    for idx2 in pos_map[cand]:
                        graph[idx1].add(idx2)
                        graph[idx2].add(idx1)

    return graph









"""

def build_causal_graph(events):
    events: list of (i, j, k, t) или (i, j, t)
    graph = defaultdict(set)
    by_t = defaultdict(list)
    for idx, ev in enumerate(events):
        by_t[ev[-1]].append((idx,) + ev[:-1])

    for t in sorted(by_t)[:-1]:
        if t + 1 not in by_t:
            continue
        for rec1 in by_t[t]:
            for rec2 in by_t[t + 1]:
                idx1, coords1 = rec1[0], rec1[1:]
                idx2, coords2 = rec2[0], rec2[1:]
                dist = sum(abs(a - b) for a, b in zip(coords1, coords2))
                if dist <= 1:
                    graph[idx1].add(idx2)
                    graph[idx2].add(idx1)
    return graph
"""

def measure_dimension(graph, n_events, label=""):
    if n_events < 30:
        return None, 0, 0

    # BFS от узла 0
    root = 0
    dist = {root: 0}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        for nb in graph[node]:
            if nb not in dist:
                dist[nb] = dist[node] + 1
                queue.append(nb)

    max_r = max(dist.values()) if dist else 0
    n_reached = len(dist)

    if max_r < 5:
        return None, max_r, n_reached

    radii  = list(range(1, max_r + 1))
    volumes = [sum(1 for d in dist.values() if d <= r) for r in radii]

    # наклон в лог-лог на средней части
    n = len(radii)
    lo, hi = n // 4, 3 * n // 4
    if hi <= lo + 2:
        return None, max_r, n_reached

    lr = [math.log(radii[i])   for i in range(lo, hi)]
    lv = [math.log(volumes[i]) for i in range(lo, hi) if volumes[i] > 0]
    if len(lr) != len(lv):
        lr = lr[:len(lv)]
    if len(lr) < 2:
        return None, max_r, n_reached

    coeffs = np.polyfit(lr, lv, 1)
    return round(coeffs[0], 3), max_r, n_reached


# ================================================================
# ДИНАМИКА: j = 1/2
# ================================================================

def vac_j05_2d(nx, ny):
    g = np.zeros((nx, ny), dtype=int)
    for i in range(nx):
        for j in range(ny):
            g[i, j] = 1 if (i + j) % 2 == 0 else -1
    return g

def vac_j05_3d(nx, ny, nz):
    g = np.zeros((nx, ny, nz), dtype=int)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                g[i, j, k] = 1 if (i + j + k) % 2 == 0 else -1
    return g

def defects_j05_3d(g):
    nx, ny, nz = g.shape
    out = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for di,dj,dk in [(1,0,0),(0,1,0),(0,0,1)]:
                    ni,nj,nk = (i+di)%nx,(j+dj)%ny,(k+dk)%nz
                    if g[i,j,k] * g[ni,nj,nk] == 1:
                        out.append((i,j,k,ni,nj,nk))
    return out

def step_j05_3d(g, frozen):
    nx,ny,nz = g.shape
    ng = g.copy()
    nf = {pos: t-1 for pos,t in frozen.items() if t > 1}

    for i,j,k,ni,nj,nk in defects_j05_3d(g):
        for ri,rj,rk in [(i,j,k),(ni,nj,nk)]:
            if (ri,rj,rk) in frozen or (ri,rj,rk) in nf:
                continue
            ng[ri,rj,rk] = -g[ri,rj,rk]
            nf[(ri,rj,rk)] = 2
            break  # Rule B: только один из пары

    return ng, nf


# ================================================================
# ДИНАМИКА: j = 1
# ================================================================

VALS_J1 = [-1, 0, 1]
V2I_J1  = {v: k for k, v in enumerate(VALS_J1)}

def vac_j1_2d(nx, ny):
    g = np.zeros((nx, ny), dtype=int)
    for i in range(nx):
        for j in range(ny):
            g[i, j] = VALS_J1[(i + j) % 3]
    return g

def vac_j1_3d(nx, ny, nz):
    g = np.zeros((nx, ny, nz), dtype=int)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                g[i, j, k] = VALS_J1[(i + j + k) % 3]
    return g

def defects_j1_3d(g):
    nx,ny,nz = g.shape
    out = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for di,dj,dk in [(1,0,0),(0,1,0),(0,0,1)]:
                    ni,nj,nk=(i+di)%nx,(j+dj)%ny,(k+dk)%nz
                    i1=V2I_J1.get(g[i,j,k],-1)
                    i2=V2I_J1.get(g[ni,nj,nk],-1)
                    if i1<0 or i2<0 or (i2-i1)%3 not in (1,2):
                        out.append((i,j,k,ni,nj,nk))
    return out

def step_j1_3d(g, frozen):
    nx,ny,nz = g.shape
    ng = g.copy()
    nf = {pos: t-1 for pos,t in frozen.items() if t > 1}

    for i,j,k,ni,nj,nk in defects_j1_3d(g):
        for ri,rj,rk in [(i,j,k),(ni,nj,nk)]:
            if (ri,rj,rk) in frozen or (ri,rj,rk) in nf:
                continue
            ci = V2I_J1.get(g[ri,rj,rk], 1)
            ni2 = min(ci+1, 2) if ci < 1 else max(ci-1, 0)
            ng[ri,rj,rk] = VALS_J1[ni2]
            nf[(ri,rj,rk)] = 2
            break
    return ng, nf


# ================================================================
# ДИНАМИКА: j = 3/2
# ================================================================

VALS_J32 = [-3, -1, 1, 3]
V2I_J32  = {v: k for k, v in enumerate(VALS_J32)}

def vac_j32_2d(nx, ny):
    g = np.zeros((nx, ny), dtype=int)
    for i in range(nx):
        for j in range(ny):
            g[i, j] = VALS_J32[(i + j) % 4]
    return g

def vac_j32_3d(nx, ny, nz):
    g = np.zeros((nx, ny, nz), dtype=int)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                g[i, j, k] = VALS_J32[(i + j + k) % 4]
    return g

def defects_j32(g, ndim):
    shape = g.shape
    dirs = [(1,0),(0,1)] if ndim==2 else [(1,0,0),(0,1,0),(0,0,1)]
    out = []
    it = np.ndindex(*shape)
    for idx in it:
        for step in dirs:
            nidx = tuple((idx[d]+step[d])%shape[d] for d in range(ndim))
            i1 = V2I_J32.get(g[idx], -1)
            i2 = V2I_J32.get(g[nidx], -1)
            if i1<0 or i2<0 or (i2-i1)%4 not in (1,3):
                out.append(idx + nidx)
    return out

def step_j32_nd(g, frozen, ndim):
    shape = g.shape
    ng = g.copy()
    nf = {pos: t-1 for pos,t in frozen.items() if t > 1}

    for pair in defects_j32(g, ndim):
        n = len(pair) // 2
        p1 = pair[:n]
        p2 = pair[n:]
        for pos in [p1, p2]:
            if pos in frozen or pos in nf:
                continue
            ci = V2I_J32.get(g[pos], 1)
            ni2 = min(ci+1, 3) if ci <= 1 else max(ci-1, 0)
            ng[pos] = VALS_J32[ni2]
            nf[pos] = 2
            break
    return ng, nf


# ================================================================
# ДИНАМИКА: j = 2
# ================================================================

VALS_J2 = [-4, -2, 0, 2, 4]
V2I_J2  = {v: k for k, v in enumerate(VALS_J2)}

def vac_j2_2d(nx, ny):
    g = np.zeros((nx, ny), dtype=int)
    for i in range(nx):
        for j in range(ny):
            g[i, j] = VALS_J2[(i + j) % 5]
    return g

def defects_j2_2d(g):
    nx,ny = g.shape
    out = []
    for i in range(nx):
        for j in range(ny):
            for di,dj in [(1,0),(0,1)]:
                ni,nj=(i+di)%nx,(j+dj)%ny
                i1=V2I_J2.get(g[i,j],-1)
                i2=V2I_J2.get(g[ni,nj],-1)
                if i1<0 or i2<0 or (i2-i1)%5 not in (1,4):
                    out.append((i,j,ni,nj))
    return out

def step_j2_2d(g, frozen):
    nx,ny=g.shape
    ng=g.copy()
    nf={pos:t-1 for pos,t in frozen.items() if t>1}

    for i,j,ni,nj in defects_j2_2d(g):
        for ri,rj in [(i,j),(ni,nj)]:
            if (ri,rj) in frozen or (ri,rj) in nf:
                continue
            ci=V2I_J2.get(g[ri,rj],2)
            ni2=min(ci+1,4) if ci<2 else max(ci-1,0)
            ng[ri,rj]=VALS_J2[ni2]
            nf[(ri,rj)]=2
            break
    return ng,nf


# ================================================================
# УНИВЕРСАЛЬНЫЙ RUNNER
# ================================================================

def run_sim(label, nx, ny, nz, j_value, n_defects, t_steps, rng_seed=42):
    print(f"\n{'='*60}")
    print(f"SIM: {label}")
    print(f"  Grid: {nx}x{ny}" + (f"x{nz}" if nz else "") +
          f"  j={j_value}  defects={n_defects}  steps={t_steps}")
    print('='*60)

    rng = np.random.default_rng(rng_seed)
    ndim = 3 if nz else 2

    # инициализация вакуума
    if j_value == 0.5:
        g = vac_j05_3d(nx,ny,nz) if nz else vac_j05_2d(nx,ny)
        vals = [-1, 1]
    elif j_value == 1.0:
        g = vac_j1_3d(nx,ny,nz) if nz else vac_j1_2d(nx,ny)
        vals = VALS_J1
    elif j_value == 1.5:
        g = vac_j32_3d(nx,ny,nz) if nz else vac_j32_2d(nx,ny)
        vals = VALS_J32
    elif j_value == 2.0:
        g = vac_j2_2d(nx,ny)
        vals = VALS_J2
    else:
        print(f"  ОШИБКА: j={j_value} не поддерживается")
        return None

    # случайные дефекты
    for _ in range(n_defects):
        if nz:
            i,j,k = rng.integers(1,nx-1), rng.integers(1,ny-1), rng.integers(1,nz-1)
            g[i,j,k] = rng.choice(vals)
        else:
            i,j = rng.integers(1,nx-1), rng.integers(1,ny-1)
            g[i,j] = rng.choice(vals)

    frozen = {}
    events = []

    for t in range(t_steps):
        # шаг динамики
        if j_value == 0.5 and nz:
            g, frozen = step_j05_3d(g, frozen)
            defs = defects_j05_3d(g)
            for d in defs: events.append(d[:3] + (t,))
        elif j_value == 1.0 and nz:
            g, frozen = step_j1_3d(g, frozen)
            defs = defects_j1_3d(g)
            for d in defs: events.append(d[:3] + (t,))
        elif j_value == 1.5:
            g, frozen = step_j32_nd(g, frozen, ndim)
            defs = defects_j32(g, ndim)
            for d in defs:
                coords = d[:ndim]
                events.append(coords + (t,))
        elif j_value == 2.0:
            g, frozen = step_j2_2d(g, frozen)
            defs = defects_j2_2d(g)
            for d in defs: events.append((d[0],d[1],t))
        elif j_value == 0.5 and not nz:
            # 2D j=1/2 — для контроля
            from collections import defaultdict as dd2
            g2 = g.copy()
            nf2 = {pos:t2-1 for pos,t2 in frozen.items() if t2>1}
            nx2,ny2=g.shape
            defs2=[]
            for ii in range(nx2):
                for jj in range(ny2):
                    for di,dj in [(1,0),(0,1)]:
                        ni2,nj2=(ii+di)%nx2,(jj+dj)%ny2
                        if g[ii,jj]*g[ni2,nj2]==1:
                            defs2.append((ii,jj,ni2,nj2))
            for ii,jj,ni2,nj2 in defs2:
                for ri,rj in [(ii,jj),(ni2,nj2)]:
                    if (ri,rj) in frozen or (ri,rj) in nf2: continue
                    g2[ri,rj]=-g[ri,rj]; nf2[(ri,rj)]=2; break
            g,frozen=g2,nf2
            defs=defs2
            for d in defs: events.append((d[0],d[1],t))
        elif j_value == 1.0 and not nz:
            g2=g.copy()
            nf2={pos:t2-1 for pos,t2 in frozen.items() if t2>1}
            nx2,ny2=g.shape
            defs2=[]
            for ii in range(nx2):
                for jj in range(ny2):
                    for di,dj in [(1,0),(0,1)]:
                        ni2,nj2=(ii+di)%nx2,(jj+dj)%ny2
                        i1=V2I_J1.get(g[ii,jj],-1)
                        i2=V2I_J1.get(g[ni2,nj2],-1)
                        if i1<0 or i2<0 or (i2-i1)%3 not in(1,2):
                            defs2.append((ii,jj,ni2,nj2))
            for ii,jj,ni2,nj2 in defs2:
                for ri,rj in [(ii,jj),(ni2,nj2)]:
                    if (ri,rj) in frozen or (ri,rj) in nf2: continue
                    ci=V2I_J1.get(g[ri,rj],1)
                    ni3=min(ci+1,2) if ci<1 else max(ci-1,0)
                    g2[ri,rj]=VALS_J1[ni3]; nf2[(ri,rj)]=2; break
            g,frozen=g2,nf2
            defs=defs2
            for d in defs: events.append((d[0],d[1],t))

        n_defs = len(defs) if 'defs' in dir() else 0
        if t % (t_steps // 4) == 0:
            print(f"  t={t:4d}  дефектов: {n_defs:6d}  событий: {len(events):8d}")

        if len(events) > 500000:
            print(f"  Обрезаем события до 500000")
            events = events[:500000]
            break

    print(f"  Итого событий: {len(events)}")

    if not events:
        print("  НЕТ СОБЫТИЙ")
        return None

    graph = build_causal_graph(events)
    d, max_r, n_reached = measure_dimension(graph, len(events), label)

    print(f"  Достигнуто узлов из корня: {n_reached}/{len(events)}")
    print(f"  Максимальный радиус: {max_r}")
    if d is not None:
        print(f"  >>> d = {d} <<<")
    else:
        print(f"  >>> d = НЕТ ДАННЫХ (мало связных узлов) <<<")

    return d


# ================================================================
# ЗАПУСК ВСЕХ СИМУЛЯЦИЙ
# ================================================================

if __name__ == "__main__":
    results = []

    # --- Sim 1: j=1/2, 3D, Rule B ---
    # Ожидание: d ≈ 4
    d = run_sim(
        label   = "Sim1: j=1/2, 3D, Rule B",
        nx=20, ny=20, nz=20,
        j_value = 0.5,
        n_defects = 80,
        t_steps   = 60,
    )
    results.append(("j=1/2, 3D", d))

    # --- Sim 2: j=1, 3D ---
    # Ожидание: d ≈ 4-5
    d = run_sim(
        label   = "Sim2: j=1, 3D",
        nx=20, ny=20, nz=20,
        j_value = 1.0,
        n_defects = 80,
        t_steps   = 60,
    )
    results.append(("j=1, 3D", d))

    # --- Sim 3: j=3/2, 2D ---
    # Ожидание: d ≈ 3-4
    d = run_sim(
        label   = "Sim3: j=3/2, 2D",
        nx=60, ny=60, nz=None,
        j_value = 1.5,
        n_defects = 60,
        t_steps   = 80,
    )
    results.append(("j=3/2, 2D", d))

    # --- Sim 4: j=3/2, 3D ---
    # Ожидание: d ≈ 4-5
    d = run_sim(
        label   = "Sim4: j=3/2, 3D",
        nx=20, ny=20, nz=20,
        j_value = 1.5,
        n_defects = 80,
        t_steps   = 60,
    )
    results.append(("j=3/2, 3D", d))

    # --- Sim 5: j=2, 2D ---
    # Ожидание: d ≈ 5
    d = run_sim(
        label   = "Sim5: j=2, 2D",
        nx=60, ny=60, nz=None,
        j_value = 2.0,
        n_defects = 60,
        t_steps   = 80,
    )
    results.append(("j=2, 2D", d))

    # ================================================================
    # ИТОГОВАЯ ТАБЛИЦА — СКОПИРУЙ И ОТПРАВЬ
    # ================================================================

    print("\n")
    print("=" * 60)
    print("ИТОГОВАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ")
    print("=" * 60)
    print(f"{'Симуляция':<20} {'d (измерено)':<15} {'2j+1 (теория)'}")
    print("-" * 60)

    theory = {
        "j=1/2, 3D": 2,
        "j=1, 3D":   3,
        "j=3/2, 2D": 4,
        "j=3/2, 3D": 4,
        "j=2, 2D":   5,
    }

    # предыдущие результаты для контекста
    prev = [
        ("j=1/2, 2D", 1.12,  2),
        ("j=1,   2D", 2.10,  3),
    ]
    for name, d_val, th in prev:
        print(f"  {name:<18} {d_val:<15.3f} {th}  (предыдущий прогон)")

    for name, d_val in results:
        th = theory.get(name, "?")
        d_str = f"{d_val:.3f}" if d_val is not None else "нет данных"
        print(f"  {name:<18} {d_str:<15} {th}")

    print("=" * 60)
    print("\nЗакон: d_spacetime = 2j + 1 (гипотеза)")
    print("Финитно-размерная поправка MM: ~+0.4-0.5")
    print("=" * 60)
