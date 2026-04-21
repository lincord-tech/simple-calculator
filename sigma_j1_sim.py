"""
Sigma Model — j=1 динамика
===========================
При j=1: три равноправных состояния {-1, 0, +1}.
0 теперь не переходное — это полноценное состояние.
Переходы: -1 ↔ 0 ↔ +1, прямой -1 ↔ +1 запрещён.

Вакуум: σ(i,j) = (i+j) mod 3 - 1  → {-1, 0, +1} чередуются.
Дефект: нарушение этого порядка.

Измеряем d методом N(r) ~ r^d из причинного графа.
Сравниваем с j=1/2 результатом.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from collections import defaultdict, deque

# ───────────────────────────── параметры ─────────────────────────────
NX, NY = 60, 60
T_STEPS = 80
N_DEFECTS = 20
J_VALUE = 1   # меняй на 0.5 для сравнения с j=1/2

# ───────────────────────── вакуум j=1 ────────────────────────────────
def make_vacuum_j1(nx, ny):
    """Трёхзначный шахматный вакуум: (i+j) mod 3 → {-1,0,+1}"""
    grid = np.zeros((nx, ny), dtype=int)
    for i in range(nx):
        for j in range(ny):
            grid[i, j] = (i + j) % 3 - 1
    return grid

def make_vacuum_j05(nx, ny):
    """Двузначный шахматный вакуум: (-1)^(i+j)"""
    grid = np.zeros((nx, ny), dtype=int)
    for i in range(nx):
        for j in range(ny):
            grid[i, j] = 1 if (i + j) % 2 == 0 else -1
    return grid

# ───────────────────────── дефекты ───────────────────────────────────
def find_defects_j1(grid):
    """Дефект = ребро где разность не равна ±1 по mod 3."""
    nx, ny = grid.shape
    defects = []
    for i in range(nx):
        for j in range(ny):
            for di, dj in [(1,0),(0,1)]:
                ni, nj = (i+di)%nx, (j+dj)%ny
                diff = (grid[ni,nj] - grid[i,j]) % 3
                if diff != 1 and diff != 2:  # не ±1 mod 3
                    defects.append((i,j,ni,nj))
    return defects

def find_defects_j05(grid):
    """Дефект = ребро где σ_i * σ_j = +1."""
    nx, ny = grid.shape
    defects = []
    for i in range(nx):
        for j in range(ny):
            for di, dj in [(1,0),(0,1)]:
                ni, nj = (i+di)%nx, (j+dj)%ny
                if grid[i,j] * grid[ni,nj] == 1:
                    defects.append((i,j,ni,nj))
    return defects

# ───────────────────────── динамика ──────────────────────────────────
def step_j1(grid, frozen):
    """j=1: переход через соседнее значение. -1→0→+1 или +1→0→-1."""
    nx, ny = grid.shape
    new_grid = grid.copy()
    new_frozen = {}

    # обновляем заморозку
    for pos, t in frozen.items():
        if t > 1:
            new_frozen[pos] = t - 1

    defects = find_defects_j1(grid)
    for i, j, ni, nj in defects:
        if (i,j) in frozen or (ni,nj) in frozen:
            continue
        # узел с большим значением уменьшается на 1 (шаг к соседу)
        if grid[i,j] > grid[ni,nj]:
            reactor = (i,j)
        else:
            reactor = (ni,nj)
        ri, rj = reactor
        if reactor not in new_frozen:
            # шаг к нулю или от нуля
            cur = grid[ri,rj]
            if cur == 1:
                new_grid[ri,rj] = 0
            elif cur == -1:
                new_grid[ri,rj] = 0
            elif cur == 0:
                # выбираем направление из соседей
                neighbors = [grid[(ri+di)%nx, (rj+dj)%ny]
                             for di,dj in [(1,0),(-1,0),(0,1),(0,-1)]]
                avg = sum(neighbors)
                new_grid[ri,rj] = 1 if avg > 0 else -1
            new_frozen[reactor] = 2  # заморозка на 2 шага

    return new_grid, new_frozen

def step_j05(grid, frozen):
    """j=1/2: стандартная динамика σ→0→-σ."""
    nx, ny = grid.shape
    new_grid = grid.copy()
    new_frozen = {}

    for pos, t in frozen.items():
        if t > 1:
            new_frozen[pos] = t - 1

    defects = find_defects_j05(grid)
    for i, j, ni, nj in defects:
        if (i,j) in frozen or (ni,nj) in frozen:
            continue
        reactor = (i,j)
        ri, rj = reactor
        if reactor not in new_frozen:
            new_grid[ri,rj] = -grid[ri,rj]
            new_frozen[reactor] = 2

    return new_grid, new_frozen

# ───────────────────── причинный граф и d ────────────────────────────
def build_causal_graph(events):
    """events: list of (i, j, t). Связь если |t2-t1|=1 и (i,j) рядом."""
    graph = defaultdict(set)
    by_t = defaultdict(list)
    for idx, (i, j, t) in enumerate(events):
        by_t[t].append((idx, i, j))

    all_t = sorted(by_t.keys())
    for k, t in enumerate(all_t[:-1]):
        t_next = all_t[k+1]
        if t_next - t != 1:
            continue
        for idx1, i1, j1 in by_t[t]:
            for idx2, i2, j2 in by_t[t_next]:
                if abs(i1-i2) + abs(j1-j2) <= 1:
                    graph[idx1].add(idx2)
                    graph[idx2].add(idx1)
    return graph

def measure_dimension(graph, n_events):
    """N(r) ~ r^d из BFS от случайного корня."""
    if n_events < 10:
        return None, [], []

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
    if max_r < 3:
        return None, [], []

    radii = list(range(1, max_r+1))
    volumes = [sum(1 for d in dist.values() if d <= r) for r in radii]

    # оценка наклона в лог-лог (средняя часть)
    n = len(radii)
    lo, hi = n//4, 3*n//4
    if hi <= lo + 1:
        return None, radii, volumes
    lr = [math.log(radii[i]) for i in range(lo, hi)]
    lv = [math.log(volumes[i]) for i in range(lo, hi) if volumes[i] > 0]
    if len(lr) != len(lv):
        lr = lr[:len(lv)]
    if len(lr) < 2:
        return None, radii, volumes
    coeffs = np.polyfit(lr, lv, 1)
    return coeffs[0], radii, volumes

# ───────────────────────────── запуск ────────────────────────────────
def run_sim(j_value):
    label = f"j={'1/2' if j_value==0.5 else '1'}"
    print(f"\n=== {label} ===")

    if j_value == 1:
        grid = make_vacuum_j1(NX, NY)
        find_def = find_defects_j1
        step_fn  = step_j1
    else:
        grid = make_vacuum_j05(NX, NY)
        find_def = find_defects_j05
        step_fn  = step_j05

    # вносим дефекты случайно
    rng = np.random.default_rng(42)
    for _ in range(N_DEFECTS):
        i, j = rng.integers(1, NX-1), rng.integers(1, NY-1)
        if j_value == 1:
            grid[i,j] = rng.choice([-1,0,1])
        else:
            grid[i,j] = rng.choice([-1,1])

    frozen = {}
    events = []  # (i, j, t)

    for t in range(T_STEPS):
        grid, frozen = step_fn(grid, frozen)
        defs = find_def(grid)
        for di,dj,ni,nj in defs:
            events.append((di, dj, t))

        if t % 20 == 0:
            print(f"  t={t:3d}  дефектов: {len(defs):4d}  "
                  f"событий всего: {len(events):6d}")

    print(f"  Итого событий: {len(events)}")
    graph = build_causal_graph(events)
    d, radii, volumes = measure_dimension(graph, len(events))

    if d is not None:
        print(f"  Оценка размерности d ≈ {d:.3f}")
    else:
        print("  Недостаточно данных для оценки d")

    return d, radii, volumes, label

# ───────────────────────────── plot ──────────────────────────────────
def plot(results):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor('#0f0f1a')
    for ax in axes:
        ax.set_facecolor('#0f0f1a')

    colors = {'j=1/2': '#44ff88', 'j=1': '#ff4466'}

    for d, radii, volumes, label in results:
        if not radii:
            continue
        col = colors.get(label, '#ffffff')
        axes[0].plot(radii, volumes, 'o-', color=col,
                     lw=1.5, ms=3, label=label, alpha=0.85)
        axes[1].plot(radii, volumes, 'o', color=col,
                     ms=2, alpha=0.6, label=label)

        if d is not None:
            r_a = np.array(radii, float)
            mid = len(volumes)//2
            sc = volumes[mid] / (radii[mid]**d + 1e-9)
            axes[1].plot(r_a, sc * r_a**d, '--', color=col,
                         alpha=0.7, label=f'{label}: d≈{d:.2f}')

    for ax in axes:
        ax.tick_params(colors='white')
        ax.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=9)

    axes[0].set_xlabel('r', color='white')
    axes[0].set_ylabel('N(r)', color='white')
    axes[0].set_title('Рост объёма N(r)', color='white')

    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].set_xlabel('r (лог)', color='white')
    axes[1].set_ylabel('N(r) (лог)', color='white')
    axes[1].set_title('Лог-лог: оценка размерности', color='white')

    plt.tight_layout()
    out = "/mnt/user-data/outputs/sigma_j_comparison.png"
    plt.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print(f"\nГрафик: {out}")

# ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    results = []
    for j in [0.5, 1.0]:
        d, radii, volumes, label = run_sim(j)
        results.append((d, radii, volumes, label))

    print("\n=== Итог ===")
    for d, _, _, label in results:
        if d is not None:
            print(f"  {label}: d ≈ {d:.3f}")
        else:
            print(f"  {label}: d не определена")

    plot(results)
