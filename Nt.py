import numpy as np

# -----------------------
# Модель (2D, ближайшие соседи)
# -----------------------
class Sigma2D:
    def __init__(self, Lx, Ly, tau=2, seed=None):
        self.shape = (Lx, Ly)
        self.tau = tau
        self.rng = np.random.default_rng(seed)

        # антиферромагнитный вакуум
        grid = np.indices(self.shape)
        parity = (grid[0] + grid[1]) % 2
        self.sigma = np.where(parity == 0, 1, -1)

        self.in_trans = np.zeros(self.shape, dtype=int)
        self.target = np.zeros(self.shape, dtype=int)

        self.events = []  # (t, (i,j))

    def neighbors(self, i, j):
        for di, dj in [(1,0),(-1,0),(0,1),(0,-1)]:
            ni, nj = i+di, j+dj
            if 0 <= ni < self.shape[0] and 0 <= nj < self.shape[1]:
                yield ni, nj

    def frozen(self, i, j):
        if self.sigma[i,j] == 0:
            return True
        for ni, nj in self.neighbors(i,j):
            if self.in_trans[ni,nj] > 0:
                return True
        return False

    def find_defects(self):
        defects = []
        Lx, Ly = self.shape
        for i in range(Lx):
            for j in range(Ly):
                for ni, nj in self.neighbors(i,j):
                    if self.sigma[i,j] != 0 and self.sigma[ni,nj] != 0:
                        if self.sigma[i,j] == self.sigma[ni,nj]:
                            defects.append(((i,j),(ni,nj)))
        return defects

    def seed_flips(self, k):
        # k случайных флипов в вакууме
        Lx, Ly = self.shape
        idxs = set()
        while len(idxs) < k:
            idxs.add((self.rng.integers(0,Lx), self.rng.integers(0,Ly)))
        for (i,j) in idxs:
            self.sigma[i,j] *= -1

    def step(self, t):
        # завершение
        Lx, Ly = self.shape
        for i in range(Lx):
            for j in range(Ly):
                if self.in_trans[i,j] == t:
                    self.sigma[i,j] = self.target[i,j]
                    self.in_trans[i,j] = 0

        # запуск
        for (a, b) in self.find_defects():
            for (i,j) in [a, b]:
                if self.in_trans[i,j] == 0 and not self.frozen(i,j):
                    self.target[i,j] = -self.sigma[i,j]
                    self.sigma[i,j] = 0
                    self.in_trans[i,j] = t + self.tau
                    self.events.append((t, (i,j)))
                    break

    def run(self, T):
        for t in range(T):
            self.step(t)

# -----------------------
# Анализ: N(t) и d(t)
# -----------------------
def Nt_series(events, T):
    # кумулятивное число событий
    Nt = np.zeros(T, dtype=float)
    for (t, _) in events:
        if t < T:
            Nt[t] += 1
    Nt = np.cumsum(Nt)
    return Nt

def local_slope_d(times, Nt, window=30, eps=1e-12):
    # d(t) как наклон ln N vs ln t на скользящем окне
    ts = np.array(times, dtype=float)
    Ns = np.array(Nt, dtype=float)
    dvals = np.full_like(ts, fill_value=np.nan, dtype=float)

    # избегаем нулей
    mask = (ts > 0) & (Ns > 0)
    ts = ts[mask]
    Ns = Ns[mask]

    for i in range(len(ts)):
        j0 = max(0, i - window//2)
        j1 = min(len(ts), i + window//2)
        if j1 - j0 < 5:
            continue
        x = np.log(ts[j0:j1] + eps)
        y = np.log(Ns[j0:j1] + eps)
        a, b = np.polyfit(x, y, 1)
        dvals[mask][i] = a

    # вернуть в исходной длине (с nan где нет данных)
    out = np.full(len(times), np.nan)
    out[mask] = dvals
    return out

def plateau_estimate(d_t, t_min_frac=0.5):
    # среднее d(t) на «поздних временах»
    n = len(d_t)
    start = int(n * t_min_frac)
    tail = d_t[start:]
    tail = tail[~np.isnan(tail)]
    if len(tail) == 0:
        return np.nan, np.nan
    return float(np.mean(tail)), float(np.std(tail))

# -----------------------
# Эксперимент (ансамбль)
# -----------------------
def run_ensemble(L=80, T=400, k_list=(3,5,10), runs=5, tau=2, window=40):
    results = {}
    for k in k_list:
        ds = []
        for r in range(runs):
            model = Sigma2D(L, L, tau=tau, seed=1234 + r + 100*k)
            model.seed_flips(k)
            model.run(T)

            Nt = Nt_series(model.events, T)
            times = np.arange(T)

            d_t = local_slope_d(times, Nt, window=window)
            d_mean, d_std = plateau_estimate(d_t, t_min_frac=0.6)

            ds.append((d_mean, d_std, Nt[-1]))

        results[k] = ds
    return results

# -----------------------
# Запуск
# -----------------------
if __name__ == "__main__":
    res = run_ensemble(L=80, T=400, k_list=(3,5,10,20), runs=6, tau=2, window=50)

    for k, vals in res.items():
        means = [v[0] for v in vals if not np.isnan(v[0])]
        stds  = [v[1] for v in vals if not np.isnan(v[1])]
        Ns    = [v[2] for v in vals]
        if len(means) == 0:
            print(f"k={k}: insufficient data")
            continue
        print(f"k={k}: d ≈ {np.mean(means):.3f} ± {np.mean(stds):.3f}  |  events ~ {int(np.mean(Ns))}")