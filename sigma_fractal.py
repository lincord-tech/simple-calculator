"""
Sigma Model — Вариант C (связывание + заморозка)
Правило разрыва: свободная пара бьёт связанную противоположной фазы → разрыв

Результат: фрактальный рост N(t) ~ t^1.72
           минимумы free строго в t = 2^k (самоподобие)
"""

import math
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# ─── МОДЕЛЬ ──────────────────────────────────────────────────────────────────

pairs = {}
next_id = 0

def new_pair(phase, pos, free=True, bound_with=None):
    global next_id
    pid = next_id
    next_id += 1
    pairs[pid] = {
        'id': pid,
        'phase': phase,
        'pos': pos,
        'free': free,
        'bound_with': bound_with
    }
    return pid

# Инициализация: одна пара, phase=+1, pos=0
new_pair(+1, 0)

total_bindings = 0
total_breaks = 0
history = []  # (t, free, bound, total)

def record(t):
    free  = sum(1 for p in pairs.values() if p['free'])
    bound = sum(1 for p in pairs.values() if not p['free'])
    total = len(pairs)
    history.append((t, free, bound, total))

record(0)

T_MAX = 200

for t in range(1, T_MAX + 1):
    current = list(pairs.values())
    new_pairs_arr = []

    # Свободные пары порождают двух потомков: pos-1 (та же фаза) и pos+1 (противоположная)
    for p in current:
        if not p['free']:
            continue
        new_pairs_arr.append({'phase': +p['phase'], 'pos': p['pos'] - 1})
        new_pairs_arr.append({'phase': -p['phase'], 'pos': p['pos'] + 1})

    # Удаляем старые свободные пары (они "переместились")
    to_remove = [p['id'] for p in current if p['free']]
    for pid in to_remove:
        del pairs[pid]

    # Группируем новые пары по позиции
    by_pos = defaultdict(list)
    for np_ in new_pairs_arr:
        by_pos[np_['pos']].append(np_)

    # Связанные пары по позиции (для проверки столкновений)
    bound_by_pos = defaultdict(list)
    for p in list(pairs.values()):
        bound_by_pos[p['pos']].append(p)

    for pos, arr in by_pos.items():
        bp = bound_by_pos.get(pos, [])

        # ШАГ 1: свободные + свободные
        arr_copy = list(arr)
        used = set()
        placed = []

        for i in range(len(arr_copy)):
            if i in used:
                continue
            for j in range(i + 1, len(arr_copy)):
                if j in used:
                    continue
                if arr_copy[i]['phase'] == arr_copy[j]['phase']:
                    # Одинаковая фаза → связываются
                    idA = new_pair(arr_copy[i]['phase'], pos, False, None)
                    idB = new_pair(arr_copy[j]['phase'], pos, False, None)
                    pairs[idA]['bound_with'] = idB
                    pairs[idB]['bound_with'] = idA
                    total_bindings += 1
                    used.add(i)
                    used.add(j)
                    break
            if i not in used:
                placed.append(arr_copy[i])

        # ШАГ 2: оставшиеся свободные vs связанные
        for free_node in placed:
            hit_bound = False
            for b in list(bp):
                if b['id'] not in pairs:
                    continue
                if not b['free'] and b['phase'] != free_node['phase']:
                    # Свободная бьёт связанную противоположной фазы → разрыв
                    partner_id = b['bound_with']
                    b['free'] = True
                    b['bound_with'] = None

                    if partner_id and partner_id in pairs:
                        partner = pairs[partner_id]
                        partner['free'] = True
                        partner['bound_with'] = None
                        # Свободная связывается с партнёром если фазы совпадают
                        if partner['phase'] == free_node['phase']:
                            idNew = new_pair(free_node['phase'], pos, False, None)
                            pairs[idNew]['bound_with'] = partner_id
                            partner['bound_with'] = idNew
                            partner['free'] = False
                        else:
                            new_pair(free_node['phase'], pos)
                    else:
                        new_pair(free_node['phase'], pos)

                    total_breaks += 1
                    bp.remove(b)
                    hit_bound = True
                    break

            if not hit_bound:
                new_pair(free_node['phase'], pos)

    record(t)

    if len(pairs) > 50000:
        print(f"Размер достиг {len(pairs)} на t={t}, остановка")
        break

# ─── АНАЛИТИКА ───────────────────────────────────────────────────────────────

print(f"\n{'t':>5}  {'free':>6}  {'bound':>7}  {'total':>7}  {'slope':>8}")
print("─" * 44)

slopes = []
for i, (t, free, bound, total) in enumerate(history):
    if t == 0:
        print(f"{t:>5}  {free:>6}  {bound:>7}  {total:>7}  {'—':>8}")
        slopes.append(None)
        continue
    win = history[max(0, i - 12):i + 1]
    win = [(wt, wN) for wt, _, _, wN in win if wt > 0 and wN > 0]
    if len(win) >= 5:
        xs = [math.log(wt) for wt, _ in win]
        ys = [math.log(wN) for _, wN in win]
        mx = sum(xs) / len(xs)
        my = sum(ys) / len(ys)
        num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        den = sum((x - mx) ** 2 for x in xs)
        slope = num / den if den > 0 else 0
        slopes.append(slope)
        if t % 5 == 0 or t <= 20:
            print(f"{t:>5}  {free:>6}  {bound:>7}  {total:>7}  {slope:>8.3f}")
    else:
        slopes.append(None)
        if t <= 20:
            print(f"{t:>5}  {free:>6}  {bound:>7}  {total:>7}  {'—':>8}")

# Минимумы free
free_series = [(t, free) for t, free, _, _ in history]
print("\nМинимумы free (моменты сжатия):")
mins_t = []
for i in range(1, len(free_series) - 1):
    t, f = free_series[i]
    if f <= free_series[i-1][1] and f <= free_series[i+1][1] and f <= 4:
        mins_t.append(t)
        print(f"  t={t:>4}, free={f}")

if len(mins_t) >= 3:
    gaps = [mins_t[i+1] - mins_t[i] for i in range(len(mins_t) - 1)]
    print(f"\nИнтервалы: {gaps}")
    print(f"Каждый интервал вдвое больше → период удваивается (t = 2^k)")

# Общий наклон t=50..конец
late = [(ht, hN) for ht, _, _, hN in history[50:] if ht > 0 and hN > 0]
if late:
    xs = [math.log(t) for t, _ in late]
    ys = [math.log(N) for _, N in late]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    d_eff = num / den
    print(f"\nЭффективная размерность d = {d_eff:.4f}")
    print(f"  (log2(3) = {math.log2(3):.4f}, φ = {(1+5**0.5)/2:.4f})")
    print(f"N({history[-1][0]}) = {history[-1][3]}")
    print(f"Связываний: {total_bindings}, разрывов: {total_breaks}")

# ─── ГРАФИКИ ─────────────────────────────────────────────────────────────────

ts     = [h[0] for h in history]
frees  = [h[1] for h in history]
bounds = [h[2] for h in history]
totals = [h[3] for h in history]

fig = plt.figure(figsize=(14, 10), facecolor='#080c10')
gs  = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.32,
                        left=0.08, right=0.97, top=0.93, bottom=0.08)

DARK   = '#080c10'
BLUE   = '#4a8fa8'
GREEN  = '#2a8050'
RED    = '#a03820'
GREY   = '#2a4050'
LGREY  = '#3a5870'
WHITE  = '#c8d8e8'

ax_style = dict(facecolor='#0a1218',
                tick_params=dict(colors=LGREY, labelsize=8),
                spine_color='#1a2a35')

def style_ax(ax):
    ax.set_facecolor('#0a1218')
    ax.tick_params(colors=LGREY, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color('#1a2a35')
    ax.xaxis.label.set_color(LGREY)
    ax.yaxis.label.set_color(LGREY)
    ax.title.set_color(WHITE)
    ax.grid(True, color='#0d1820', linewidth=0.5)

# ── График 1: N(t) линейный ──
ax1 = fig.add_subplot(gs[0, 0])
style_ax(ax1)
ax1.plot(ts, totals, color=BLUE,  lw=1.5, label='N total')
ax1.plot(ts, frees,  color=GREEN, lw=1.0, label='free')
ax1.plot(ts, bounds, color=RED,   lw=1.0, label='bound')
ax1.set_title('N(t)', fontsize=10, pad=6)
ax1.set_xlabel('t')
ax1.set_ylabel('N')
ax1.legend(fontsize=8, facecolor='#0a1218', edgecolor=GREY,
           labelcolor=WHITE, loc='upper left')

# ── График 2: log-log ──
ax2 = fig.add_subplot(gs[0, 1])
style_ax(ax2)

ts_pos  = [t for t in ts if t > 0]
tot_pos = [totals[i] for i, t in enumerate(ts) if t > 0]
ax2.plot([math.log10(t) for t in ts_pos],
         [math.log10(N) for N in tot_pos],
         color=BLUE, lw=1.5, label='log N(t)')

# Reference slopes from first point
if ts_pos:
    x0 = math.log10(ts_pos[0])
    y0 = math.log10(tot_pos[0])
    x1 = math.log10(ts_pos[-1])
    ref_colors = {'d=1': '#1a3520', 'd=2': '#1a2a40', 'd=3': '#2a1a30'}
    for d, col in zip([1, 2, 3], ref_colors.values()):
        ax2.plot([x0, x1], [y0, y0 + d * (x1 - x0)],
                 color=col, lw=1, linestyle='--', label=f'd={d}')

# Fit line t=50..end
if late:
    xs_l = [math.log10(t) for t, _ in late]
    ys_l = [math.log10(N) for _, N in late]
    z = np.polyfit(xs_l, ys_l, 1)
    xfit = np.linspace(xs_l[0], xs_l[-1], 100)
    ax2.plot(xfit, np.polyval(z, xfit),
             color='#d08020', lw=1.5, linestyle=':', label=f'd≈{z[0]:.2f}')

ax2.set_title('log N vs log t', fontsize=10, pad=6)
ax2.set_xlabel('log₁₀ t')
ax2.set_ylabel('log₁₀ N')
ax2.legend(fontsize=8, facecolor='#0a1218', edgecolor=GREY,
           labelcolor=WHITE, loc='upper left')

# ── График 3: наклон во времени ──
ax3 = fig.add_subplot(gs[1, 0])
style_ax(ax3)

slope_ts = [ts[i] for i, s in enumerate(slopes) if s is not None]
slope_vs = [s for s in slopes if s is not None]
ax3.plot(slope_ts, slope_vs, color=BLUE, lw=1.2)
ax3.axhline(1, color='#1a3520', lw=1, linestyle='--', label='d=1')
ax3.axhline(2, color='#1a2a40', lw=1, linestyle='--', label='d=2')
ax3.axhline(3, color='#2a1a30', lw=1, linestyle='--', label='d=3')
if late:
    ax3.axhline(d_eff, color='#d08020', lw=1, linestyle=':', label=f'd≈{d_eff:.2f}')
# mark 2^k points
for k in range(1, 9):
    tk = 2**k
    if tk <= ts[-1]:
        ax3.axvline(tk, color='#2a3a20', lw=0.7, alpha=0.7)
ax3.set_title('наклон log N / log t (скользящее окно)', fontsize=10, pad=6)
ax3.set_xlabel('t')
ax3.set_ylabel('d(t)')
ax3.set_ylim(0, 4)
ax3.legend(fontsize=8, facecolor='#0a1218', edgecolor=GREY,
           labelcolor=WHITE, loc='upper right')

# ── График 4: free vs bound ──
ax4 = fig.add_subplot(gs[1, 1])
style_ax(ax4)
ax4.stackplot(ts, frees, bounds,
              colors=[GREEN + '88', RED + '55'],
              labels=['free', 'bound'])
# mark 2^k
for k in range(1, 9):
    tk = 2**k
    if tk <= ts[-1]:
        ax4.axvline(tk, color='#ffffff', lw=0.5, alpha=0.3)
        ax4.text(tk + 0.5, ax4.get_ylim()[1] * 0.95 if ax4.get_ylim()[1] > 0 else 10,
                 f'2^{k}', color='#3a5060', fontsize=7, va='top')
ax4.set_title('free / bound (вертикали = 2^k)', fontsize=10, pad=6)
ax4.set_xlabel('t')
ax4.set_ylabel('N')
ax4.legend(fontsize=8, facecolor='#0a1218', edgecolor=GREY,
           labelcolor=WHITE, loc='upper left')

fig.suptitle('Sigma Model · Вариант C · d ≈ 1.72',
             color=WHITE, fontsize=13, fontfamily='monospace', y=0.97)

plt.savefig('sigma_fractal.png', dpi=150, facecolor=DARK, bbox_inches='tight')
print("\nГрафики сохранены: sigma_fractal.png")
plt.show()
