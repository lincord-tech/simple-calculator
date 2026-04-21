"""
Sigma Model v3
==============
Ключевое изменение: нет флага in_transition.
Узел может порождать нескольких потомков — по одному на каждый конфликт.
Это и есть ветвление → объём.

Правило: конфликт (i,j) с одинаковым sigma → оба порождают потомка
с противоположным sigma через T_MIN. Каждый потомок — новый узел.
"""

import heapq, math
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

T_MIN = 1

class Node:
    def __init__(self, nid, sigma, t, parent_id=None):
        self.id      = nid
        self.sigma   = sigma
        self.t       = t
        self.parent_id = parent_id
        self.children  = []

class Simulation:
    def __init__(self, max_nodes=2000):
        self.max_nodes   = max_nodes
        self.nodes       = {}
        self.graph       = nx.Graph()
        self.event_queue = []
        self.next_id     = 0
        self.scheduled   = set()   # пары уже запланированных конфликтов

    def add_node(self, sigma, t, parent_id=None):
        nid  = self.next_id; self.next_id += 1
        node = Node(nid, sigma, t, parent_id)
        self.nodes[nid] = node
        self.graph.add_node(nid, sigma=sigma, t=t)
        if parent_id is not None:
            self.graph.add_edge(parent_id, nid)
            self.nodes[parent_id].children.append(nid)
        return node

    def try_schedule(self, id_a, id_b, t):
        key = (min(id_a,id_b), max(id_a,id_b))
        if key not in self.scheduled:
            self.scheduled.add(key)
            heapq.heappush(self.event_queue, (t, id_a, id_b))

    def causal_neighbors(self, node):
        """Все узлы в окне ±T_MIN по времени."""
        return [nb for nid, nb in self.nodes.items()
                if nid != node.id and abs(nb.t - node.t) <= T_MIN]

    def spawn_child(self, parent, t_birth):
        """Родитель порождает потомка с -sigma."""
        if len(self.nodes) >= self.max_nodes:
            return None
        child = self.add_node(-parent.sigma, t_birth, parent_id=parent.id)
        # регистрируем связи и новые конфликты
        for nb in self.causal_neighbors(child):
            if not self.graph.has_edge(child.id, nb.id):
                self.graph.add_edge(child.id, nb.id)
            if nb.sigma == child.sigma:
                self.try_schedule(child.id, nb.id, t_birth + T_MIN)
        return child

    def run(self):
        # старт: антисимметричная пара
        a = self.add_node(-1, 0)
        b = self.add_node(+1, 0)
        self.graph.add_edge(a.id, b.id)

        # первый импульс: b порождает потомка
        c = self.spawn_child(b, T_MIN)               # c.sigma = -1
        # c конфликтует с a (оба -1)
        if c:
            self.try_schedule(c.id, a.id, T_MIN + T_MIN)

        while self.event_queue and len(self.nodes) < self.max_nodes:
            t_ev, id_a, id_b = heapq.heappop(self.event_queue)
            if id_a not in self.nodes or id_b not in self.nodes:
                continue
            na, nb = self.nodes[id_a], self.nodes[id_b]
            if na.sigma != nb.sigma:
                continue   # конфликт уже разрешён другим путём

            # ОБА порождают потомков — это и есть ветвление
            self.spawn_child(na, t_ev + T_MIN)
            self.spawn_child(nb, t_ev + T_MIN)

        n = len(self.nodes)
        e = self.graph.number_of_edges()
        print(f"Узлов: {n}  Рёбер: {e}  Средняя степень: {2*e/max(n,1):.2f}")

    # ---------- анализ ----------

    def measure_volume(self):
        if 0 not in self.graph: return [], []
        lengths = nx.single_source_shortest_path_length(self.graph, 0)
        max_r   = max(lengths.values())
        radii, volumes = [], []
        for r in range(1, max_r+1):
            radii.append(r)
            volumes.append(sum(1 for d in lengths.values() if d <= r))
        return radii, volumes

    def estimate_dim(self, radii, volumes):
        if len(radii) < 6: return None
        n    = len(radii)
        lo, hi = n//4, 3*n//4
        lr   = [math.log(radii[i])   for i in range(lo, hi)]
        lv   = [math.log(volumes[i]) for i in range(lo, hi)]
        return np.polyfit(lr, lv, 1)[0]

    # ---------- визуализация ----------

    def plot(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        fig.patch.set_facecolor('#0f0f1a')
        for ax in (ax1, ax2): ax.set_facecolor('#0f0f1a')

        # граф
        try:    pos = nx.kamada_kawai_layout(self.graph)
        except: pos = nx.spring_layout(self.graph, seed=42, k=0.5)

        colors = ['#ff4466' if self.nodes[n].sigma==1 else '#4488ff'
                  for n in self.graph.nodes()]
        nx.draw_networkx_edges(self.graph, pos, ax=ax1,
                               edge_color='#ffffff18', width=0.4)
        nx.draw_networkx_nodes(self.graph, pos, ax=ax1,
                               node_color=colors, node_size=15, alpha=0.9)
        ax1.set_title("Причинный граф  (красный=+1, синий=−1)",
                      color='white'); ax1.axis('off')

        # объём
        radii, volumes = self.measure_volume()
        d = self.estimate_dim(radii, volumes)
        if radii:
            ax2.plot(radii, volumes, 'o-', color='#44ff88', lw=2, ms=4)
            if d is not None:
                r_a  = np.array(radii, float)
                mid  = len(volumes)//2
                sc   = volumes[mid] / (radii[mid]**d + 1e-9)
                ax2.plot(r_a, sc*r_a**d, '--', color='#ffaa00',
                         alpha=.7, label=f'~ r^{d:.2f}')
                ax2.legend(facecolor='#1a1a2e', labelcolor='white')
                ax2.set_title(f"N(r)   d ≈ {d:.2f}", color='white')
            else:
                ax2.set_title("N(r)", color='white')
            ax2.set_xscale('log'); ax2.set_yscale('log')
            ax2.set_xlabel('r', color='white')
            ax2.set_ylabel('N(r)', color='white')
            ax2.tick_params(colors='white')

        plt.tight_layout()
        plt.savefig("/mnt/user-data/outputs/sigma_sim_result.png",
                    dpi=150, facecolor=fig.get_facecolor())
        print("График сохранён.")

if __name__ == "__main__":
    sim = Simulation(max_nodes=2000)
    sim.run()
    sim.plot()
