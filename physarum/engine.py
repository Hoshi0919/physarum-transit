"""
Physarum Transit · 黏菌自适应网络动力学核心 (Physarum Polycephalum Dynamics)
基于 Tero et al. (Science 2010) 与 Nakagaki et al. (Nature 2000) 经典理论：
- 泊肃叶流体通量 (Poiseuille Tube Flux)
- 节点质量守恒压力泊松方程 (Pressure Poisson Relaxation)
- 非线性正反馈管径自适应与萎缩修剪 (Nonlinear Tube Adaptation)
"""

import math
from typing import Dict, List, Tuple, Optional, Set

class Node:
    def __init__(self, node_id: str, x: float, y: float, injection: float = 0.0, label: str = "", kind: str = "relay"):
        self.id = node_id
        self.x = x
        self.y = y
        self.injection = float(injection)  # 正值为源 (Source)，负值为汇 (Sink)
        self.label = label or node_id
        self.kind = kind  # "source", "sink", "relay", "hub", "food"
        self.pressure = 0.0

    def to_dict(self):
        return {
            "id": self.id,
            "x": self.x,
            "y": self.y,
            "injection": self.injection,
            "label": self.label,
            "kind": self.kind,
            "pressure": self.pressure
        }


class Edge:
    def __init__(self, u: str, v: str, length: float, conductivity: float = 1.0, d_min: float = 0.005, d_max: float = 10.0):
        self.u = min(u, v)
        self.v = max(u, v)
        self.key = f"{self.u}--{self.v}"
        self.length = max(0.001, float(length))
        self.conductivity = max(d_min, min(d_max, float(conductivity)))
        self.d_min = d_min
        self.d_max = d_max
        self.flux = 0.0  # 正表示 u -> v，负表示 v -> u
        self.active = True

    @property
    def resistance(self) -> float:
        if not self.active or self.conductivity <= self.d_min * 1.05:
            return float('inf')
        return self.length / self.conductivity

    def to_dict(self):
        return {
            "u": self.u,
            "v": self.v,
            "key": self.key,
            "length": self.length,
            "conductivity": self.conductivity,
            "flux": self.flux,
            "active": self.active
        }


class PhysarumNetwork:
    def __init__(self):
        self.nodes: Dict[str, Node] = {}
        self.edges: Dict[str, Edge] = {}
        self.adj: Dict[str, Set[str]] = {}
        self.step_count: int = 0

    def add_node(self, node_id: str, x: float, y: float, injection: float = 0.0, label: str = "", kind: str = "relay") -> Node:
        node = Node(node_id, x, y, injection, label, kind)
        self.nodes[node_id] = node
        if node_id not in self.adj:
            self.adj[node_id] = set()
        return node

    def add_edge(self, u: str, v: str, conductivity: float = 1.0, length: Optional[float] = None) -> Edge:
        if u not in self.nodes or v not in self.nodes:
            raise ValueError(f"Nodes {u} and {v} must exist in the network")
        if u == v:
            raise ValueError("Self-loops are not allowed")

        n1, n2 = self.nodes[u], self.nodes[v]
        if length is None:
            length = math.hypot(n1.x - n2.x, n1.y - n2.y)
        
        edge = Edge(u, v, length, conductivity)
        self.edges[edge.key] = edge
        self.adj[u].add(v)
        self.adj[v].add(u)
        return edge

    def sever_edge(self, u: str, v: str) -> bool:
        """用手术刀切断管道"""
        key = f"{min(u, v)}--{max(u, v)}"
        if key in self.edges:
            self.edges[key].active = False
            self.edges[key].conductivity = self.edges[key].d_min
            return True
        return False

    def solve_pressures(self, max_iters: int = 400, tol: float = 1e-5, omega: float = 1.3) -> int:
        """
        超松弛迭代法 (Successive Over-Relaxation, SOR) 求解泊松网络方程：
        对于节点 i: sum_j [ (D_ij / L_ij) * (p_i - p_j) ] = I_i
        => p_i = ( I_i + sum_j [ (D_ij / L_ij) * p_j ] ) / sum_j (D_ij / L_ij)
        """
        if not self.nodes:
            return 0

        node_keys = list(self.nodes.keys())
        anchor_id = node_keys[0]
        # 若存在 sink 节点，优先选 sink 为零电位参考点
        for nid, n in self.nodes.items():
            if n.injection < 0:
                anchor_id = nid
                break

        # 归一化注入电流以保证电荷守恒 sum(I) == 0
        total_inj = sum(n.injection for n in self.nodes.values())
        non_anchor_count = max(1, len(self.nodes) - 1)
        inj_correction = total_inj / len(self.nodes)

        iters = 0
        for it in range(max_iters):
            max_diff = 0.0
            for nid in node_keys:
                if nid == anchor_id:
                    self.nodes[nid].pressure = 0.0
                    continue

                node = self.nodes[nid]
                sum_cond = 0.0
                sum_flow = 0.0

                for neighbor_id in self.adj[nid]:
                    ekey = f"{min(nid, neighbor_id)}--{max(nid, neighbor_id)}"
                    edge = self.edges.get(ekey)
                    if edge and edge.active:
                        conductance = edge.conductivity / edge.length
                        sum_cond += conductance
                        sum_flow += conductance * self.nodes[neighbor_id].pressure

                if sum_cond > 1e-7:
                    # 修正注入
                    eff_inj = node.injection - inj_correction
                    target_p = (eff_inj + sum_flow) / sum_cond
                    # SOR 超松弛
                    new_p = (1.0 - omega) * node.pressure + omega * target_p
                    diff = abs(new_p - node.pressure)
                    if diff > max_diff:
                        max_diff = diff
                    node.pressure = new_p

            iters += 1
            if max_diff < tol:
                break

        self.compute_fluxes()
        return iters

    def compute_fluxes(self):
        """计算每根管道内的流体通量 Q_uv"""
        for edge in self.edges.values():
            if not edge.active:
                edge.flux = 0.0
                continue
            pu = self.nodes[edge.u].pressure
            pv = self.nodes[edge.v].pressure
            # 泊肃叶通量：Q = D / L * (p_u - p_v)
            edge.flux = (edge.conductivity / edge.length) * (pu - pv)

    def step_adaptation(self, dt: float = 0.1, alpha: float = 0.15, gamma: float = 1.3, decay: float = 0.25) -> float:
        """
        黏菌管径自适应演化方程 (Tero et al. 2010):
        dD/dt = alpha * |Q|^gamma - decay * D
        高通量管道超线性正反馈增粗，低通量劣势管道自然萎缩凋亡。
        """
        self.solve_pressures()
        max_delta = 0.0

        for edge in self.edges.values():
            if not edge.active:
                continue

            abs_q = abs(edge.flux)
            growth = alpha * math.pow(abs_q, gamma)
            shrink = decay * edge.conductivity
            dD = (growth - shrink) * dt

            new_d = edge.conductivity + dD
            new_d = max(edge.d_min, min(edge.d_max, new_d))
            delta = abs(new_d - edge.conductivity)
            if delta > max_delta:
                max_delta = delta
            edge.conductivity = new_d

        self.step_count += 1
        return max_delta

    def metrics(self) -> Dict[str, float]:
        """计算全网络健康度与拓扑指标"""
        total_d = sum(e.conductivity for e in self.edges.values())
        total_flux = sum(abs(e.flux) for e in self.edges.values())
        active_edges = sum(1 for e in self.edges.values() if e.active and e.conductivity > e.d_min * 1.5)
        
        # 能量损耗 (Power Dissipation): E = sum (L / D * Q^2)
        dissipation = 0.0
        for e in self.edges.values():
            if e.active and e.conductivity > e.d_min:
                dissipation += (e.length / e.conductivity) * (e.flux ** 2)

        # 连通性冗余度 (Biconnectivity / Robustness proxy)
        robustness = (active_edges / max(1, len(self.nodes))) if self.nodes else 0.0

        return {
            "active_edges": active_edges,
            "total_edges": len(self.edges),
            "total_conductivity": total_d,
            "total_flux": total_flux,
            "energy_dissipation": dissipation,
            "network_robustness": robustness,
            "steps": self.step_count
        }

    @classmethod
    def create_braess_network(cls) -> "PhysarumNetwork":
        """
        布雷斯悖论 (Braess's Paradox) 场景：
        4 个节点：S (源), A (上中继), B (下中继), T (汇)
        初始两条经典路径：S->A->T 与 S->B->T
        中途连通零成本高速快捷路 A->B。
        """
        net = cls()
        net.add_node("S", 100, 300, injection=10.0, label="起点 S", kind="source")
        net.add_node("A", 300, 150, injection=0.0, label="中继 A", kind="relay")
        net.add_node("B", 300, 450, injection=0.0, label="中继 B", kind="relay")
        net.add_node("T", 500, 300, injection=-10.0, label="终点 T", kind="sink")

        # 初始四条主干
        net.add_edge("S", "A", conductivity=1.0)
        net.add_edge("S", "B", conductivity=1.0)
        net.add_edge("A", "T", conductivity=1.0)
        net.add_edge("B", "T", conductivity=1.0)
        # 快捷交叉路 (可手动或默认低导通)
        net.add_edge("A", "B", conductivity=0.01)
        return net

    @classmethod
    def create_tokyo_network(cls) -> "PhysarumNetwork":
        """
        东京铁路网经典复现 (Tero et al. Science 2010):
        东京站为中心源 (Tokyo Central)，周边 10 个主要卫星核心节点 (横滨、千叶、埼玉、八王子等)，
        并在节点间建立 Delaunay / 紧密候选网格，经黏菌演化形成山手线环路与放射状干线！
        """
        net = cls()
        # 中心枢纽 (Tokyo Central)
        net.add_node("Tokyo", 400, 300, injection=30.0, label="东京 (Tokyo)", kind="hub")

        # 卫星主要城市与微弱汇流
        satellites = [
            ("Yokohama", 360, 430, -5.0, "横滨"),
            ("Kawasaki", 380, 370, -3.0, "川崎"),
            ("Chiba", 530, 310, -4.0, "千叶"),
            ("Funabashi", 470, 270, -3.0, "船桥"),
            ("Matsudo", 450, 220, -3.0, "松户"),
            ("Saitama", 370, 190, -5.0, "埼玉"),
            ("Kawagoe", 310, 180, -2.0, "川越"),
            ("Tokorozawa", 320, 240, -2.0, "所泽"),
            ("Hachioji", 260, 320, -3.0, "八王子"),
            ("Machida", 320, 360, -3.0, "町田"),
        ]
        for nid, x, y, inj, lbl in satellites:
            net.add_node(nid, x, y, injection=inj, label=lbl, kind="sink")

        # 连接放射线与环状候选线
        for nid, _, _, _, _ in satellites:
            net.add_edge("Tokyo", nid, conductivity=0.5)

        # 环状相邻连线 (山手线与外环线候选)
        ring_order = ["Yokohama", "Kawasaki", "Machida", "Hachioji", "Tokorozawa", "Kawagoe", "Saitama", "Matsudo", "Funabashi", "Chiba"]
        for i in range(len(ring_order)):
            u = ring_order[i]
            v = ring_order[(i + 1) % len(ring_order)]
            net.add_edge(u, v, conductivity=0.3)

        # 交叉连线
        net.add_edge("Kawasaki", "Machida", conductivity=0.2)
        net.add_edge("Saitama", "Funabashi", conductivity=0.2)
        net.add_edge("Tokorozawa", "Saitama", conductivity=0.2)

        return net

    @classmethod
    def create_factory_ring(cls) -> "PhysarumNetwork":
        """
        工场自愈物流环 (Fault-Tolerant Factory Ring):
        模拟 Minecraft 自动化或工业生产中的 6 个功能集群：
        - 矿物输入 (Ore Input, Source)
        - 高炉冶炼 (Smelting)
        - 晶片晶圆制造 (Wafer Fab)
        - 装配装配中心 (Assembly)
        - 发电与储能 (Power Grid)
        - 发射井交付仓 (Delivery Silo, Sink)
        具有内外双环与交叉互备管道，验证主干被切断时的自愈分流。
        """
        net = cls()
        net.add_node("Ore", 150, 300, injection=15.0, label="矿物总进料", kind="source")
        net.add_node("Smelt", 280, 180, injection=0.0, label="高炉冶炼", kind="relay")
        net.add_node("Fab", 440, 180, injection=0.0, label="晶片制造", kind="relay")
        net.add_node("Power", 280, 420, injection=0.0, label="裂变电网", kind="relay")
        net.add_node("Assem", 440, 420, injection=0.0, label="精密装配", kind="relay")
        net.add_node("Silo", 570, 300, injection=-15.0, label="发射井交付", kind="sink")

        # 拓扑连接
        net.add_edge("Ore", "Smelt", conductivity=1.0)
        net.add_edge("Ore", "Power", conductivity=1.0)
        net.add_edge("Smelt", "Fab", conductivity=1.0)
        net.add_edge("Power", "Assem", conductivity=1.0)
        net.add_edge("Fab", "Silo", conductivity=1.0)
        net.add_edge("Assem", "Silo", conductivity=1.0)
        # 交叉备份管道
        net.add_edge("Smelt", "Assem", conductivity=0.2)
        net.add_edge("Power", "Fab", conductivity=0.2)
        net.add_edge("Ore", "Silo", conductivity=0.05)  # 远距低容备用

        return net
