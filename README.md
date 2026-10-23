# Physarum Transit · 黏菌自适应网络与拟态物流模拟器

> 基于多头绒泡黏菌（*Physarum polycephalum*）自组织流体动力学构建的自适应路网演化与无闭锁自愈物流模拟器（Python 数学仿真 + HTML5 Canvas & Web Audio 实时交互双栈）。

[![Tests](https://img.shields.io/badge/tests-15%20passed-brightgreen.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)]()
[![Live Demo](https://img.shields.io/badge/demo-GitHub%20Pages-success.svg)](https://hoshi0919.github.io/physarum-transit/)

---

## 🌌 为什么是黏菌？（从生物智能到城市交通与工业物流）

自然界中的多头绒泡黏菌没有大脑、没有中枢神经系统，更没有全局路网地图。然而在经典生物实验（Tero et al., *Science* 2010 与 Nakagaki et al., *Nature* 2000）中，黏菌不仅在培养皿中几乎完美复现了**东京乃至关东地区的 JR 轨道交通网络拓扑**，还能高效求解复杂迷宫的最短测地线（Geodesic），并在多目标帕累托前沿（成本、传输效率与容错冗余）上超越了传统人工规划。

在游戏与系统工程中，这一机制直击两大核心领域：
1. **《城市天际线》(Cities: Skylines) 与城市路网规划**：
   - 传统自私路径选择（Selfish Routing）极易诱发布雷斯悖论（Braess's Paradox，新增道路反使全城拥堵加剧）；
   - 黏菌流体拓扑具有天然的全局压力均衡与自适应管径调节机制，能动态抑制恶性支线。
2. **Minecraft 科技整合包与自动化工业物流（Tech Modpacks）**：
   - 复杂装配线与多原料配比在面临单机停摆或局部断管时，极易发生背压死锁（Backpressure Deadlock）；
   - 黏菌拟态管道具备**瞬时自愈重路由能力**：当主干管路被切断时，管网压力场毫秒级重构，次级休眠管径迅速膨胀分流，实现零死锁运行。

---

## 📐 核心动力学与数学方程

模拟器忠实复现了生物力学经典模型：

1. **泊肃叶流体通量方程 (Poiseuille Tube Flux)**：
   每根黏菌管道 $(i, j)$ 具有欧氏几何长度 $L_{ij}$ 与动态导电/导通度 $D_{ij}$。管道两端节点压力差驱动流体通量：
   $$Q_{ij} = \frac{D_{ij}}{L_{ij}} (p_i - p_j)$$

2. **节点质量连续性方程与泊松松弛 (Mass Conservation & Poisson Relaxation)**：
   在任意非源汇中继节点处，流入与流出净通量严格守恒：
   $$\sum_{j \in \mathcal{N}(i)} Q_{ij} = I_i$$
   其中 $I_i > 0$ 为生产营养源（Source），$I_i < 0$ 为消耗汇（Sink），$I_i = 0$ 为中继枢纽。模拟器采用超松弛迭代法（Successive Over-Relaxation, SOR，松弛因子 $\omega = 1.35$）实时求解非对称压力势场。

3. **非线性管径自适应演化方程 (Adaptive Feedback Dynamics)**：
   管径随流体剪切应力与通量正反馈扩张，随代谢消耗持续自然萎缩：
   $$\frac{d D_{ij}}{dt} = \alpha |Q_{ij}|^\gamma - \mu D_{ij}$$
   - **$\gamma = 1.0$（线性适应）**：保留较多冗余环路，兼顾高度容错鲁棒性；
   - **$\gamma > 1.2$（超线性增益）**：高通量路径以超线性速度掠夺资源，促使次要迂回死胡同以指数速度萎缩，自动收敛为极简高速干线与最短测地线！

---

## 🎮 经典动力学情景 (Scenarios)

1. **🇯🇵 东京城铁复现 (Tokyo Railway Grid)**：
   - 以东京站为核心源点，环布横滨、川崎、千叶、埼玉、八王子等 10 大卫星城镇；
   - 预埋致密候选网格，在 10 秒演化内自发涌现**山手线环状骨架**与主要放射性干线！
2. **🧭 迷宫拓扑寻优 (Labyrinth Geodesic)**：
   - 迷宫入口与出口分别注入正负流体通量；
   - 初始全域探索所有通道，随着演化推进，盲端死胡同因无净通量流经而完全萎缩消失，仅留下唯一的全局最短路径。
3. **⚖️ 布雷斯悖论消解 (Braess's Paradox)**：
   - 经典的四节点对称双路拓扑，中间架设一条看似极速的交叉捷径；
   - 交互式观察黏菌系统如何在全局能耗优化下自主权衡捷径承载力。
4. **🏭 工场自愈物流环 (Fault-Tolerant Factory Ring)**：
   - 包含原料进料、高炉冶炼、芯片制造、核电储能、精密装配与发射交付 6 大工业集群；
   - 支持使用手术刀切断任意骨干管线，验证备用旁路在 3 秒内的自发膨胀与自愈接管。
5. **🧪 自由培养皿 (Freeform Sandbox)**：
   - 零初始约束，支持在画布任意位置播种营养孢子，实时扩展网络。

---

## 🛠️ 交互工具与原生 Web Audio 声学

- **全键盘与鼠标交互**：
  - `[1]` **检视 / 移位**：拖拽任意节点实时重塑几何拓扑；
  - `[2]` **投喂建站**：在空白处点击播种新站点，自动接入临近候选管网；
  - `[3]` **手术刀**：精准点击切断指定管道，诱发断管急救演化；
  - `[4]` **洪峰压力测试**：瞬间激增 250% 通量负载，检验网络吞吐弹性；
  - `[Space]` **暂停 / 继续** 演化。
- **原生 Web Audio 动态声学合成**：
  - **Sub-Bass 55Hz 生物原浆基频**：频率随全网总导通通量动态微调；
  - **断管急切剪断音**：高频锯齿波向下急跌的穿透性反馈；
  - **洪峰激波共振音**：三角波扫描过载音效。

---

## 🧪 自动化测试套件 (Python)

项目配备了严格的物理动力学与拓扑数学单元测试，覆盖 15 项关键力学机制：

```bash
PYTHONPATH=. uv run --with pytest pytest tests/
```

- `test_node_and_edge_creation`：节点与边缘拓扑构建与阻抗折算
- `test_poiseuille_pressure_and_flux`：泊肃叶通量与压力单调性
- `test_kirchhoff_mass_conservation`：中继节点基尔霍夫质量守恒
- `test_adaptation_superlinear_growth_and_decay`：超线性正反馈与低通量萎缩
- `test_parallel_shortest_path_geodesic_selection`：两点平行长短路径测地线竞争选择（>15倍优势）
- `test_sever_edge_and_fault_reroute`：断管手术刀切断与备用管线自愈重路由
- `test_braess_paradox_network`：布雷斯拓扑压力与平衡通量
- `test_tokyo_network_stability`：东京城铁多步演化数值稳定性（无 NaN，导通度下界保护）
- `test_factory_ring_self_healing`：工业六节点环网断管自愈验证
- `test_energy_dissipation_decrease`：Murray 能耗与网络总导通度收敛
- `test_invalid_inputs`：自环与未知节点边界校验
- `test_sever_non_existent_edge`：无效切断容错
- `test_custom_conductance_bounds`：大通量下导通度刚性上限与下限保护
- `test_pressure_relaxation_empty_network`：空网络鲁棒性
- `test_maze_channel_convergence`：迷宫通道死胡同管线萎缩收敛

---

## 🚀 启动与体验

- **在线即玩**：[https://hoshi0919.github.io/physarum-transit/](https://hoshi0919.github.io/physarum-transit/)
- **本地运行**：直接使用任意现代浏览器打开 `index.html` 即可，零外部依赖，纯原生 Canvas 与 Web Audio。
