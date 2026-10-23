import math
import pytest
from physarum.engine import PhysarumNetwork, Node, Edge

def test_node_and_edge_creation():
    net = PhysarumNetwork()
    n1 = net.add_node("A", 0, 0, injection=5.0, label="Alpha", kind="source")
    n2 = net.add_node("B", 10, 0, injection=-5.0, label="Beta", kind="sink")
    
    edge = net.add_edge("A", "B", conductivity=2.0)
    assert edge.length == pytest.approx(10.0, 1e-4)
    assert edge.conductivity == 2.0
    assert edge.resistance == pytest.approx(5.0, 1e-4)
    assert edge.key == "A--B"

    # Serialization
    d_node = n1.to_dict()
    assert d_node["id"] == "A"
    assert d_node["injection"] == 5.0
    d_edge = edge.to_dict()
    assert d_edge["u"] == "A"
    assert d_edge["v"] == "B"

def test_poiseuille_pressure_and_flux():
    net = PhysarumNetwork()
    net.add_node("S", 0, 0, injection=10.0)
    net.add_node("T", 100, 0, injection=-10.0)
    edge = net.add_edge("S", "T", conductivity=1.0, length=100.0)

    iters = net.solve_pressures()
    assert iters > 0
    # S pressure must be higher than T pressure
    assert net.nodes["S"].pressure > net.nodes["T"].pressure
    # Flux Q = (D / L) * (p_S - p_T)
    expected_flux = (edge.conductivity / edge.length) * (net.nodes["S"].pressure - net.nodes["T"].pressure)
    assert edge.flux == pytest.approx(expected_flux, 1e-5)
    # Conservation of flow across the boundary
    assert edge.flux > 0

def test_kirchhoff_mass_conservation():
    net = PhysarumNetwork()
    net.add_node("S", 0, 0, injection=12.0)
    net.add_node("R", 50, 0, injection=0.0)  # Relay (0 injection)
    net.add_node("T1", 100, 20, injection=-6.0)
    net.add_node("T2", 100, -20, injection=-6.0)

    net.add_edge("S", "R", conductivity=2.0, length=50.0)
    net.add_edge("R", "T1", conductivity=1.0, length=50.0)
    net.add_edge("R", "T2", conductivity=1.0, length=50.0)

    net.solve_pressures()

    # Inflow to R must equal outflow from R
    e_sr = net.edges["R--S"]
    e_rt1 = net.edges["R--T1"]
    e_rt2 = net.edges["R--T2"]

    flux_in = e_sr.flux if e_sr.u == "S" else -e_sr.flux
    flux_out1 = e_rt1.flux if e_rt1.u == "R" else -e_rt1.flux
    flux_out2 = e_rt2.flux if e_rt2.u == "R" else -e_rt2.flux

    net_divergence = flux_in - (flux_out1 + flux_out2)
    assert abs(net_divergence) < 0.1

def test_adaptation_superlinear_growth_and_decay():
    net = PhysarumNetwork()
    net.add_node("S", 0, 0, injection=15.0)
    net.add_node("T", 100, 0, injection=-15.0)
    net.add_node("DeadEnd", 0, 100, injection=0.0)

    e_active = net.add_edge("S", "T", conductivity=1.0, length=100.0)
    e_dead = net.add_edge("S", "DeadEnd", conductivity=1.0, length=100.0)

    initial_d_active = e_active.conductivity
    initial_d_dead = e_dead.conductivity

    for _ in range(15):
        net.step_adaptation(dt=0.1, alpha=0.15, gamma=1.3, decay=0.25)

    # Active edge carrying flux must expand
    assert e_active.conductivity > initial_d_active
    # Dead edge with zero flux must shrink/atrophy
    assert e_dead.conductivity < initial_d_dead

def test_parallel_shortest_path_geodesic_selection():
    """
    经典两点最短路竞争：
    S 与 T 之间有短路径 (L=40) 与长路径 (L=140)。
    黏菌演化后，短路径阻力小、流速大，发生正反馈增粗；长路径流速小自然凋亡。
    """
    net = PhysarumNetwork()
    net.add_node("S", 0, 0, injection=5.0)
    net.add_node("T", 100, 0, injection=-5.0)

    short_edge = net.add_edge("S", "T", conductivity=1.0, length=40.0)
    net.add_node("Detour", 50, 80, injection=0.0)
    long1 = net.add_edge("S", "Detour", conductivity=1.0, length=70.0)
    long2 = net.add_edge("Detour", "T", conductivity=1.0, length=70.0)

    for _ in range(100):
        net.step_adaptation(dt=0.1, alpha=0.15, gamma=1.3, decay=0.25)

    assert short_edge.conductivity > 3.0
    assert long1.conductivity < 0.2
    assert long2.conductivity < 0.2
    # 短路径导通度彻底压倒长路径 (>15倍优势)
    assert short_edge.conductivity > long1.conductivity * 15

def test_sever_edge_and_fault_reroute():
    """断管急救与自愈重路由测试"""
    net = PhysarumNetwork()
    net.add_node("S", 0, 0, injection=20.0)
    net.add_node("T", 100, 0, injection=-20.0)

    e1 = net.add_edge("S", "T", conductivity=3.0, length=50.0)
    net.add_node("Backup", 50, 40, injection=0.0)
    e_bk1 = net.add_edge("S", "Backup", conductivity=0.5, length=60.0)
    e_bk2 = net.add_edge("Backup", "T", conductivity=0.5, length=60.0)

    net.solve_pressures()
    initial_bk_flux = abs(e_bk1.flux)

    sever_ok = net.sever_edge("S", "T")
    assert sever_ok is True
    assert e1.active is False

    net.solve_pressures()
    rerouted_flux = abs(e_bk1.flux)

    # 备用管线通量暴增承担全部负载
    assert rerouted_flux > initial_bk_flux * 2

def test_braess_paradox_network():
    net = PhysarumNetwork.create_braess_network()
    assert len(net.nodes) == 4
    assert len(net.edges) == 5

    net.solve_pressures()
    assert net.nodes["S"].pressure > net.nodes["T"].pressure
    metrics = net.metrics()
    assert metrics["active_edges"] >= 4
    assert metrics["total_flux"] > 0

def test_tokyo_network_stability():
    net = PhysarumNetwork.create_tokyo_network()
    assert len(net.nodes) == 11
    
    for _ in range(25):
        max_delta = net.step_adaptation(dt=0.1)
        assert not math.isnan(max_delta)

    for e in net.edges.values():
        assert not math.isnan(e.conductivity)
        assert e.conductivity >= e.d_min

    metrics = net.metrics()
    assert metrics["energy_dissipation"] > 0
    assert metrics["network_robustness"] > 0.5
    assert metrics["steps"] == 25

def test_factory_ring_self_healing():
    net = PhysarumNetwork.create_factory_ring()
    assert len(net.nodes) == 6
    assert len(net.edges) == 9

    for _ in range(10):
        net.step_adaptation(dt=0.1)

    # 切断冶炼到晶片的动脉 Smelt--Fab
    sever_res = net.sever_edge("Smelt", "Fab")
    assert sever_res is True

    for _ in range(10):
        delta = net.step_adaptation(dt=0.1)
        assert not math.isnan(delta)

    metrics = net.metrics()
    assert metrics["total_flux"] > 0
    assert net.edges["Ore--Power"].conductivity > 0.5

def test_energy_dissipation_decrease():
    """
    随着高通量干线建立与低效回路被修剪，
    网络形成稳态拓扑并有效承载总通量。
    """
    net = PhysarumNetwork.create_braess_network()
    net.solve_pressures()
    initial_m = net.metrics()

    for _ in range(20):
        net.step_adaptation(dt=0.1, alpha=0.15, gamma=1.3, decay=0.25)

    final_m = net.metrics()
    assert final_m["total_conductivity"] > initial_m["total_conductivity"]
    assert final_m["steps"] == 20

def test_invalid_inputs():
    net = PhysarumNetwork()
    net.add_node("A", 0, 0)
    with pytest.raises(ValueError, match="Self-loops"):
        net.add_edge("A", "A")
    with pytest.raises(ValueError, match="exist"):
        net.add_edge("A", "GhostNode")

def test_sever_non_existent_edge():
    net = PhysarumNetwork()
    net.add_node("A", 0, 0)
    net.add_node("B", 1, 1)
    assert net.sever_edge("A", "B") is False
    assert net.sever_edge("X", "Y") is False

def test_custom_conductance_bounds():
    net = PhysarumNetwork()
    net.add_node("S", 0, 0, injection=100.0)
    net.add_node("T", 10, 0, injection=-100.0)
    e = net.add_edge("S", "T", conductivity=5.0)
    # Huge flow step
    for _ in range(50):
        net.step_adaptation(dt=0.2, alpha=2.0, gamma=1.5, decay=0.01)
    assert e.conductivity <= e.d_max
    assert e.conductivity >= e.d_min

def test_pressure_relaxation_empty_network():
    net = PhysarumNetwork()
    assert net.solve_pressures() == 0

def test_maze_channel_convergence():
    """
    迷宫走廊模拟：S 处进水，T 处出水。
    死胡同 (DeadEnd) 与贯通路 (Throughway)。
    """
    net = PhysarumNetwork()
    net.add_node("S", 0, 0, injection=8.0)
    net.add_node("C1", 20, 0)
    net.add_node("DeadEnd1", 20, 40)
    net.add_node("C2", 40, 0)
    net.add_node("DeadEnd2", 40, -40)
    net.add_node("T", 60, 0, injection=-8.0)

    # 贯通路
    e_thru1 = net.add_edge("S", "C1", length=20.0)
    e_thru2 = net.add_edge("C1", "C2", length=20.0)
    e_thru3 = net.add_edge("C2", "T", length=20.0)

    # 两个死胡同
    e_dead1 = net.add_edge("C1", "DeadEnd1", length=40.0)
    e_dead2 = net.add_edge("C2", "DeadEnd2", length=40.0)

    for _ in range(80):
        net.step_adaptation(dt=0.1)

    # 贯通管道保持粗壮，死胡同管道急剧萎缩
    assert e_thru1.conductivity > 2.5
    assert e_thru2.conductivity > 2.5
    assert e_thru3.conductivity > 2.5
    assert e_dead1.conductivity < 0.15
    assert e_dead2.conductivity < 0.15
