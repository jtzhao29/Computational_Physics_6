import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# 势能函数与其导数
# ----------------------------
def V(q: float) -> float:
    return 0.5 * q**2

def dV_dq(q: float) -> float:
    return q

# ----------------------------
# Heun 方法模拟 Langevin 动力学，记录能量演化
# ----------------------------
def heun_for_kinetic_and_potential(q0: float, v0: float) -> tuple[np.ndarray, np.ndarray]:
    q = np.zeros((N_particles, N_steps + 1))
    v = np.zeros((N_particles, N_steps + 1))
    q[:, 0] = q0
    v[:, 0] = v0

    kinetic = np.zeros(N_steps + 1)
    potential = np.zeros(N_steps + 1)
    kinetic[0] = np.sum(0.5 * v[:, 0]**2) / N_particles
    potential[0] = np.sum(0.5 * q[:, 0]**2) / N_particles

    for i in range(N_steps):
        dW = np.random.normal(0, np.sqrt(dt), N_particles)

        v_tilde = v[:, i] + (-v[:, i] - dV_dq(q[:, i])) * dt + D * dW
        q_tilde = q[:, i] + v[:, i] * dt

        v[:, i + 1] = v[:, i] + 0.5 * ((-v[:, i] - dV_dq(q[:, i])) + (-v_tilde - dV_dq(q_tilde))) * dt + D * dW
        q[:, i + 1] = q[:, i] + 0.5 * (v[:, i] + v_tilde) * dt

        kinetic[i + 1] = np.sum(0.5 * v[:, i + 1]**2) / N_particles
        potential[i + 1] = np.sum(0.5 * q[:, i + 1]**2) / N_particles

    return kinetic, potential

# ----------------------------
# 画出 A 与 B 两个初态下能量演化
# ----------------------------
def plot_kinetic_and_potential_energy(E_k_over_time_A, V_over_time_A, E_k_over_time_B, V_over_time_B):
    plt.figure(figsize=(10, 6))
    plt.plot(t, E_k_over_time_A, label='Kinetic Energy (A)', color='blue')
    plt.plot(t, V_over_time_A, label='Potential Energy (A)', color='blue', linestyle='--')
    plt.plot(t, E_k_over_time_B, label='Kinetic Energy (B)', color='red')
    plt.plot(t, V_over_time_B, label='Potential Energy (B)', color='red', linestyle='--')
    plt.xlabel('Time $t$', fontsize=20)
    plt.ylabel('Energy', fontsize=20)
    plt.title('Evolution of Ensemble-Averaged Kinetic and Potential Energy', fontsize=20)
    plt.legend()
    plt.grid(True)
    plt.savefig('./figures/A_1.png', bbox_inches='tight')
    plt.show()

# ----------------------------
# 改变温度，计算充分弛豫后的平均能量
# ----------------------------
def scan_temperature(T_list: list[float], q0: float = 0.0, v0: float = 0.0):
    avg_Ek = []
    avg_V = []

    for T in T_list:
        global D
        D = np.sqrt(2 * T)

        kinetic, potential = heun_for_kinetic_and_potential(q0, v0)

        # 弛豫后取最后 20% 做时间平均
        start_idx = int(0.8 * len(kinetic))
        Ek_inf = np.mean(kinetic[start_idx:])
        V_inf = np.mean(potential[start_idx:])

        avg_Ek.append(Ek_inf)
        avg_V.append(V_inf)

    # 作图
    plt.figure(figsize=(8, 6))
    plt.plot(T_list, avg_Ek, 'o-', label=r'$\langle E_k(\infty) \rangle$')
    plt.plot(T_list, avg_V, 's--', label=r'$\langle V(\infty) \rangle$')
    plt.plot(T_list, [0.5*T for T in T_list], 'k--', label=r'$0.5T$')
    plt.xlabel('Temperature $T$', fontsize=16)
    plt.ylabel('Average Energy', fontsize=16)
    plt.title('Long-Time Average Energy vs Temperature', fontsize=18)
    plt.legend()
    plt.grid(True)
    plt.savefig('./figures/temperature_scan.png', bbox_inches='tight')
    plt.show()

# ----------------------------
# 主程序入口
# ----------------------------
if __name__ == "__main__":
    # 参数设置
    T_total = 10.0
    dt = 0.01
    N_steps = int(T_total / dt)
    N_particles = 10000
    D = np.sqrt(2)  # 初始温度 T=1 的扩散强度
    t = np.linspace(0, T_total, N_steps + 1)

    # (a) 不同初态下能量演化
    E_k_A_time, V_A_time = heun_for_kinetic_and_potential(q0=0, v0=1)
    E_k_B_time, V_B_time = heun_for_kinetic_and_potential(q0=4, v0=0)
    plot_kinetic_and_potential_energy(E_k_A_time, V_A_time, E_k_B_time, V_B_time)

    # (b) 改变温度，分析充分弛豫后能量
    T_list = np.linspace(0.5, 5.0, 10)
    scan_temperature(T_list)
