import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# 新势能 V(x) = 0.5 * x^4 及其导数
# ----------------------------
def V(q: float) -> float:
    return 0.5 * q**4

def dV_dq(q: float) -> float:
    return 2 * q**3

def heun_for_kinetic_and_potential(q0: float, v0: float) -> tuple[np.ndarray, np.ndarray]:
    q = np.zeros((N_particles, N_steps + 1))
    v = np.zeros((N_particles, N_steps + 1))
    q[:, 0] = q0
    v[:, 0] = v0

    kinetic = np.zeros(N_steps + 1)
    potential = np.zeros(N_steps + 1)
    kinetic[0] = np.sum(0.5 * v[:, 0]**2) / N_particles
    potential[0] = np.sum(V(q[:, 0])) / N_particles

    for i in range(N_steps):
        dW = np.random.normal(0, np.sqrt(dt), N_particles)

        v_tilde = v[:, i] + (-v[:, i] - dV_dq(q[:, i])) * dt + D * dW
        q_tilde = q[:, i] + v[:, i] * dt

        v[:, i + 1] = v[:, i] + 0.5 * ((-v[:, i] - dV_dq(q[:, i])) + (-v_tilde - dV_dq(q_tilde))) * dt + D * dW
        q[:, i + 1] = q[:, i] + 0.5 * (v[:, i] + v_tilde) * dt

        kinetic[i + 1] = np.sum(0.5 * v[:, i + 1]**2) / N_particles
        potential[i + 1] = np.sum(V(q[:, i + 1])) / N_particles

    return kinetic, potential

def scan_temperature(T_list: list[float], q0: float = 0.0, v0: float = 0.0):
    avg_Ek = []
    avg_V = []

    for T in T_list:
        global D
        D = np.sqrt(2 * T)

        kinetic, potential = heun_for_kinetic_and_potential(q0, v0)

        start_idx = int(0.8 * len(kinetic))
        Ek_inf = np.mean(kinetic[start_idx:])
        V_inf = np.mean(potential[start_idx:])

        avg_Ek.append(Ek_inf)
        avg_V.append(V_inf)

    plt.figure(figsize=(8, 6))
    plt.plot(T_list, avg_Ek, 'o-', label=r'$\langle E_k(\infty) \rangle$')
    plt.plot(T_list, avg_V, 's--', label=r'$\langle V(\infty) \rangle$')
    plt.plot(T_list, [0.5*T for T in T_list], 'k--', label=r'$0.5 T$ ')
    plt.plot(T_list, [0.25*T for T in T_list], 'g--', label=r'$0.25 T$ ')
    plt.xlabel('Temperature $T$', fontsize=20)
    plt.ylabel('Average Energy', fontsize=20)
    plt.title('Quartic Potential: Long-Time Average Energy vs Temperature', fontsize=20)
    plt.legend()
    plt.grid(True)
    plt.savefig('./figures/temperature_scan_quartic.png', bbox_inches='tight')
    plt.show()

# ----------------------------
# 主程序入口
# ----------------------------
if __name__ == "__main__":
    T_total = 10.0
    dt = 0.01
    N_steps = int(T_total / dt)
    N_particles = 10000
    D = np.sqrt(2)  # 初始温度 T=1
    t = np.linspace(0, T_total, N_steps + 1)

    T_list = np.linspace(0.5, 5.0, 10)
    scan_temperature(T_list)
