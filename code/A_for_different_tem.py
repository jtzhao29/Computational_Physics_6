import numpy as np
import matplotlib.pyplot as plt
import os

def V(q: np.ndarray) -> np.ndarray:
    """势能 V(q) = 0.5 * q^2"""
    return 0.5 * q**2

def dV_dq(q: np.ndarray) -> np.ndarray:
    """势能的导数 dV/dq = q"""
    return q

def heun_for_kinetic_and_potential(q0: float, v0: float, D: float, N_particles: int, N_steps: int, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Heun方法模拟Langevin动力学，记录每一步的系统平均动能和势能
    """
    q = np.zeros((N_particles, N_steps + 1))
    v = np.zeros((N_particles, N_steps + 1))
    q[:, 0] = q0
    v[:, 0] = v0

    kinetic = np.zeros(N_steps + 1)
    potential = np.zeros(N_steps + 1)
    kinetic[0] = np.mean(0.5 * v0**2)
    potential[0] = np.mean(0.5 * q0**2)

    for i in range(N_steps):
        dW = np.random.normal(0, np.sqrt(dt), N_particles)

        v_tilde = v[:, i] + (-v[:, i] - dV_dq(q[:, i])) * dt + D * dW
        q_tilde = q[:, i] + v[:, i] * dt

        v[:, i+1] = v[:, i] + 0.5 * ((-v[:, i] - dV_dq(q[:, i])) + (-v_tilde - dV_dq(q_tilde))) * dt + D * dW
        q[:, i+1] = q[:, i] + 0.5 * (v[:, i] + v_tilde) * dt

        kinetic[i+1] = np.mean(0.5 * v[:, i+1]**2)
        potential[i+1] = np.mean(0.5 * q[:, i+1]**2)

    return kinetic, potential

def plot_energy_vs_time(t: np.ndarray, kinetic_A: np.ndarray, potential_A: np.ndarray,kinetic_B: np.ndarray, potential_B: np.ndarray, temperature: float):
    """
    画出在某一温度下，系统随时间演化的动能和势能
    """
    plt.figure(figsize=(10, 6))
    plt.plot(t, kinetic_A, label="Kinetic Energy of A", color='blue')
    plt.plot(t, potential_A, label="Potential Energy of B", color='blue', linestyle='--')
    plt.plot(t, kinetic_B, label="Kinetic Energy of B", color='red')
    plt.plot(t, potential_B, label="Potential Energy of B", color='red', linestyle='--')
    plt.xlabel('Time $t$', fontsize=16)
    plt.ylabel('Energy', fontsize=16)
    plt.title(f'Average Energy vs Time at T={temperature}', fontsize=18)
    plt.legend()
    plt.grid(True)
    os.makedirs('./figures', exist_ok=True)
    plt.savefig(f'./figures/energy_vs_time_T={temperature}_time={t[-1]}.png', bbox_inches='tight')
    plt.show()

def simulate_temperature(T: float, N_particles: int, N_steps: int, dt: float):
    """
    在给定温度 T 下进行模拟，并输出平均动能和势能
    """
    D = np.sqrt(2 * T)
    q0, v0 = 0.0, 1
    q0_1,v0_1 = 4.0,0.0
    kinetic_A, potential_A = heun_for_kinetic_and_potential(q0, v0, D, N_particles, N_steps, dt)
    kinetic_B, potential_B = heun_for_kinetic_and_potential(q0_1, v0_1, D, N_particles, N_steps, dt)

    # 只取最后一半时间的数据作为充分弛豫后的平衡态
    steady_k_A = np.mean(kinetic_A[N_steps // 2:])
    steady_v_A= np.mean(potential_A[N_steps // 2:])
    stead_k_B = np.mean(kinetic_B[N_steps // 2:])
    stead_v_B = np.mean(potential_B[N_steps // 2:])

    print(f"A:T={T:.2f}: <E_k>={steady_k_A:.4f}, <V>={steady_v_A:.4f}")
    print(f"B:T={T:.2f}: <E_k>={stead_k_B:.4f}, <V>={stead_v_B:.4f}")

    t = np.linspace(0, N_steps * dt, N_steps + 1)
    plot_energy_vs_time(t, kinetic_A, potential_A,kinetic_B,potential_B, temperature=T)

if __name__ == "__main__":
    T_total = 10.0
    dt = 0.01
    N_steps = int(T_total / dt)
    N_particles = 10000

    # 多个温度下分别模拟
    temperatures = [ 0.5,1.0,2.0,4.0,8.0,16.0]
    for T in temperatures:
        simulate_temperature(T, N_particles, N_steps, dt)
