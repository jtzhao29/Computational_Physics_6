import numpy as np
import matplotlib.pyplot as plt


def V(q:float)->float:
    """
    势能
    """
    return 0.5 * q**2

def dV_dq(q:float)->float:
    """
    势能的导数
    """
    return q

# Heun method for Langevin dynamics
def heun_langevin(q0:float, v0:float)->tuple[np.ndarray,np.ndarray]:
    """
    heun算法
    """
    q = np.zeros((N_particles, N_steps + 1)) #粒子数-时间数
    v = np.zeros((N_particles, N_steps + 1))
    q[:, 0] = q0
    v[:, 0] = v0

    for i in range(N_steps):
        dW = np.random.normal(0, np.sqrt(dt), N_particles)
        
        v_tilde = v[:, i] + (-v[:, i] - dV_dq(q[:, i])) * dt + D * dW
        q_tilde = q[:, i] + v[:, i] * dt
        
        # Corrector step
        v[:, i+1] = v[:, i] + 0.5 * ((-v[:, i] - dV_dq(q[:, i])) + (-v_tilde - dV_dq(q_tilde))) * dt + D * dW
        q[:, i+1] = q[:, i] + 0.5 * (v[:, i] + v_tilde) * dt

    return q, v
def heun_for_kinetic_and_potential(q0:float, v0:float,N_particles:int)->tuple[np.ndarray,np.ndarray]:
    """
    heun算法
    """
    q = np.zeros((N_particles, N_steps + 1)) #粒子数-时间数
    v = np.zeros((N_particles, N_steps + 1))
    q[:, 0] = q0
    v[:, 0] = v0

    kinetic = np.zeros( N_steps + 1)
    potential = np.zeros( N_steps + 1)
    kinetic[0] = 0.5 * v0**2
    potential[0] = 0.5 * q0**2

    for i in range(N_steps):
        dW = np.random.normal(0, np.sqrt(dt), N_particles)
        
        v_tilde = v[:, i] + (-v[:, i] - dV_dq(q[:, i])) * dt + D * dW
        q_tilde = q[:, i] + v[:, i] * dt
        
        # Corrector step
        v[:, i+1] = v[:, i] + 0.5 * ((-v[:, i] - dV_dq(q[:, i])) + (-v_tilde - dV_dq(q_tilde))) * dt + D * dW
        q[:, i+1] = q[:, i] + 0.5 * (v[:, i] + v_tilde) * dt
        kinetic[i+1] = np.sum(0.5*v[:,i+1]**2)/N_particles
        potential[i+1] = np.sum(0.5*q[:,i+1]**2)/N_particles


    return kinetic,potential


# Plotting
def plot_kinetic_and_potential_energy(E_k_over_time_A:np.ndarray,V_over_time_A:np.ndarray,E_k_over_time_B:np.ndarray,V_over_time_B:np.ndarray):
    """
    画图的函数
    """
    plt.figure(figsize=(10, 6))
    plt.plot(t, E_k_over_time_A, label='Kinetic Energy (A)', color='blue')
    plt.plot(t, V_over_time_A, label='Potential Energy (A)', color='blue', linestyle='--')
    plt.plot(t, E_k_over_time_B, label='Kinetic Energy (B)', color='red')
    plt.plot(t, V_over_time_B, label='Potential Energy (B)', color='red', linestyle='--')
    plt.xlabel('Time $t$',fontsize=20)
    plt.ylabel('Energy',fontsize=20)
    plt.title('Evolution of Ensemble-Averaged Kinetic and Potential Energy',fontsize=20)
    plt.legend()
    plt.grid(True)
    # plt.tight_layout()
    path = f"./figures/A_1.png"
    plt.savefig(path, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    
    T_total = 10.0    
    dt = 0.01         
    N_steps = int(T_total / dt)
    N_particles = 10000  
    D = np.sqrt(2)      
    t = np.linspace(0, T_total, N_steps + 1)

    # q_A, v_A = heun_langevin(q0=0, v0=1)
    # q_B, v_B = heun_langevin(q0=4, v0=0)
        
    # E_k_A = 0.5 * np.mean(v_A**2, axis=0)
    # V_A = 0.5 * np.mean(q_A**2, axis=0)

    # E_k_B = 0.5 * np.mean(v_B**2, axis=0)
    # V_B = 0.5 * np.mean(q_B**2, axis=0)

    E_k_A_time, v_A_time = heun_for_kinetic_and_potential(0, 1,N_particles)
    E_k_B_time, V_B_time = heun_for_kinetic_and_potential(4, 0,N_particles)

    plot_kinetic_and_potential_energy(E_k_A_time,v_A_time,E_k_B_time,V_B_time)



