import numpy as np
import matplotlib.pyplot as plt
import os

# 创建图像保存文件夹
os.makedirs("./figures", exist_ok=True)

# 模拟参数
T = 1000  # 最大时间步
N = 10000  # 粒子数量

def simulate_comb_walk():
    x_positions = np.zeros(N, dtype=int)
    y_positions = np.zeros(N, dtype=int)
    avg_abs_x = []

    for t in range(1, T+1):
        for i in range(N):
            x, y = x_positions[i], y_positions[i]
            if y == 0:
                direction = np.random.choice(['up', 'down', 'left', 'right'])
            else:
                direction = np.random.choice(['up', 'down'])

            if direction == 'up':
                y_positions[i] += 1
            elif direction == 'down':
                y_positions[i] -= 1
            elif direction == 'left':
                x_positions[i] -= 1
            elif direction == 'right':
                x_positions[i] += 1

        avg_abs_x.append(np.mean(np.abs(x_positions)))

    return avg_abs_x

def simulate_1d_walk():
    positions = np.zeros(N, dtype=int)
    avg_abs_x = []

    for t in range(1, T+1):
        steps = np.random.choice([-1, 1], size=N)
        positions += steps
        avg_abs_x.append(np.mean(np.abs(positions)))

    return avg_abs_x

# 模拟并绘图
comb_x = simulate_comb_walk()
walk1d_x = simulate_1d_walk()

# 图1：梳子结构随机行走
plt.figure(figsize=(8, 6))
plt.plot(range(1, T+1), comb_x, label='Comb Lattice')
plt.xlabel("Time", fontsize=20)
plt.ylabel(r"$\langle |x| \rangle$", fontsize=20)
plt.title("Average Horizontal Displacement on Comb", fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig("./figures/comb_walk_avg_abs_x.png")
plt.close()

# 图2：与一维对比
plt.figure(figsize=(8, 6))
plt.plot(range(1, T+1), comb_x, label="Comb Lattice")
plt.plot(range(1, T+1), walk1d_x, label="1D Random Walk")
plt.xlabel("Time", fontsize=20)
plt.ylabel(r"$\langle |x| \rangle$", fontsize=20)
plt.title("Comparison: Comb vs 1D Random Walk", fontsize=20)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig("./figures/comparison_1d_vs_comb.png")
plt.close()