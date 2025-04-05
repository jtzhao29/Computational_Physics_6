import numpy as np
import matplotlib.pyplot as plt
import random as rand

## 先设定可选的移动组合
choices = [
    np.array([0, 1]),
    np.array([0, -1]),
    np.array([1, 0]),
    np.array([-1, 0])
]

def move(position:np.ndarray,choices:list)->np.ndarray:
    if position[1]==0:
        choice = rand.randint(0,4)
        position+=choices[choice]
    else:
        choice = rand.randint(0,2)
        position+=choices[choice]
    return position


def partical_random_walk(t:int)->np.ndarray:
    position = np.ndarray[0,0]
    time_position = np.zeros(t,2)
    for i in range(t):
        position = move(position,choices)
        time_position[i,:] = position
    return time_position

def calculate_x_abs_mean(t0:int,tf: int, num: int, choices: list) -> np.ndarray:
    time_ndarray = np.arange(t0, tf)
    x_mean = np.zeros(tf)
    for j in range(num):
        position = np.array([0, 0])  # 修复 np.ndarray 的错误用法
        for i in range(tf):
            x_mean[i] += abs(position[0])  # 修复索引访问错误
            if position[1] == 0:
                choice = rand.randint(0, 3)  # 修复索引范围错误
                position += choices[choice]
            else:
                choice = rand.randint(0, 1)
                position += choices[choice]
    return x_mean / num


def plot_x_abs_mean_log(x_mean:np.ndarray):
    plt.plot(x_mean)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log(t)',fontsize=20)
    plt.ylabel(r'log(<|x|>)',fontsize=20)
    plt.title(r'resolution of <|x|>',fontsize=20)
    plt.grid(True)
    path = r"./figures/B_log.png"
    plt.savefig(path)
    
    plt.show()


if __name__ == '__main__':
    t = 1000
    num = 10000
    x_mean = calculate_x_abs_mean(0,t,num,choices)
    # plt.plot(x_mean)
    # plt.xlabel('time',fontsize=20)
    # plt.ylabel('x_mean',fontsize=20)
    # plt.title(r'resolution of <|x|>',fontsize=20)
    # plt.grid(True)
    # path = r"./figures/B.png"
    # plt.savefig(path)
    # plt.show()
    plot_x_abs_mean_log(x_mean)
