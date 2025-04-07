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
    x_mean = x_mean[t0:tf]
    return x_mean / num,time_ndarray


def calculate_x_abs_mean_one_D(t0:int,tf: int,N_particle:int) -> np.ndarray:
    time_ndarray = np.arange(tf)
    x_abs_mean = np.zeros(tf)
    position = np.zeros((tf,N_particle) )
    position[0,:] = 0
    x_abs_mean[0] = np.sum(np.abs(position[0,:]))
    for time in range(tf-1):
        for particle in range(N_particle):
            position[time+1,particle] = position[time,particle] + np.random.choice([-1,1])
        x_abs_mean[time+1] = np.sum(np.abs(position[time+1,:]))
    x_abs_mean = x_abs_mean [t0:tf]
    return x_abs_mean/N_particle,time_ndarray[t0:tf]

def plot_x_abs_mean_log(x_mean:np.ndarray,t:int):
    plt.plot(x_mean)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log(t)',fontsize=20)
    plt.ylabel(r'log(<|x|>)',fontsize=20)
    plt.title(r'resolution of <|x|>',fontsize=20)
    plt.grid(True)
    path = f"./figures/B_log_t={t}.png"
    plt.savefig(path, bbox_inches='tight')
    
    plt.show()

def plot_x_abs_mean(x_mean:np.ndarray,t:int):
    plt.plot(x_mean)
    plt.xlabel('time',fontsize=20)
    plt.ylabel('x_mean',fontsize=20)
    plt.title(r'resolution of <|x|>',fontsize=20)
    plt.grid(True)
    path = f"./figures/B_t={t}.png"
    plt.savefig(path, bbox_inches='tight')
    plt.show()

# def plot_x_abs_mean_log_fit(x_mean:np.ndarray,t:int,time_ndarry:np.ndarray):
#     x_mean = np.log(x_mean)
#     time_ndarry = np.log(time_ndarry)
#     plt.plot(time_ndarry,x_mean,label='data')
#     parmas = np.polyfit(time_ndarry,x_mean,1)
#     slope,intercept = parmas
#     plt.plot(time_ndarry,np.polyval(parmas,time_ndarry),'r',label =f'log<|x|> = {slope}*log(t)+{intercept}')
#     plt.xlabel('log(t)',fontsize=20)
#     plt.ylabel(r'log(<|x|>)',fontsize=20)
#     plt.title(r'resolution of <|x|>',fontsize=20)
#     plt.grid(True)
#     path = f"./figures/B_log_t={t}_fit.png"
#     plt.savefig(path, bbox_inches='tight')
#     plt.show()
    

def plot_x_abs_mean_log_fit(x_mean: np.ndarray, t: int,t0:int, time_ndarray: np.ndarray):
    # 对 x_mean 和 time_ndarray 取对数
    log_x_mean = np.log(x_mean)
    log_time_ndarray = np.log(time_ndarray)

    # 进行线性拟合
    params = np.polyfit(log_time_ndarray, log_x_mean, 1)
    slope, intercept = params

    # 绘制数据点和拟合直线
    plt.figure(figsize=(8, 6))
    plt.plot(log_time_ndarray, log_x_mean, 'o', label='Data (log-log)', markersize=5)
    plt.plot(log_time_ndarray, np.polyval(params, log_time_ndarray), 'r-', 
             label=f'Fit: log(<|x|>) = {slope:.6f} * log(t) + {intercept:.6f}')

    # 设置图形属性
    plt.xlabel('log(t)', fontsize=20)
    plt.ylabel(r'log(<|x|>)', fontsize=20)
    plt.title(r'Log-Log Fit of <|x|> vs t', fontsize=20)
    plt.legend(fontsize=12)
    plt.grid(True)

    # 保存图像
    path = f"./figures/B_log_from_{t0}_to_{t}_fit.png"
    plt.savefig(path, bbox_inches='tight')
    plt.show()
def plot_x_abs_mean_logfit(x_mean: np.ndarray, t: int, t0: int, time_ndarray: np.ndarray):
    # 对 time_ndarray 取对数
    log_time_ndarray = np.log(time_ndarray)

    # 进行线性拟合，拟合形式为 x = a * ln(t) + b
    params = np.polyfit(log_time_ndarray, x_mean, 1)
    a, b = params

    # 绘制数据点和拟合直线
    plt.figure(figsize=(8, 6))
    plt.plot(log_time_ndarray, x_mean, 'o', label='Data (log-linear)', markersize=5)
    plt.plot(log_time_ndarray, np.polyval(params, log_time_ndarray), 'r-', 
             label=f'Fit: x = {a:.6f} * ln(t) + {b:.6f}')

    # 设置图形属性
    plt.xlabel('ln(t)', fontsize=20)
    plt.ylabel(r'<|x|>', fontsize=20)
    plt.title(r'Log-Linear Fit of <|x|> vs ln(t)', fontsize=20)
    plt.legend(fontsize=12)
    plt.grid(True)

    # 保存图像
    path = f"./figures/B_loglinear_from_{t0}_to_{t}_fit.png"
    plt.savefig(path, bbox_inches='tight')
    plt.show()

def plot_one_D(x_mean:np.ndarray,t:int,t0:int,time_ndarry:np.ndarray):
    plt.plot(time_ndarry,x_mean)
    plt.xlabel('time',fontsize=20)
    plt.ylabel('x_mean',fontsize=20)
    plt.title(r'resolution of <|x|> in 1D random walk',fontsize=20)
    plt.grid(True)
    path = f"./figures/B_t={t}_from_{t0}_1D.png"
    plt.savefig(path, bbox_inches='tight')
    plt.show()


def plot_one_D_fit(x_mean:np.ndarray,t:int,t0:int,time_ndarry:np.ndarray):
    x_mean_log = np.log(x_mean)
    time_ndarry_log = np.log(time_ndarry)

    params = np.polyfit(time_ndarry_log, x_mean_log, 1)
    slope, intercept = params
    plt.figure(figsize=(8, 6))
    plt.plot(time_ndarry_log, x_mean_log, 'o', label='Data', markersize=4)
    plt.plot(time_ndarry_log, np.polyval(params, time_ndarry_log), 'r-', 
             label=f'Fit: log(<|x|>) = {slope:.6f} * log(t) + {intercept:.6f}') 
    plt.xlabel('log(t)', fontsize=20)
    plt.ylabel(r'log(<|x|>)', fontsize=20)
    plt.title(r'Log-Log Fit of <|x|> vs t in 1D random walk', fontsize=20)
    plt.legend(fontsize=12)
    plt.grid(True)
    path = f"./figures/B_log_from_{t0}_to_{t}_fit_1D.png"   
    plt.savefig(path, bbox_inches='tight')
    plt.show()

def plot_both(x_mean:np.ndarray,x_mean_1D:np.ndarray,t:int,t0:int,time_ndarry:np.ndarray):
    plt.plot(time_ndarry,x_mean,label='chomb random walk')
    plt.plot(time_ndarry_1D,x_mean_1D,label='1D random walk')
    plt.xlabel('time',fontsize=20)
    plt.ylabel('x_mean',fontsize=20)
    plt.title(r'resolution of <|x|>',fontsize=20)
    plt.grid(True)
    plt.legend()
    path = f"./figures/B_both_t={t}_from_{t0}_both.png"
    plt.savefig(path, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    t = 1000
    num = 10000
    t0=1
    x_mean,time_ndarry = calculate_x_abs_mean(t0,t,num,choices)
    # # plot_x_abs_mean(x_mean,t)
    # # plot_x_abs_mean_log(x_mean)
    # # print(x_mean[:20])
    plot_x_abs_mean_logfit(x_mean,t,t0,time_ndarry)
    x_mean_1D,time_ndarry_1D = calculate_x_abs_mean_one_D(t0,t,num)
    print(x_mean_1D[1])
    print(time_ndarry_1D[1])
    plot_one_D_fit(x_mean_1D,t,t0,time_ndarry_1D)
    plot_both(x_mean,x_mean_1D,t,t0,time_ndarry)
    