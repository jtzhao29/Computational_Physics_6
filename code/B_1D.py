import numpy as np
import matplotlib.pyplot as plt

def onr_D_random_walk(t:int)->np.ndarray:
    position = np.zeros(t+1)
    choices = [-1,1]
    position[0] = 0
    for i in range(t):
        position[i+1] = position[i] + np.random.choice(choices)
    return position

def calculate_x_abs_mean_one_D(t0:int,tf: int,N_particle:int) -> np.ndarray:
    time_ndarray = np.arange(t0, tf)
    x_abs_mean = np.zeros(tf)
    position = np.zeros(tf,N_particle) 
    position[0,:] = 0
    x_abs_mean[0] = np.sum(np.abs(position[0,:]))
    for time in range(tf):
        for particle in range(N_particle):
            position[time+1,particle] = position[time,particle] + np.random.choice([-1,1])
        x_abs_mean[time+1] = np.sum(np.abs(position[time+1,:]))
    return x_abs_mean/N_particle,time_ndarray

def plot_x_abs_mean(x_mean:np.ndarray,t:int):
    plt.plot(x_mean)
    plt.xlabel("time")
    plt.ylabel("x_abs_mean")
    plt.title("x_abs_mean vs time")
    plt.show()
    


