import numpy as np
import matplotlib.pyplot as plt

# functions needed:
# 1. mean of internodal distance to nearest 3 neighbours (average)
# 2. std deviation of internodal distance to nearest 3 neighbours (average)
# 3. mean of (d_max - d_min)

def stats(nodes):
    avgs_di = np.array([])
    differences = np.array([])
    min_di = np.array([])
    for i in range(len(nodes)):
        remaining_nodes = np.delete(nodes, [i], axis=0)
        distance_xy_sq = ((remaining_nodes) - (nodes[i]))**2
        distance_sq = np.array(distance_xy_sq[:,0] + distance_xy_sq[:,1])
        distance = np.sqrt(distance_sq)
        # arg_three_nearest = np.argpartition(distance, 3)
        # di = distance[arg_three_nearest] # this is just sorted
        di = np.sort(distance)
        avg_di = (di[0] + di[1] + di[2])/3
        avgs_di = np.append(avgs_di, avg_di)

        difference = np.abs(di[2] - di[0])
        differences = np.append(differences, difference)

        min_di = np.append(min_di, di[0])
    
        
    # Creating histogram
    fig, ax = plt.subplots(figsize =(10, 7))
    # ax.hist(min_di, bins = np.linspace(0.03, 0.06, 21), color= "darkturquoise", edgecolor= "black") # fornberg perfect
    # ax.hist(min_di, bins = np.linspace(0.01, 0.06, 21), color= "darkturquoise", edgecolor= "black")
    ax.hist(min_di, bins = np.linspace(0.005, 0.04, 21), color= "darkturquoise", edgecolor= "black")

    # Show plot
    plt.show()

    mean = np.average(avgs_di)
    std_var = np.std(avgs_di)
    diff = np.average(differences)


    return mean, std_var, diff
