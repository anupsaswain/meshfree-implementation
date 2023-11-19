# imports
import numpy as np
import matplotlib.pyplot as plt
from fornberg_rectangle import *        

# input: boundary nodes directly? equation of shape? y = f(x) form? what if implicit?

# discretise boundary ("if needed, use root finding algorithms": fix later)

# input:
# h (distance bw nodes); equivalent to rbf for uniform
# shape of boundary

def convert_radial_to_cartesian(r, theta):
    x = r*np.cos(theta)
    y = r*np.sin(theta)

    point = np.zeros([len(x), 2])

    point[:, 0] = x
    point[:, 1] = y

    return point

def boundary_nodes_gen(r_of_theta, r_dash_theta, h):
    # take a string and place particles h distance apart. how to calculate perimeter? do we have x and y? 
    # start at theta = 0
    theta = 0
    thetas_array = np.array([])
    
    while theta < 2*np.pi:
        thetas_array = np.append(thetas_array, np.array([theta]), axis=0)
        dtheta = np.sqrt(h**2/(r_of_theta(theta)**2 + r_dash_theta(theta)**2))
        theta = theta + dtheta

    # after getting final list of thetas,
    rs = r_of_theta(thetas_array)

    nodes = convert_radial_to_cartesian(rs, thetas_array)

    return nodes

# discretise domain
def discretize_domain(boundary_nodes, h):
    left_end = np.min(boundary_nodes[:, 0])
    right_end = np.max(boundary_nodes[:, 0])
    upper_end = np.max(boundary_nodes[:, 1])
    lower_end = np.min(boundary_nodes[:, 1])

    h_lines = np.arange(lower_end - h*4, upper_end + h*4, h)
    v_lines = np.arange(left_end - h*4, right_end + h*4, h)

    domain_nodes = np.array([[v_lines[i], h_lines[j]] for i in range(len(v_lines)) for j in range(len(h_lines))])

    return domain_nodes

# find frame nodes
def separate_frame_nodes(curve_equation, domain_nodes):
    frame_nodes = np.zeros([0, 2])
    interior_nodes = np.zeros([0, 2])
    for node in domain_nodes:
        x_node = node[0]
        y_node = node[1]
        if curve_equation(x_node, y_node) > 0:
            frame_nodes = np.append(frame_nodes, [[x_node, y_node]], axis=0)
        else: interior_nodes = np.append(interior_nodes, [[x_node, y_node]], axis=0)
    
    return frame_nodes, interior_nodes
            
# removing unwanted frame nodes:
def clean_frame(boundary_nodes, frame_nodes): 
    mesh_length = radius(boundary_nodes[0], boundary_nodes[1])

    points_out_of_frame = np.zeros(0, dtype=int)

    for i in range(len(frame_nodes)):
        nodes_in_ref = boundary_nodes - frame_nodes[i]
        distance_sq = (nodes_in_ref[:, 0])**2 + (nodes_in_ref[:, 1])**2
        distance = np.sqrt(distance_sq)

        if np.min(distance) > mesh_length*4: 
            points_out_of_frame = np.append(points_out_of_frame, [i], axis=0)
    
    frame_nodes = np.delete(frame_nodes, points_out_of_frame, axis=0)

    return frame_nodes

# fix frame + boundary nodes

# add more layers to frame

# discard points too close to boundary
def clear_boundary(interior_nodes, thick_boundary_nodes, h):
    points_to_delete = np.zeros(0, dtype=int)

    for i in range(len(interior_nodes)):
        nodes_in_ref = thick_boundary_nodes - interior_nodes[i]
        distance_sq = (nodes_in_ref[:, 0])**2 + (nodes_in_ref[:, 1])**2
        distance = np.sqrt(distance_sq)

        if np.min(distance) < h*.1 : 
            points_to_delete = np.append(points_to_delete, [i], axis=0)
    
    interior_nodes = np.delete(interior_nodes, points_to_delete, axis=0)

    return interior_nodes

# move the rest by anti-gravity equation; n iterations

def adjust(interior_nodes, fixed_nodes, n_iterations, n_proximity, h, tolerance):
    no_of_nodes = len(interior_nodes)
    for a in range(n_iterations):                                                   # iterate the process n times                        
        for i in range(len(interior_nodes)):                                        # for every node,
            interior_nodes_temp = np.delete(interior_nodes, [i], axis=0)            # excluding the node itself
            force_nodes = np.append(interior_nodes_temp, fixed_nodes, axis=0)       # all nodes that _could_ affect node to be adjusted
            nodes_in_ref = force_nodes - interior_nodes[i]
            distance_sq = (nodes_in_ref[:, 0])**2 + (nodes_in_ref[:, 1])**2
            distance = np.sqrt(distance_sq)                                         # define distance from every node

            effective_points = np.argpartition(distance, n_proximity)
            effective_points = effective_points[:n_proximity]
            force_sum = np.array([0, 0])

            if a >=1 and np.min(distance) >= h*(1 - tolerance) and np.min(distance) <= h*(1 + tolerance): break

            for k in range(n_proximity):
                point_arg = effective_points[k]
                r = distance[point_arg]
                r_vector = interior_nodes[i] - force_nodes[point_arg]

                # force_eff = r_vector*.5/(a+1)
                force_eff = r_vector/r**3

                # force_value = np.sqrt(force_eff[0]**2 + force_eff[1]**2)
                # breakpoint()

                # if force_value < 1e-2 or force_value > 1e-1 : break
                force_sum = force_sum + force_eff
            
            force_sum_mod = np.sqrt(force_sum[0]**2 + force_sum[1]**2)

            if force_sum_mod == 0: force_direction = np.array([0, 0]) 
            else: force_direction = force_sum/force_sum_mod                  # force_direction is unit!
            interior_nodes[i] = interior_nodes[i] + force_direction*(.2*h + .3*h*np.random.random(1))
        # plot_scatter(interior_nodes)

    return interior_nodes

# different adjust function:
def adjust_mod1(interior_nodes, fixed_nodes, n_iterations, n_proximity):
    # as n_iterations increases, uniformity should increase.
    for a in range(n_iterations):                                               # iterate the process n times
        force_nodes = np.append(interior_nodes, fixed_nodes, axis=0)            # all nodes that _could_ affect node to be adjusted
        for i in range(len(interior_nodes)):                                    # for every node,
            nodes_in_ref = force_nodes - interior_nodes[i]
            distance_sq = (nodes_in_ref[:, 0])**2 + (nodes_in_ref[:, 1])**2     # define distance from every node
            distance = np.sqrt(distance_sq)

            effective_points = np.argpartition(distance, n_proximity)
            effective_points = effective_points[:n_proximity]

            force_sum = np.array([0, 0])

            for k in range(n_proximity):
                point_arg = effective_points[k]
                r = distance[point_arg]
                if r == 0: break
                r_vector = interior_nodes[i] - force_nodes[point_arg]

                force_eff = r_vector/r**3     

                force_sum = force_sum + force_eff    # use only the direction

            force_sum_mod = np.sqrt(force_sum[0]**2 + force_sum[1]**2)

            if force_sum_mod == 0: force_direction = np.array([0, 0])
            else: force_direction = force_sum/force_sum_mod

            interior_nodes[i] = interior_nodes[i] + force_direction*.001

    return interior_nodes


# keep interior + boundary. discard frame

# plot
def plot_scatter(dot_locations):
    plt.scatter(dot_locations[:, 0], dot_locations[:, 1], s = 0.6)
    plt.show()

