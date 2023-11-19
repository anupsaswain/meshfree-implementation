import numpy as np
from bayona_method import *
import time

# circle
def curve_equation_circle1(x, y):
    return x*x + y*y - 1

def curve_by_theta1(thetas, r):
    x0 = 0
    y0 = 0

    xs = x0 + r*np.cos(thetas)
    ys = y0 + r*np.sin(thetas)

    points = np.zeros([len(thetas), 2])
    points[:, 0] = xs
    points[:, 1] = ys
    return points

# circle:
def circle_boundary(curve_by_theta, r, h):
    perimeter = 2*np.pi*r
    n = int(perimeter/h) # number of points in boundary

    thetas = np.linspace(0, 2*np.pi, n)
    thetas = thetas[:-1]

    boundary_nodes = curve_by_theta(thetas, r)

    return boundary_nodes

def ghost_nodes(curve_by_theta, boundary_nodes, h, layers, r):
    ghost_nodes = boundary_nodes
    for i in range(layers):
        r_layer = r + h*(i+1)
        perimeter_layer = 2*np.pi*r_layer

        n_layer = int(perimeter_layer/h)

        layer_thetas = np.linspace(0, 2*np.pi, n_layer)
        layer_thetas = layer_thetas[:-1]
        
        boundary_layer = curve_by_theta(layer_thetas, r_layer)

        ghost_nodes = np.append(ghost_nodes, boundary_layer, axis=0)

    return ghost_nodes


if __name__ == "__main__":
    start = time.time()
    # circle
    r = 1
    circle_area = np.pi*r*r

    h = .025
    layers = 3 # ghost node layers
    n_iterations = 100
    n_proximity = 20
    tolerance = 1e-11

    # breakpoint()

    boundary_nodes2 = circle_boundary(curve_by_theta1, r, h)
    # plot_scatter(boundary_nodes2)
    ghost_nodes2 = ghost_nodes(curve_by_theta1, boundary_nodes2, h, layers, r)
    domain_nodes2 = discretize_domain(boundary_nodes2, h)
    frame_nodes2, interior_nodes2 = separate_frame_nodes(curve_equation_circle1, domain_nodes2)
    # frame_nodes2 = clean_frame(boundary_nodes2, frame_nodes2)
    # plot_scatter(frame_nodes2)
    thick_boundary2 = np.append(boundary_nodes2, ghost_nodes2, axis=0)
    # plot_scatter(thick_boundary2)
    fixed_points = np.append(frame_nodes2, thick_boundary2, axis=0)
    # plot_scatter(interior_nodes2)
    interior_nodes2 = clear_boundary(interior_nodes2, thick_boundary2, h)
    # plot_scatter(interior_nodes2)
    # interior_nodes_adjusted2 = adjust(interior_nodes2, fixed_points, n_iterations, n_proximity, h, tolerance)
    interior_nodes_adjusted2 = adjust(interior_nodes2, boundary_nodes2, n_iterations, n_proximity, h, tolerance)
    # interior_nodes_adjusted2 = adjust_mod1(interior_nodes2, fixed_points, n_iterations, n_proximity)

    # point to note: domain nodes = frame nodes + interior nodes

    frame_nodes2, interior_nodes_adjusted2 = separate_frame_nodes(curve_equation_circle1, interior_nodes_adjusted2)

    temp_points = np.append(boundary_nodes2, interior_nodes_adjusted2, axis=0)
    # print(len(temp_points))
    end = time.time()
    plot_scatter(temp_points)
    print("Time taken: " + str(end - start) + "s")

    mean, std_var, diff = stats(temp_points)
    print(mean)
    print(std_var)
    print(diff)