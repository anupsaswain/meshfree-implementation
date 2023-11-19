import numpy as np
from bayona_method import *
from get_stats import stats

def rectangle_boundary(rectangle, n_per_unit_length):
    left_end  = rectangle.centre[0] - rectangle.length/2
    right_end = rectangle.centre[0] + rectangle.length/2
    upper_end = rectangle.centre[1] + rectangle.height/2
    lower_end = rectangle.centre[1] - rectangle.height/2

    n_horizontal = int((right_end - left_end)*n_per_unit_length)
    n_vertical = int((upper_end - lower_end)*n_per_unit_length)

    h_points_x_coords = np.linspace(left_end, right_end, n_horizontal)
    h_points_x_coords = h_points_x_coords[1:-1]
    h_points_y_coords_lower = lower_end*np.ones_like(h_points_x_coords)
    h_points_y_coords_upper = upper_end*np.ones_like(h_points_x_coords)

    v_points_y_coords = np.linspace(lower_end, upper_end, n_vertical)
    v_points_x_coords_left = left_end*np.ones_like(v_points_y_coords)
    v_points_x_coords_right = right_end*np.ones_like(v_points_y_coords)

    lower_line = np.zeros([len(h_points_x_coords), 2])
    lower_line[:, 0] = h_points_x_coords
    lower_line[:, 1] = h_points_y_coords_lower

    upper_line = np.zeros([len(h_points_x_coords), 2])
    upper_line[:, 0] = h_points_x_coords
    upper_line[:, 1] = h_points_y_coords_upper
    points = np.append(lower_line, upper_line, axis=0)

    left_line = np.zeros([len(v_points_y_coords), 2])
    left_line[:, 0] = v_points_x_coords_left
    left_line[:, 1] = v_points_y_coords
    points = np.append(points, left_line, axis=0)

    right_line = np.zeros([len(v_points_y_coords), 2])
    right_line[:, 0] = v_points_x_coords_right
    right_line[:, 1] = v_points_y_coords
    points = np.append(points, right_line, axis=0)

    # plot_scatter(points)

    return points

def curve_equation_rectangle1(x, y):
    a = (x-2.5)**2 - 6.25
    b = (y-1.5)**2 - 2.25

    if a > 0 and b > 0:
        return a*b
    elif a > 0 and b == 0: return a
    elif b > 0 and a == 0 : return b
    else: return -a*b

if __name__ == "__main__":
    
    h = 0.05
    a = 1e-11
    n_proximity = 20
    n_iterations = 100

    # rectangle
    box1 = Rectangle(5, 3, 2.5, 1.5)
    rectangle_area = 5*3
    boundary_nodes1 = rectangle_boundary(box1, int(1/h))
    domain_nodes1 = discretize_domain(boundary_nodes1, h)
    frame_nodes1, interior_nodes1 = separate_frame_nodes(curve_equation_rectangle1, domain_nodes1)
    # frame_nodes1 = clean_frame(boundary_nodes1, frame_nodes1)
    fixed_points = np.append(frame_nodes1, boundary_nodes1, axis=0)
    interior_nodes1 = clear_boundary(interior_nodes1, boundary_nodes1, h)
    interior_nodes_adjusted1 = adjust(interior_nodes1, fixed_points, n_iterations, n_proximity, h, a)

    # point to note: domain nodes = frame nodes + interior nodes

    temp_points = np.append(boundary_nodes1, interior_nodes_adjusted1, axis=0)
    # print(len(temp_points))
    plot_scatter(temp_points)

    mean, std_var, diff = stats(temp_points)
    print(mean)
    print(std_var)
    print(diff)