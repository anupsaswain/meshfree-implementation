import numpy as np
import matplotlib.pyplot as plt
from get_stats import stats

# for a rectangle with variable density function, density(x, y)
# rectangle definitions of length, width, centre

class Rectangle:
    def __init__(self, length, height, centre_x, centre_y):
        self.length = length
        self.height = height
        self.centre = np.array([centre_x, centre_y])
        self.left_end = self.centre[0] - self.length/2
        self.right_end = self.centre[0] + self.length/2


    # def boundary(self, n_per_unit_length):

# defining density function:
def density1(x, y):
    r = np.sqrt(x*x + y*y)
    e = np.e
    return (x+1)

def fornberg_distribution(density_function, rectangle, n_per_length):
    left_end  = rectangle.centre[0] - rectangle.length/2
    right_end = rectangle.centre[0] + rectangle.length/2
    upper_end = rectangle.centre[1] + rectangle.height/2
    lower_end = rectangle.centre[1] - rectangle.height/2

    l = rectangle.length

    # defining lowest line
    lowest_line = np.linspace(left_end, right_end, int(n_per_length*l))
    lower_ys = np.ones_like(lowest_line)*lower_end
    upper_ys = np.ones_like(lowest_line)*upper_end

    # defining pdp (potential dot positions)
    pdp = np.zeros([len(lowest_line), 2])
    pdp[:, 0] = lowest_line
    pdp[:, 1] = lower_ys

    # adding other boundary points as pdps
    # upper line
    pdp_append = np.zeros([len(lowest_line), 2])
    pdp_append[:, 0] =  lowest_line
    pdp_append[:, 1] = upper_ys

    left_line = np.linspace(lower_end, upper_end, int(n_per_length*l))
    left_xs = np.ones_like(left_line)*left_end
    right_xs = np.ones_like(left_line)*right_end

    pdp_left = np.zeros([len(left_line), 2])
    pdp_left[:, 1] = left_line
    pdp_left[:, 0] = left_xs

    pdp_right = np.zeros([len(left_line), 2])
    pdp_right[:, 1] = left_line
    pdp_right[:, 0] = right_xs

    pdp = np.append(pdp, pdp_append, axis=0)
    pdp = np.append(pdp, pdp_left, axis=0)
    pdp = np.append(pdp, pdp_right, axis=0)

    # defining dot locations:
    dot_locations = np.zeros([0, 2])
    points = 0

    # algorithm
    while(True):
        lowest_pdp_i = np.argmin(pdp[:, 1])
        lowest_pdp = pdp[lowest_pdp_i]
        dot_locations = np.append(dot_locations, [lowest_pdp], axis=0)
        points = points + 1

        if lowest_pdp[1] > upper_end : break
        if points == 9000 : break

        r0 = density_function(lowest_pdp[0], lowest_pdp[1])
        points_in_circle = np.zeros(0, dtype=int)
        for i in range(len(pdp)):
            if radius(pdp[i], lowest_pdp) < r0 :
                points_in_circle = np.append(points_in_circle, [i], axis=0)
            elif pdp[i, 0] < left_end or pdp[i, 0] > right_end or pdp[i, 1] < lower_end: # or pdp[i, 1] > upper_end:
                points_in_circle = np.append(points_in_circle, [i], axis=0)
        pdp = np.delete(pdp, points_in_circle, 0)

        closest_point_left = np.array([left_end, upper_end])
        closest_point_right = np.array([right_end, upper_end])

        for i in range(len(pdp)):
            r1 = radius(closest_point_left, lowest_pdp)
            r2 = radius(closest_point_right, lowest_pdp)
            ri = radius(pdp[i], lowest_pdp)

            if pdp[i, 0] - lowest_pdp[0] < 0 and ri < r1:
                closest_point_left = pdp[i]
            elif pdp[i, 0] - lowest_pdp[0] > 0 and ri < r2:
                closest_point_right = pdp[i]

        # we have closest points 1 and 2 now
        if closest_point_left[0] - lowest_pdp[0] == 0: theta1 = np.pi/2
        else: 
            if closest_point_left[0] - lowest_pdp[0] < 0:
                if closest_point_left[1] - lowest_pdp[1] < 0:
                    theta1 = np.pi + np.arctan((closest_point_left[1] - lowest_pdp[1])/(closest_point_left[0] - lowest_pdp[0]))
                else:
                    theta1 = np.pi + np.arctan((closest_point_left[1] - lowest_pdp[1])/(closest_point_left[0] - lowest_pdp[0]))
            else:
                if closest_point_left[1] - lowest_pdp[1] < 0:
                    theta1 = np.arctan((closest_point_left[1] - lowest_pdp[1])/(closest_point_left[0] - lowest_pdp[0]))
                else:
                    theta1 = np.arctan((closest_point_left[1] - lowest_pdp[1])/(closest_point_left[0] - lowest_pdp[0]))

        if closest_point_right[0] - lowest_pdp[0] == 0: theta2 = np.pi/2
        else: 
            if closest_point_right[0] - lowest_pdp[0] < 0:
                if closest_point_right[1] - lowest_pdp[1] < 0:
                    theta2 = np.pi + np.arctan((closest_point_right[1] - lowest_pdp[1])/(closest_point_right[0] - lowest_pdp[0]))
                else:
                    theta2 = np.pi + np.arctan((closest_point_right[1] - lowest_pdp[1])/(closest_point_right[0] - lowest_pdp[0]))
            else:
                if closest_point_right[1] - lowest_pdp[1] < 0:
                    theta2 = np.arctan((closest_point_right[1] - lowest_pdp[1])/(closest_point_right[0] - lowest_pdp[0]))
                else:
                    theta2 = np.arctan((closest_point_right[1] - lowest_pdp[1])/(closest_point_right[0] - lowest_pdp[0]))

        theta_news = np.linspace(theta1, theta2, 5)

        x_news = lowest_pdp[0] + r0*np.cos(theta_news)
        y_news = lowest_pdp[1] + r0*np.sin(theta_news)

        pdp_news = np.zeros([5, 2])

        pdp_news[:, 0] = x_news
        pdp_news[:, 1] = y_news

        pdp = np.append(pdp, pdp_news, 0)

    print(len(dot_locations))
    plt.scatter(dot_locations[:, 0], dot_locations[:, 1], s = 0.6)
    plt.show()

    mean, std_var, diff = stats(dot_locations)
    print(mean)
    print(std_var)
    print(diff)

def radius(point1, point2):
    x1 = point1[0]
    x2 = point2[0]
    y1 = point1[1]
    y2 = point2[1]

    return np.sqrt((x1-x2)**2 + (y1-y2)**2)

def uniform(x, y):
    return 0.05


if __name__ == "__main__":
    box1 = Rectangle(5, 3, 2.5, 1.5)
    fornberg_distribution(uniform, box1, 100)

