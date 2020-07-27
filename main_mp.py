from operator import inv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize, basinhopping
import mpmath

R = 1
C = np.array([0, -2])


acc = []


def fullreflection(d, s):
    bounces = []

    v = s - C

    e = np.square(v.dot(d)) - (v.dot(v) - R * R)

    if e < 0:
        return None, None, []

    # p = -v.dot(d) + np.sqrt(e)
    # m = -v.dot(d) - np.sqrt(e)

    # if np.linalg.norm(p - d) < np.linalg.norm(m - d):
    #     print("pos")
    #     t = -v.dot(d) + np.sqrt(e)
    # else:
    #     t = -v.dot(d) - np.sqrt(e)
    #     print("neg")

    t = -v.dot(d) - np.sqrt(e)
    intersection_point = s + t * d

    bounces.append(intersection_point)

    n = (intersection_point - C) / np.linalg.norm(intersection_point - C)

    reflection_dir = d - 2 * (n.dot(d)) * n

    # print("ref", d, n, intersection_point, reflection_dir)
    if reflection_dir[1] < 0:
        return None, None, bounces

    plane_n = np.array([0, -1])
    mu = (np.array([0, 0]) - intersection_point).dot(plane_n) / (reflection_dir.dot(plane_n))

    plane_intersection_point = intersection_point + mu * reflection_dir
    plane_reflection_dir = -(
        np.array([1, -1])
        * (intersection_point - plane_intersection_point)
        / np.linalg.norm(intersection_point - plane_intersection_point)
    )
    bounces.append(plane_intersection_point)
    return plane_intersection_point, plane_reflection_dir, bounces


def num_reflections(x):
    d = np.array([x + 2, -2]) / np.linalg.norm(np.array([x + 2, 2]))

    s = np.array([x, 0])

    bounces = []

    while d is not None:
        s, d, b = fullreflection(d, s)
        # if d is not None:
        # if np.isnan(d[0]):
        #     break

        bounces.extend(b)

    acc.append([x, len(bounces)])
    return np.array([[-2, -2]] + [[x, 0]] + bounces)


def shittyoptimize(min_x, max_x, n):
    step = (max_x - min_x) / n
    x = min_x
    result = []
    for i in range(n):
        x = i * step + min_x
        result.append([x, num_reflections(mpmath.mpf(x)).shape[0]])

    max_bounces = max(i[1] for i in result)

    min_idx = min(i for i, res in enumerate(result) if res[1] == max_bounces)
    max_idx = max(i for i, res in enumerate(result) if res[1] == max_bounces)

    print(max_bounces, min_idx, max_idx, min_idx * step + min_x)
    if min_idx == max_idx:
        shittyoptimize(min_x + step * (min_idx - 1), min_x + step * (max_idx + 1), n)
    else:
        shittyoptimize(min_x + step * (min_idx - 1), min_x + step * (max_idx + 1), n * 2)


def main():
    mpmath.mp.dps = 100
    # print(fullreflection(np.array([mpmath.mpf(1), mpmath.mpf(-1)]), np.array([mpmath.mpf(-1), 0])))
    # print(num_reflections(mpmath.mpf(-0.822486324943389801589432863693)))

    # invbounces = lambda x: 1 / max(len(num_reflections(x[0])), 1)
    # a = mpmath.calculus.optimization.Secant(
    #     ctx=1, f=invbounces, x0=[mpmath.mpf(-0.822486324943389801589432863693)], verbose=True
    # )
    # print(num_reflections(mpmath.mpf(-0.822486324943389801589432863693)).shape)

    min_x = mpmath.mpf("-0.8224863249433905787455501013027969747781753540039062500000000000")
    max_x = mpmath.mpf("-0.8224863249433896905671304011775646358728408813476562500000000000")
    mpmath.nprint(min_x, 100)
    shittyoptimize(min_x, max_x, 100)
    # circle = np.linspace(-R / 2, R / 2, 20)
    # xs = np.linspace(min_x, max_x, 1000000)

    # # for x in xs:
    # # num_reflections(x)

    # # x = -1.0

    # # x0 = -0.822486324943389801589432863693

    # # f = minimize_scalar(invbounces, bracket=[-0.9, -0.8, -0.7], options={"disp": True, "maxiter": 500000}, tol=1e-40)
    # f = minimize(
    #     invbounces, -0.8224863230969456, method="Nelder-Mead", options={"maxiter": 10000, "xtol": 1e-60, "ftol": 1e-60}
    # )

    # a = np.array(acc)
    # # print()
    # # a = a[a[:, 1] == 44, :]
    # # print(format(a, ".60g"))
    # print(a)
    # # print(format(np.max(a[a[:, 1] > 43, 0]), ".60g"))


if __name__ == "__main__":
    main()
