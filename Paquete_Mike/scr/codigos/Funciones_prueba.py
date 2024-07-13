import math

def himmelblau(x):
    return (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2

def testfunction(x):
    return x[0]**2 + x[1]**2

def sphere(x):
    return sum(xi**2 for xi in x)

def rastrigin(x, A=10):
    return A * len(x) + sum(xi**2 - A * math.cos(2 * math.pi * xi) for xi in x)

def rosenbrock(x):
    return sum(100 * (x[i+1] - x[i]**2)**2 + (x[i] - 1)**2 for i in range(len(x) - 1))

def beale(x):
    return ((1.5 - x[0] + x[0] * x[1])**2 +
            (2.25 - x[0] + x[0] * x[1]**2)**2 +
            (2.625 - x[0] + x[0] * x[1]**3)**2)

def goldstein(x):
    part1 = (1 + (x[0] + x[1] + 1)**2 * 
            (19 - 14 * x[0] + 3 * x[0]**2 - 14 * x[1] + 6 * x[0] * x[1] + 3 * x[1]**2))
    part2 = (30 + (2 * x[0] - 3 * x[1])**2 * 
            (18 - 32 * x[0] + 12 * x[0]**2 + 48 * x[1] - 36 * x[0] * x[1] + 27 * x[1]**2))
    return part1 * part2

def boothfunction(x):
    return (x[0] + 2 * x[1] - 7)**2 + (2 * x[0] + x[1] - 5)**2

def bunkinn6(x):
    return 100 * math.sqrt(abs(x[1] - 0.001 * x[0]**2)) + 0.01 * abs(x[0] + 10)

def matyas(x):
    return 0.26 * (x[0]**2 + x[1]**2) - 0.48 * x[0] * x[1]

def levi(x):
    part1 = math.sin(3 * math.pi * x[0])**2
    part2 = (x[0] - 1)**2 * (1 + math.sin(3 * math.pi * x[1])**2)
    part3 = (x[1] - 1)**2 * (1 + math.sin(2 * math.pi * x[1])**2)
    return part1 + part2 + part3

def threehumpcamel(x):
    return 2 * x[0]**2 - 1.05 * x[0]**4 + (x[0]**6) / 6 + x[0] * x[1] + x[1]**2

def easom(x):
    return -math.cos(x[0]) * math.cos(x[1]) * math.exp(-(x[0] - math.pi)**2 - (x[1] - math.pi)**2)

def crossintray(x):
    op = abs(math.sin(x[0]) * math.sin(x[1]) * math.exp(abs(100 - math.sqrt(x[0]**2 + x[1]**2) / math.pi)))
    return -0.0001 * (op + 1)**0.1

def eggholder(x):
    op1 = -(x[1] + 47) * math.sin(math.sqrt(abs(x[0] / 2 + (x[1] + 47))))
    op2 = -x[0] * math.sin(math.sqrt(abs(x[0] - (x[1] + 47))))
    return op1 + op2

def holdertable(x):
    op = abs(math.sin(x[0]) * math.cos(x[1]) * math.exp(abs(1 - math.sqrt(x[0]**2 + x[1]**2) / math.pi)))
    return -op

def mccormick(x):
    return math.sin(x[0] + x[1]) + (x[0] - x[1])**2 - 1.5 * x[0] + 2.5 * x[1] + 1

def schaffern2(x):
    numerator = math.sin(x[0]**2 - x[1]**2)**2 - 0.5
    denominator = (1 + 0.001 * (x[0]**2 + x[1]**2))**2
    return 0.5 + numerator / denominator

def schaffern4(x):
    num = math.cos(math.sin(abs(x[0]**2 - x[1]**2))) - 0.5
    den = (1 + 0.001 * (x[0]**2 + x[1]**2))**2
    return 0.5 + num / den

def styblinskitang(x):
    return sum((xi**4 - 16 * xi**2 + 5 * xi) / 2 for xi in x)

def shekel(x, a=None, c=None):
    if a is None:
        a = [
            [4.0, 4.0, 4.0, 4.0],
            [1.0, 1.0, 1.0, 1.0],
            [8.0, 8.0, 8.0, 8.0],
            [6.0, 6.0, 6.0, 6.0],
            [3.0, 7.0, 3.0, 7.0],
            [2.0, 9.0, 2.0, 9.0],
            [5.0, 5.0, 3.0, 3.0],
            [8.0, 1.0, 8.0, 1.0],
            [6.0, 2.0, 6.0, 2.0],
            [7.0, 3.6, 7.0, 3.6]
        ]
    if c is None:
        c = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5]
        
    m = len(c)
    s = 0
    for i in range(m):
        s -= 1 / (sum((x[j] - a[i][j])**2 for j in range(2)) + c[i])
    return s
