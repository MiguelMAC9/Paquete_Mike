import random
import math

# Random Walk
def random_walk_2d(n_steps, initial_position, mu=0, sigma=1):
    x, y = [initial_position[0]], [initial_position[1]]  # Posición inicial
    
    for _ in range(n_steps):
        step_x = random.gauss(mu, sigma)  # Incremento aleatorio en x
        step_y = random.gauss(mu, sigma)  # Incremento aleatorio en y
        x.append(x[-1] + step_x)
        y.append(y[-1] + step_y)
    return x, y

# Nelder-Mead Simplex
def NelderMead_Simplex(funcion, x0, tol=1e-6, iteraciones=100, alpha=1, beta=0.5, gamma=2):
    n = len(x0)
    vectores = [x0[:]]  # Crear una lista de listas
    for i in range(n):
        nuevo_vector = x0[:]
        nuevo_vector[i] += 0.05
        vectores.append(nuevo_vector)
    
    valores_vect = [funcion(v) for v in vectores]

    for _ in range(iteraciones):
        ordenar_vect = sorted(range(len(valores_vect)), key=lambda k: valores_vect[k])
        vectores = [vectores[i] for i in ordenar_vect]
        valores_vect = [valores_vect[i] for i in ordenar_vect]

        xc = [sum(vectores[j][i] for j in range(n)) / n for i in range(n)]

        xr = [xc[i] + alpha * (xc[i] - vectores[-1][i]) for i in range(n)]
        fr = funcion(xr)
        if valores_vect[-2] > fr >= valores_vect[0]:
            vectores[-1], valores_vect[-1] = xr, fr
        elif fr < valores_vect[0]:
            expansion = [xc[i] + gamma * (xr[i] - xc[i]) for i in range(n)]
            fe = funcion(expansion)
            if fe < fr:
                vectores[-1], valores_vect[-1] = expansion, fe
            else:
                vectores[-1], valores_vect[-1] = xr, fr
        else:
            contraccion = [xc[i] + beta * (vectores[-1][i] - xc[i]) for i in range(n)]
            fc = funcion(contraccion)
            if fc < valores_vect[-1]:
                vectores[-1], valores_vect[-1] = contraccion, fc
            else:
                for i in range(1, len(vectores)):
                    for j in range(n):
                        vectores[i][j] = (vectores[0][j] + vectores[i][j]) / 2
                    valores_vect[i] = funcion(vectores[i])

    return vectores[0]

# Hooke-Jeeves
def Busqueda(x, d, funcion, limite=1e10):
    x_i = x[:]
    for i in range(len(x)):
        for direction in [-1, 1]:
            x_t = x_i[:]
            x_t[i] += direction * d
            if abs(x_t[i]) > limite:
                continue
            if funcion(x_t) < funcion(x_i):
                x_i = x_t
    return x_i

def hooke_jeeves(x_i, delta, alpha, e, n_iter, funcion, limite=1e10):
    x_b = x_i[:]
    x_m = x_b[:]
    iter_c = 0
    resul = [x_b[:]]

    while delta > e and iter_c < n_iter:
        x_n = Busqueda(x_b, delta, funcion, limite)
        if funcion(x_n) < funcion(x_m):
            x_b = [2 * x_n[i] - x_m[i] for i in range(len(x_n))]
            x_m = x_n[:]
        else:
            delta *= alpha
            x_b = x_m[:]
        resul.append(x_b[:])
        iter_c += 1

    return x_m, resul

# Cauchy
def gradiente(f, x, deltaX=0.001):
    grad = []
    for i in range(len(x)):
        xp = x[:]
        xn = x[:]
        xp[i] += deltaX
        xn[i] -= deltaX
        grad.append((f(xp) - f(xn)) / (2 * deltaX))
    return grad

def cauchy(funcion, x0, epsilon1, epsilon2, M, optimizador_univariable):
    terminar = False
    xk = x0[:]
    k = 0
    while not terminar:
        grad = gradiente(funcion, xk)

        if math.sqrt(sum(g**2 for g in grad)) < epsilon1 or k >= M:
            terminar = True
        else:
            def alpha_funcion(alpha):
                return funcion([xk[i] - alpha * grad[i] for i in range(len(xk))])

            alpha = optimizador_univariable(alpha_funcion, epsilon2, a=0.0, b=1.0)
            x_k1 = [xk[i] - alpha * grad[i] for i in range(len(xk))]
            print(xk, alpha, grad, x_k1)

            if math.sqrt(sum((x_k1[i] - xk[i])**2 for i in range(len(xk)))) / (math.sqrt(sum(xi**2 for xi in xk)) + 0.00001) <= epsilon2:
                terminar = True
            else:
                k += 1
                xk = x_k1
    return xk

# Método de Fletcher-Reeves
def gradient(f, x, deltaX=1e-5):
    grad = [0] * len(x)
    for i in range(len(x)):
        x1 = x[:]
        x2 = x[:]
        x1[i] += deltaX
        x2[i] -= deltaX
        grad[i] = (f(x1) - f(x2)) / (2 * deltaX)
    return grad

def fletcher_reeves(f, x0, tol=1e-6, max_iter=1000):
    x = x0[:]
    grad_f = gradient(f, x)
    r = [-g for g in grad_f]
    d = r[:]
    rsold = sum(ri**2 for ri in r)
    
    for i in range(max_iter):
        alpha = line_search(f, x, d)
        x = [x[j] + alpha * d[j] for j in range(len(x))]
        grad_f = gradient(f, x)
        r = [-g for g in grad_f]
        rsnew = sum(ri**2 for ri in r)
        
        if math.sqrt(rsnew) < tol:
            print(f'Convergencia alcanzada en {i+1} iteraciones')
            break
        
        beta = rsnew / rsold
        d = [r[j] + beta * d[j] for j in range(len(r))]
        rsold = rsnew
        
    return x

def line_search(f, x, d, alpha0=1, c=1e-4, tau=0.9):
    alpha = alpha0
    while f([x[j] + alpha * d[j] for j in range(len(x))]) > f(x) + c * alpha * sum(gradient(f, x)[j] * d[j] for j in range(len(x))):
        alpha *= tau
    return alpha

# Método de Newton
def hessian_matrix(f, x, deltaX):
    fx = f(x)
    N = len(x)
    H = []
    for i in range(N):
        hi = []
        for j in range(N):
            if i == j:
                xp = x[:]
                xn = x[:]
                xp[i] += deltaX
                xn[i] -= deltaX
                hi.append((f(xp) - 2 * fx + f(xn)) / (deltaX ** 2))
            else:
                xpp = x[:]
                xpn = x[:]
                xnp = x[:]
                xnn = x[:]
                xpp[i] += deltaX
                xpp[j] += deltaX
                xpn[i] += deltaX
                xpn[j] -= deltaX
                xnp[i] -= deltaX
                xnp[j] += deltaX
                xnn[i] -= deltaX
                xnn[j] -= deltaX
                hi.append((f(xpp) - f(xpn) - f(xnp) + f(xnn)) / (4 * deltaX ** 2))
        H.append(hi)
    return H

def newton_method(f, x0, epsilon=1e-6, max_iter=100):
    x = x0[:]
    for k in range(max_iter):
        grad = gradient(f, x)
        H = hessian_matrix(f, x, 1e-5)
        H_inv = invert_matrix(H)
        d = [-sum(H_inv[i][j] * grad[j] for j in range(len(grad))) for i in range(len(H))]
        alpha = 1.0  
        while f([x[j] + alpha * d[j] for j in range(len(x))]) > f(x): 
            alpha *= 0.5
        x = [x[j] + alpha * d[j] for j in range(len(x))]
        if math.sqrt(sum(g**2 for g in grad)) < epsilon:
            break
    return x

def invert_matrix(matrix):
    n = len(matrix)
    identity = [[float(i == j) for i in range(n)] for j in range(n)]
    for i in range(n):
        factor = matrix[i][i]
        for j in range(n):
            matrix[i][j] /= factor
            identity[i][j] /= factor
        for k in range(n):
            if k != i:
                factor = matrix[k][i]
                for j in range(n):
                    matrix[k][j] -= factor * matrix[i][j]
                    identity[k][j] -= factor * identity[i][j]
    return identity
