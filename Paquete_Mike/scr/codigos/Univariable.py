import math

# Método de división de intervalos por la mitad
def bisection_method(f, a, b, tol=1e-6, max_iter=1000):
    if f(a) * f(b) >= 0:
        raise ValueError("La función no cambia de signo en el intervalo dado [a, b].")

    # Inicialización de los extremos del intervalo
    left = a
    right = b

    # Iteración hasta alcanzar la convergencia o el número máximo de iteraciones
    for i in range(max_iter):
        # Punto medio del intervalo
        midpoint = (left + right) / 2.0

        # Valor de la función en el punto medio
        f_mid = f(midpoint)

        # Verifica si se alcanza la tolerancia
        if abs(f_mid) < tol:
            print(f'Convergencia alcanzada en {i+1} iteraciones')
            return midpoint

        # Decide en qué mitad del intervalo continuar
        if f(left) * f_mid < 0:
            right = midpoint  # La raíz está en la mitad izquierda
        else:
            left = midpoint  # La raíz está en la mitad derecha

    raise ValueError(f'El método de bisección no converge después de {max_iter} iteraciones.')

# Búsqueda de Fibonacci
def fibonacci_search(f, a, b, n, tol=1e-6):
    # Números de Fibonacci hasta n
    fib = [0, 1]
    for i in range(2, n+1):
        fib.append(fib[i-1] + fib[i-2])
    
    # Longitud inicial del intervalo
    L = b - a
    
    for k in range(1, n):
        # Puntos intermedios basados en los números de Fibonacci
        x1 = a + (fib[n-k-1] / fib[n-k+1]) * L
        x2 = a + (fib[n-k] / fib[n-k+1]) * L
        
        # Evaluar la función en los puntos intermedios
        fx1 = f(x1)
        fx2 = f(x2)
        
        # Reducir el intervalo
        if fx1 < fx2:
            b = x2
        else:
            a = x1
        
        L = b - a
        
        # Condición de convergencia
        if L < tol:
            return (a + b) / 2
    
    return (a + b) / 2

# Método de la sección dorada
def golden_section_search(f, a, b, tol=1e-6):
    # Proporción áurea
    phi = (math.sqrt(5) - 1) / 2
    
    # Puntos iniciales
    x1 = a + (1 - phi) * (b - a)
    x2 = a + phi * (b - a)
    
    # Valores de la función en los puntos iniciales
    fx1 = f(x1)
    fx2 = f(x2)
    
    # Iteración hasta convergencia
    while abs(b - a) > tol:
        if fx1 < fx2:
            b = x2
            x2 = x1
            fx2 = fx1
            x1 = a + (1 - phi) * (b - a)
            fx1 = f(x1)
        else:
            a = x1
            x1 = x2
            fx1 = fx2
            x2 = a + phi * (b - a)
            fx2 = f(x2)
    
    # Devolver el punto medio del intervalo final como aproximación del mínimo
    return (a + b) / 2

# Metodo de Newton-Raphson
def newton_raphson(f, df, x0, tol=1e-6, max_iter=1000):
    x = x0
    for i in range(max_iter):
        fx = f(x)
        if abs(fx) < tol:
            print(f'Convergencia alcanzada en {i+1} iteraciones')
            return x
        dfx = df(x)
        if dfx == 0:
            raise ValueError("Derivada de la función es cero. No se puede continuar.")
        x = x - fx / dfx
    raise ValueError("El método de Newton-Raphson no converge después de {max_iter} iteraciones.")

# Metodo de Biseccion
def bisection(f, a, b, tol=1e-6, max_iter=1000):
    if f(a) * f(b) >= 0:
        raise ValueError("La función no cambia de signo en el intervalo dado.")
    
    for i in range(max_iter):
        c = (a + b) / 2.0
        fc = f(c)
        
        if abs(fc) < tol:
            print(f'Convergencia alcanzada en {i+1} iteraciones')
            return c
        
        if f(a) * fc < 0:
            b = c
        else:
            a = c
    
    raise ValueError("El método de bisección no converge después de {max_iter} iteraciones.")

# Metodo de la Secante
def secant_method(f, x0, x1, tol=1e-6, max_iter=1000):
    fx0 = f(x0)
    fx1 = f(x1)
    
    for i in range(max_iter):
        if abs(fx1) < tol:
            print(f'Convergencia alcanzada en {i+1} iteraciones')
            return x1
        
        # Calcula la siguiente aproximación usando el método de la secante
        x_next = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        
        x0 = x1
        x1 = x_next
        fx0 = fx1
        fx1 = f(x_next)
    
    raise ValueError("El método de la secante no converge después de {max_iter} iteraciones.")
