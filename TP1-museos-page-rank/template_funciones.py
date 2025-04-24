import numpy as np

def construye_adyacencia(D,m): 
    # Función que construye la matriz de adyacencia del grafo de museos
    # D matriz de distancias, m cantidad de links por nodo
    # Retorna la matriz de adyacencia como un numpy.
    D = D.copy()
    l = [] # Lista para guardar las filas
    for fila in D: # recorriendo las filas, anexamos vectores lógicos
        l.append(fila<=fila[np.argsort(fila)[m]] ) # En realidad, elegimos todos los nodos que estén a una distancia menor o igual a la del m-esimo más cercano
    A = np.asarray(l).astype(int) # Convertimos a entero
    np.fill_diagonal(A,0) # Borramos diagonal para eliminar autolinks
    return(A)

def calculaLU(matriz):
    # matriz es una matriz de NxN
    # Retorna la factorización LU a través de una lista con dos matrices L y U de NxN.
    # Completar! Have fun
    m=matriz.shape[0]
    n=matriz.shape[1]
    
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    L = np.eye(matriz.shape[0])
    U = matriz.copy()
    for j in range(m):
        for i in range(j+1, n):
            L[i,j]= U[i,j]/U[j,j] 
            U[i,:] = U[i,:] - L[i,j]*U[j,:]   
    
    return L, U

def resolucion_sistema(L, U, b):
    m = L.shape[0]

    # Sustitución hacia adelante: Ly = b
    y = np.zeros_like(b, dtype=float)
    for i in range(m):
        y[i] = b[i] - np.dot(L[i, :i], y[:i])

    # Sustitución hacia atrás: Ux = y
    x = np.zeros_like(b, dtype=float)
    for i in range(m-1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i+1:], x[i+1:])) / U[i, i]

    return x


def solucion_final(A, b):
    L, U = calculaLU(A)
    x_res = resolucion_sistema(L, U, b)
    return x_res


def calcula_matriz_C(A): 
    # Función para calcular la matriz de trancisiones C
    # A: Matriz de adyacencia
    # Retorna la matriz C
       
    AT = A.T

    grados = np.sum(AT, axis=1)

    Kinv = np.diag([1/g if g != 0 else 0 for g in grados]) # Calcula inversa de la matriz K, que tiene en su diagonal la suma por filas de A
    C = AT@Kinv # Calcula C multiplicando Kinv y A
    return C

    
def calcula_pagerank(A,alfa):
    # Función para calcular PageRank usando LU
    # A: Matriz de adyacencia
    # d: coeficientes de damping
    # Retorna: Un vector p con los coeficientes de page rank de cada museo
    C = calcula_matriz_C(A)
    N = A.shape[0] # Obtenemos el número de museos N a partir de la estructura de la matriz A
    M = (N/alfa) * ( np.eye(N) - ((1 - alfa) * C))
    L, U = calculaLU(M) # Calculamos descomposición LU a partir de C y d

    b =  np.ones(N) # Vector de 1s, multiplicado por el coeficiente correspondiente usando d y N. CHECKEAR!

    Up = scipy.linalg.solve_triangular(L,b,lower=True) # Primera inversión usando L
    p = scipy.linalg.solve_triangular(U,Up) # Segunda inversión usando U
    return p

def calcula_matriz_C_continua(D): 
    # Función para calcular la matriz de trancisiones C
    # A: Matriz de adyacencia
    # Retorna la matriz C en versión continua
    D = D.copy()
    F = 1/D
    np.fill_diagonal(F,0)

    suma_filas = np.sum(F, axis=1)

    Kinv = np.diag([1/s if s != 0 else 0 for s in suma_filas]) # Calcula inversa de la matriz K, que tiene en su diagonal la suma por filas de F 
    C = Kinv @ F # Calcula C multiplicando Kinv y F
    return C

def calcula_B(C,cantidad_de_visitas):
    # Recibe la matriz T de transiciones, y calcula la matriz B que representa la relación entre el total de visitas y el número inicial de visitantes
    # suponiendo que cada visitante realizó cantidad_de_visitas pasos
    # C: Matirz de transiciones
    # cantidad_de_visitas: Cantidad de pasos en la red dado por los visitantes. Indicado como r en el enunciado
    # Retorna:Una matriz B que vincula la cantidad de visitas w con la cantidad de primeras visitas v
    B = np.eye(C.shape[0])
    for i in range(cantidad_de_visitas-1):
        # Sumamos las matrices de transición para cada cantidad de pasos
    return B