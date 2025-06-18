import numpy as np
import funciones.template_funciones as f1
# Matriz A de ejemplo
A_ejemplo = np.array([
   [0, 1, 1, 1, 0, 0, 0, 0],
   [1, 0, 1, 1, 0, 0, 0, 0],
   [1, 1, 0, 1, 0, 1, 0, 0],
   [1, 1, 1, 0, 1, 0, 0, 0],
   [0, 0, 0, 1, 0, 1, 1, 1],
   [0, 0, 1, 0, 1, 0, 1, 1],
   [0, 0, 0, 0, 1, 1, 0, 1],
   [0, 0, 0, 0, 1, 1, 1, 0]
])

def calcula_L(A):
    # La función recibe la matriz de adyacencia A y calcula la matriz laplaciana
    K = np.diag(np.sum(A_ejemplo,axis=0))
    L = K - A_ejemplo
    return L

def calcula_P(A):
    k_array = np.sum(A,axis=0)
    E = np.sum(k_array) / 2

    P = np.zeros((A.shape[0], A.shape[1]))
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            P[i, j] = (k_array[i] * k_array[j]) / (2 * E)
    
    return P

def calcula_R(A):
    # La funcion recibe la matriz de adyacencia A y calcula la matriz de modularidad
    P = calcula_P(A)
    R = A - P
    return R

# AUX
def calcula_s(v):
    # Usamos sugerencia del enunciado: s = signo de cada elemento de v
    s = np.sign(v) # signos de v, excepto ceros
    s[s == 0] = 1 # tomamos al 0 como signo positivo
    
    return s

def calcula_lambda(L,v):
    # Recibe L y v y retorna el corte asociado
    
    # Primero de v armo s y s traspuesta
    s = calcula_s(v)
    s_t = s.T
    
    lambda_corte = ( 1 / 4 ) * ( s_t @ L @ s )
    
    return lambda_corte

def calcula_Q(R,v):
    # La funcion recibe R y s y retorna la modularidad (a menos de un factor 2E)
    
    # Primero con "v" armamos "s" y "s" traspuesta
    s = calcula_s(v)
    s_t = s.T
    
    # R = A - P
    Q = ( s_t @ R @ s ) # * (1 / (4 * E ) ), que lo omitimos por sugerencia del enunciado
    
    return Q

# NOTA: Probé esta matriz

# A = np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
# A tiene AVAS: {3, 1, 0}
# asociados a AVES respectivos: # (1, 1, 1)
                                # (-1, 0, 1)
                                # (1, -2, 1)
                                
# metpot de A debe ser: (1, 1, 1) o normalizado?, 3, true

def metpot1(A,tol=1e-16,maxrep=np.Inf):
    # Recibe una matriz A y calcula su autovalor de mayor módulo, con un error relativo menor a tol y-o haciendo como mucho maxrep repeticiones
    rng = np.random.default_rng(1148)
    v = rng.random(A.shape[1]) # Generamos un vector de partida aleatorio, entre -1 y 1
    v = v / f1.norma_2_matriz(v) # Lo normalizamos

    v1 = A @ v # Aplicamos la matriz una vez
    v1 = v1 / f1.norma_2_matriz(v1) # normalizamos

    l = np.dot(v, A @ v) # Calculamos el autovector estimado
    l1 = np.dot(v1, A @ v1) # Y el estimado en el siguiente paso

    nrep = 0 # Contador
    while np.abs(l1-l)/np.abs(l) > tol and nrep < maxrep: # Si estamos por debajo de la tolerancia buscada 
        v = v1 # actualizamos v y repetimos
        l = l1
        v1 = A @ v # Calculo nuevo v1
        v1 = v1 / f1.norma_2_matriz(v1) # Normalizo
        l1 = np.dot(v1, A @ v1) # Calculo autovector
        nrep += 1 # Un pasito mas
    if not nrep < maxrep:
        print('MaxRep alcanzado')
    l = np.dot(v1, A @ v1) # Calculamos el autovalor
    return v1,l,nrep<maxrep

# np.outer( v1, v2 ) multiplica a v2 por cada elemento de v1 y lo devuelve en listas diferentes
def deflaciona(A,tol=1e-16):
    # M - lambda * (v * vT)/(vT * v)
    # Recibe la matriz A, una tolerancia para el método de la potencia, y un número máximo de repeticiones
    v1,l,_ = metpot1(A,tol) # Buscamos primer autovalor (dominante) con el método de la potencia
    
    deflA = A - l * np.outer(v1, v1) # Sugerencia, usar la funcion outer de numpy
    return deflA

def metpotI(A,mu,tol=1e-16,maxrep=np.Inf):
    # Retorna el primer autovalor de la inversa de A + mu * I, junto a su autovector y si el método convergió.
    A_shifted = A + mu * np.eye(A.shape[0])
    A_shifted_inv = f1.inversa_por_LU(A_shifted)
    return metpot1(A_shifted_inv,tol=tol,maxrep=maxrep)

def metpotI2(A,mu,tol=1e-16,maxrep=np.Inf):
   X = A + mu * np.eye(A.shape[0]) # Calculamos la matriz A shifteada en mu
   iX = f1.inversa_por_LU(X) # La invertimos
   defliX = deflaciona(iX, tol=tol) # La deflacionamos
   v,l,_ =  metpot1(defliX, tol=tol, maxrep=maxrep) # Buscamos su segundo autovector
   print("l: ", l)
   l = 1/l # Reobtenemos el autovalor correcto
   l -= mu
   return v,l


def laplaciano_iterativo(A,niveles,nombres_s=None):
    # Recibe una matriz A, una cantidad de niveles sobre los que hacer cortes, y los nombres de los nodos
    # Retorna una lista con conjuntos de nodos representando las comunidades.
    # La función debe, recursivamente, ir realizando cortes y reduciendo en 1 el número de niveles hasta llegar a 0 y retornar.
    if nombres_s is None: # Si no se proveyeron nombres, los asignamos poniendo del 0 al N-1
        nombres_s = range(A.shape[0])
    if A.shape[0] == 1 or niveles == 0: # Si llegamos al último paso, retornamos los nombres en una lista
        return([nombres_s])
    else: # Sino:
        L = calcula_L(A) # Recalculamos el L
        v,l,_ = ... # Encontramos el segundo autovector de L
        # Recortamos A en dos partes, la que está asociada a el signo positivo de v y la que está asociada al negativo
        Ap = ... # Asociado al signo positivo
        Am = ... # Asociado al signo negativo
        
        return(
                laplaciano_iterativo(Ap,niveles-1,
                                     nombres_s=[ni for ni,vi in zip(nombres_s,v) if vi>0]) +
                laplaciano_iterativo(Am,niveles-1,
                                     nombres_s=[ni for ni,vi in zip(nombres_s,v) if vi<0])
                )        


def modularidad_iterativo(A=None,R=None,nombres_s=None):
    # Recibe una matriz A, una matriz R de modularidad, y los nombres de los nodos
    # Retorna una lista con conjuntos de nodos representando las comunidades.

    if A is None and R is None:
        print('Dame una matriz')
        return(np.nan)
    if R is None:
        R = calcula_R(A)
    if nombres_s is None:
        nombres_s = range(R.shape[0])
    # Acá empieza lo bueno
    if R.shape[0] == 1: # Si llegamos al último nivel
        return(...)
    else:
        v,l,_ = ... # Primer autovector y autovalor de R
        # Modularidad Actual:
        Q0 = np.sum(R[v>0,:][:,v>0]) + np.sum(R[v<0,:][:,v<0])
        if Q0<=0 or all(v>0) or all(v<0): # Si la modularidad actual es menor a cero, o no se propone una partición, terminamos
            return(...)
        else:
            ## Hacemos como con L, pero usando directamente R para poder mantener siempre la misma matriz de modularidad
            Rp = ... # Parte de R asociada a los valores positivos de v
            Rm = ... # Parte asociada a los valores negativos de v
            vp,lp,_ = ...  # autovector principal de Rp
            vm,lm,_ = ... # autovector principal de Rm
        
            # Calculamos el cambio en Q que se produciría al hacer esta partición
            Q1 = 0
            if not all(vp>0) or all(vp<0):
               Q1 = np.sum(Rp[vp>0,:][:,vp>0]) + np.sum(Rp[vp<0,:][:,vp<0])
            if not all(vm>0) or all(vm<0):
                Q1 += np.sum(Rm[vm>0,:][:,vm>0]) + np.sum(Rm[vm<0,:][:,vm<0])
            if Q0 >= Q1: # Si al partir obtuvimos un Q menor, devolvemos la última partición que hicimos
                return([[ni for ni,vi in zip(nombres_s,v) if vi>0],[ni for ni,vi in zip(nombres_s,v) if vi<0]])
            else:
                # Sino, repetimos para los subniveles
                return(...)

def test_metpot1():
    import numpy as np
    A = np.array([
   [0, 1, 1, 1, 0, 0, 0, 0],
   [1, 0, 1, 1, 0, 0, 0, 0],
   [1, 1, 0, 1, 0, 1, 0, 0],
   [1, 1, 1, 0, 1, 0, 0, 0],
   [0, 0, 0, 1, 0, 1, 1, 1],
   [0, 0, 1, 0, 1, 0, 1, 1],
   [0, 0, 0, 0, 1, 1, 0, 1],
   [0, 0, 0, 0, 1, 1, 1, 0]
])
    
    autovector, autovalor = metpotI2(A, 0.1)
    print("Autovector: ", autovector)
    print("Autovalor: ", autovalor)

def test_deflaciona():
    A = np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    autovalores, autovectores = np.linalg.eig(A)
    print("Autovalores: ", np.round(autovalores, 10))
    print("Autovectores: ", np.round(autovectores, 10))
    deflA = deflaciona(A)
    autovalores_defl, autovectores_defl = np.linalg.eig(deflA)
    print("Autovalores defl: ", np.round(autovalores_defl, 10))
    print("Autovectores defl: ", np.round(autovectores_defl, 10))

if __name__ == "__main__":
    test_deflaciona()