from sage.all import *  # Importa todas las funcionalidades de SageMath

#####################################################################################

def vector_puntos(P,Q):
    
    # Input: puntos P y Q
    # Output: vector PQ
    
    return (P.augment(Q))*matrix([[-1],[1]])

def producto_escalar(U,V):
    
    # Input: vectores U y V
    # Output: escalar U*V = U^t·V    
    
    return (transpose(U)*V)[0,0]

def producto_vectorial(U,V):
    
    # Input: vectores U y V 
    # Output: vector producto vectorial UxV   
    
    W0 = U[1,0]*V[2,0] - U[2,0]*V[1,0]
    W1 = -U[0,0]*V[2,0] + U[2,0]*V[0,0]
    W2 = U[0,0]*V[1,0] - U[1,0]*V[0,0]
    return matrix([[W0],[W1],[W2],[0]])

def modulo(U):
    
    # Input: vector U 
    # Output: escalar ||U||=sqrt(U^t·U)    
    
    return sqrt(producto_escalar(U,U))

#####################################################################################

def dibujar_poligono(vertices,**kwds):
    
    # Input: matriz con los vértices del polígono por columnas y en coordenadas homogéneas
    # Output: gráfica del polígono con el borde negro de grosor 1
    
    V = list(vertices[0:3].transpose())
    dib = polygon(V,**kwds)
    V.extend(vertices[0:3,0].transpose())
    dib += line3d(V, color ='black', thickness=1)
    return dib

def dibujar_mallado_poligonal(vertices,caras,**kwds):
    
    # Input: matriz con los vértices del polígono por columnas y en coordenadas homogéneas y matriz de caras.
    # Output: gráfica del mallado.
    
    mallado = 0
    for cara in caras:        
        indices = [n for n in cara if n>-1]
        vertices_de_la_cara = vertices[:,indices]
        mallado = mallado + dibujar_poligono(vertices_de_la_cara,**kwds)
        
    return mallado

#####################################################################################

def matriz_M(vertices, caras):
    
    num_vertices = vertices.ncols()
    M = matrix(ZZ,num_vertices,num_vertices,0) 
    for c in caras:
        indices = [v for v in c if v>-1] 
        for n in range(len(indices)):
            i = indices[n-1]
            j = indices[n]
            M[i,j]+=1 #sumamos 1 en M(i,j) ya que (i,j) es arista del mallado            
    return M

def aristas_recorrido_repetido(M):
    
    aux = M - M.transpose()
    no_nula = [(i, j) for i in range(aux.nrows()) for j in range(aux.ncols()) if aux[i, j] != 0 and i<j]
    
    return no_nula

#####################################################################################

def num_v_a_c(vertices,caras):
    
    # Input: matrices de vértices y caras de un mallado consistente que encierra una cavidad.
    # Output: vértices, aristas, caras y género de la superficie cerrada que delimita.
    
    num_v = vertices.ncols()
    num_c = caras.nrows()
    num_a = 0
    for cara in caras:
        indices = [n for n in cara if n>-1]
        num_a += len(indices)
    num_a = num_a/2
    return [num_v,num_a,num_c]

#####################################################################################

def leer_objeto_obj(path):

    # Input: Ruta donde se encuentra almacenado el archivo .obj
    # Output: Matrices vertices = coordenadas de los vértices por columna y 
    #                  caras = relaciones de los vértices que forman las caras
    
    vertices = []
    caras = []
    
    # Abrir el archivo .obj
    with open(path, 'r') as f:
        for line in f:
            # Leer los vértices
            if line.startswith('v '):
                # Extraemos las coordenadas del vértice (x, y, z)
                v = line.strip().split()[1:]
                # Convertir cada valor de la coordenada flotante a un número racional
                v_racional = [Rational(float(coord)) for coord in v]  # Convertimos a flotante primero y luego a racional
                v_racional.append(1)  # Añadimos la coordenada homogénea (1)
                vertices.append(v_racional)
            
            # Leer las caras
            elif line.startswith('f '):
                # Extraemos los índices de los vértices (restamos 1 porque los índices en .obj empiezan en 1)
                cara = [int(i.split('/')[0]) - 1 for i in line.strip().split()[1:]]
                caras.append(cara)
                
    M = max([len(c) for c in  caras])            
    caras = [c + [-1] * (M - len(c)) for c in caras]
    
    # Convertir las listas a matrices de SageMath
    caras = matrix(len(caras), M,[item for sublist in caras for item in sublist]) 
    vertices = matrix(len(vertices),4,[vertices[i][j] for i in range(len(vertices)) for j in range(4)]).transpose()

    return vertices, caras
    
#####################################################################################

def baricentro(vert): 

    B = matrix([[sum(vert.row(0))],[sum(vert.row(1))],[sum(vert.row(2))],[sum(vert.row(3))]])/sum(vert.row(3))
    
    return B

def caras_arista(i,j,C):
    
    return [k for k in range(C.nrows()) if i in C.row(k) and j in C.row(k)]

def caras_vertice(i,C):
    
    return [k for k in range(C.nrows()) if i in C.row(k)]

def vertices_adyacentes(i, C):
    adyacentes = set()
    for cara in C:
        for k in range(len(cara)):
            if cara[k] == i:
                adyacentes.add(cara[(k-1) % len(cara)])
                adyacentes.add(cara[(k+1) % len(cara)])
    adyacentes.discard(i)
    return list(adyacentes)

#####################################################################################

def calcular_incidentes(matriz_caras):
    
    # Diccionario para almacenar aristas incidentes en cada vértice
    aristas_por_vertice = {}

    # Diccionario para almacenar caras incidentes en cada vértice
    caras_por_vertice = {}

    # Diccionario para almacenar caras incidentes en cada arista
    caras_por_arista = {}

    # Recorrer cada cara
    for i, cara in enumerate(matriz_caras):
        cara = [n for n in cara if n>-1]
        n = len(cara)
        for j in range(n):
            # Obtener los vértices de la arista (v1, v2)
            v1, v2 = sorted((cara[j], cara[(j + 1) % n]))  # Ordenar para evitar duplicados
            arista = (v1, v2)

            # Añadir la arista al vértice v1 y v2
            for v in (v1, v2):
                if v not in aristas_por_vertice:
                    aristas_por_vertice[v] = set()
                aristas_por_vertice[v].add(arista)

            # Añadir la cara a la arista
            if arista not in caras_por_arista:
                caras_por_arista[arista] = []
            caras_por_arista[arista].append(i)

        # Añadir la cara a cada vértice en la cara
        for v in cara:
            if v not in caras_por_vertice:
                caras_por_vertice[v] = set()
            caras_por_vertice[v].add(i)

    # Convertir los conjuntos en listas para facilidad de uso
    for vertice in aristas_por_vertice:
        aristas_por_vertice[vertice] = list(aristas_por_vertice[vertice])
    for vertice in caras_por_vertice:
        caras_por_vertice[vertice] = list(caras_por_vertice[vertice])

    return aristas_por_vertice, caras_por_vertice, caras_por_arista

def puntos_cara(V,C):
    
    punto_por_cara = {}
    
    for i,cara in enumerate(C):        
        indices = [n for n in cara if n>-1]
        vertices_de_la_cara = V[:,indices]
        B = baricentro(vertices_de_la_cara)
        punto_por_cara[i] = B
        
    return punto_por_cara

def puntos_arista(V,caras_por_arista,pc):
    
    punto_por_arista = {}
    
    for arista, caras in caras_por_arista.items(): # caras_por_arista = [(0, 1): [0, 2]]
        v1 = V[:,arista[0]]
        v2 = V[:,arista[1]]
        v3 = pc[caras[0]]
        v4 = pc[caras[1]]
        M = block_matrix([[v1, v2, v3, v4]])
        B = baricentro(M)
        punto_por_arista[arista] = B
    
    return punto_por_arista

def vertices_modificados(V,pc,aristas_por_vertice,caras_por_vertice):
    
    #Input: V = matriz de vértices, pc = output de la función puntos_cara, 
    #       aristas_por_vertice,caras_por_vertice = outputs de la función calcular_incidentes
    #Output: matriz donde la columna i-ésima representa el vértice i-ésimo modificado
    
    Vmod = {}
    
    for v in aristas_por_vertice.keys(): 
        
        # Obtener las aristas incidentes en el vértice
        aristas_v = aristas_por_vertice[v]
        n = len(aristas_v)
        puntos_medios = Matrix(QQ, 4, 0)
        for e in aristas_v:
            e1 = V[:,e[0]]
            e2 = V[:,e[1]]
            pm = (e1+e2)/2
            puntos_medios = puntos_medios.augment(pm)
        m2 = baricentro(puntos_medios)
    
        # Obtener las caras incidentes en el vértice
        caras = caras_por_vertice[v]
        puntos_cara = [pc[c] for c in caras]
        puntos_cara = matrix(QQ, [[col[i, 0] for col in puntos_cara] for i in range(puntos_cara[0].nrows())])
             
        m1 = baricentro(puntos_cara)
        
        Vmod[v] = (m1+2*m2+(n-3)*V[:,v])/n
        
    return Vmod

#suponemos caras cuadrangulares

def caras_nuevas(C,pc,pa,vm):
       
    C2 = []

    for i, cara in enumerate(C):
        p_cara = pc[i] #punto-cara de la caraza
        cara = [n for n in cara if n>-1]
        
        #vértices de la cara modificados
        vm0 = vm[cara[0]]
        vm1 = vm[cara[1]]
        vm2 = vm[cara[2]]
        vm3 = vm[cara[3]]      
        
        #puntos-arista de las aristas adyacentes a la cara
        u1, u2 = sorted((cara[0], cara[1]))
        arista = (u1, u2)
        pa0 = pa[arista]
        
        u1, u2 = sorted((cara[1], cara[2]))
        arista = (u1, u2)
        pa1 = pa[arista]
        
        u1, u2 = sorted((cara[2], cara[3]))
        arista = (u1, u2)
        pa2 = pa[arista]
        
        u1, u2 = sorted((cara[3], cara[0]))
        arista = (u1, u2)
        pa3 = pa[arista]
        
        #nuevas caras
        c0 = Matrix(QQ, 4, 0)
        c0 = c0.augment(vm0).augment(pa0).augment(p_cara).augment(pa3)
        c1 = Matrix(QQ, 4, 0)
        c1 = c1.augment(pa0).augment(vm1).augment(pa1).augment(p_cara)
        c2 = Matrix(QQ, 4, 0)
        c2 = c2.augment(p_cara).augment(pa1).augment(vm2).augment(pa2)
        c3 = Matrix(QQ, 4, 0)
        c3 = c3.augment(pa3).augment(p_cara).augment(pa2).augment(vm3)          
              
        C2.append(c0)
        C2.append(c1)
        C2.append(c2)
        C2.append(c3)

    return C2

def vertices_nuevos(C2):
    
    # Combinar todas las columnas en un conjunto para evitar duplicados    
    vertices_unicos = set(tuple(col) for M in C2 for col in M.columns())
    
    # Crear la matriz final con las columnas únicas
    V2 = matrix(vertices_unicos).transpose()  
        
    return V2

def catmull_clark(V,C):
    
    aristas_por_vertice, caras_por_vertice, caras_por_arista = calcular_incidentes(C)
    pc = puntos_cara(V,C)
    pa = puntos_arista(V,caras_por_arista,pc)
    vm = vertices_modificados(V,pc,aristas_por_vertice,caras_por_vertice)
    
    C2 = caras_nuevas(C,pc,pa,vm) #lista donde cada elemento es una matriz con los vértices de la nueva cara
    V2 = vertices_nuevos(C2) #matriz de vértices nuevos
    
    Caras_indices = []
    for c in C2:
        posiciones = []
        for col in c.columns():
            pos = V2.columns().index(col)
            posiciones.append(pos)
        Caras_indices.append(posiciones)            
    
    return V2, matrix(QQ,Caras_indices)

def catmull_clark_it(V,C,n):
    
    V2 = V
    C2 = C
    
    i = 0
    
    while i<n:
        V2,C2 = catmull_clark(V2,C2)
        i = i+1
        
    return V2,C2 

###################################################################################

def dibuja_segmento_curva_parametrica(C, a, b, **kwds):
    
    #Input: Curva parametrica C, extremos del intervalo de definición
    #Õutput: Dibujo de la curva finita

    graph = parametric_plot3d((C[0,0], C[1,0], C[2,0]), (a, b), **kwds)  
    
    return graph

###################################################################################

def dibuja_vector(P, V,**kwds):
    
    P_origen = (P[0,0], P[1,0], P[2,0])
    P_final = (P[0,0] + V[0,0], P[1,0] + V[1,0], P[2,0] + V[2,0])
    
    graf = arrow(P_origen, P_final, **kwds)
    
    return graf

###################################################################################

def dibuja_cubica_bezier(PC, **kwds):
    
    t = var('t')
    B = matrix([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])    
    T = matrix([[t^3],[t^2],[t],[1]])
    C = PC*B*T
    graph = dibuja_segmento_curva_parametrica(C, 0, 1, **kwds)
    
    return graph

###################################################################################

def dibuja_poligono_control(PC,**kwds):
    
    
    P = [n(PC.column(i)[0:3]) for i in range(4)]
    
    return points(P,**kwds) + line3d(P,**kwds)

###################################################################################