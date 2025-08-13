# [De la web a la ciudad: ranking de museos con PageRank](/01-Museos-en-red/)

Imaginemos por un momento que somos un turista curioso durante “La Noche de los Museos” en la Ciudad de Buenos Aires. Con tantas opciones por delante, ¿cómo decidir cuál visitar primero? ¿El más famoso? ¿El más cercano? ¿O el que nos recomendaron más personas?
En este [trabajo](/01-Museos-en-red/TP1_template.ipynb) vamos a usar una herramienta nacida en el corazón de los buscadores web para responder a una pregunta muy parecida: ¿cómo ordenar lugares según su “importancia”?

A partir de la idea detrás del algoritmo [PageRank](https://es.wikipedia.org/wiki/PageRank) —sí, el mismo que alguna vez hizo famoso a Google— vamos a reemplazar páginas web por museos, y enlaces por cercanía geográfica, para así crear un ranking que nos ayude a explorar la ciudad de manera más “matemática” (y un poco más divertida).


# [De rankings a comunidades: explorando la red de museos de CABA](/02-Comunidades-museos/)

En la exploración de la red de museos de Buenos Aires, pasamos de conocer cuáles son los más destacados a buscar algo más profundo: descubrir comunidades, es decir, grupos de museos más conectados entre sí que con el resto.

Queremos identificar las “pandillas” de museos que comparten vínculos más estrechos. Para lograrlo, nos apoyaremos en dos enfoques clásicos y muy distintos: el corte mínimo, que intenta separar la red rompiendo la menor cantidad posible de conexiones, y la modularidad, que compara la estructura observada con la que esperaríamos por azar, premiando las divisiones que concentran más enlaces dentro de los grupos.

Detrás de estos conceptos hay herramientas matemáticas potentes: la matriz laplaciana, la matriz de modularidad y sus autovalores y autovectores. En otras palabras, utilizaremos álgebra lineal para identificar, de forma automática, las fronteras naturales entre comunidades. Sin embargo, aunque las matemáticas sean precisas, las comunidades reales no siempre lo son; a veces ambos métodos coincidirán y otras veces diferirán, y ahí es donde el análisis y la interpretación cobran protagonismo.

No nos limitaremos a una simple división en dos mitades. También aplicaremos la partición iterativa, que permite detectar múltiples comunidades repitiendo el proceso dentro de cada grupo. De esta forma, podremos descubrir estructuras más complejas y comparar cómo varían según el criterio utilizado.

Finalmente, pondremos todo en práctica sobre la red real de museos, probando distintos parámetros y observando cómo cambian las comunidades detectadas. El [resultado](/02-Comunidades-museos/TP1_template.ipynb) será una radiografía matemática de la ciudad, donde cada grupo estará unido no por afinidad artística o temática, sino por la estructura invisible de sus conexiones.