#Script tarea 4

#Cargar la gráfica de la red del club de Karate de Zachary
library(igraph)
karate<-make_graph("Zachary")
plot(karate)

#numero de nodos
nodos<- vcount(karate)
nodos

#calcular la densidad de la red
edge_density(karate)

#matriz de adyacencia
matriz<- as_adjacency_matrix(karate, sparse = FALSE)
matriz

#distribución de conectividades
degree<- degree(karate)
dd<- degree_distribution(karate)
tabla<- table (degree)
sum(tabla)
frecrelative<- tabla/sum(tabla)
barplot(frecrelative, xlab = "Grado", ylab = "Frecuencia relativa", 
        main ="Distribución de conectividades",
        col = "pink", 
        ylim = c(0, max(frecrelative+0.1)))

#diametro, matriz de distancias y distancia promedio
diameter(karate)

distances(karate)
matrizdist<- as.matrix(distances(karate))
matrizdist

mean_distance(karate)

#trayectoria de los nodos mas lejanos
max<- max(distances(karate))
alejados<- which(matrizdist == max, arr.ind = TRUE)
alejados

#medidas de centralidad
centr_degree(karate) #centralidad de grado
centr_betw(karate) #centralidad de intermediación
centr_clo(karate) #centralida de cercanía

#función que calcule el número de amigos que tiene cualquier persona arbitraria 
#y lo compare con el número de los amigos de estos.

g1<-read.csv("C:/Users/xymen/OneDrive/Documentos/GeneticaF/adjacency_matrix.csv")
View(g1)
rownames(g1)<-g1[,1]
View(g1)
g1<-g1[,-1]
g1<-as.matrix(g1)
diag(g1)<-0
red<-graph_from_adjacency_matrix(g1)
plot(red,vertex.size=15,vertex.size=5,
     edge.arrow.size=0.25,layout=layout_nicely,vertex.size.label=0.25)

amigos <- function(red) {
  nodo <- sample(V(red)$name, 1)  
  aminodo <- neighbors(red, nodo, mode = "out")  
  tamañoamigos <- length(aminodo)
  amiamigos <- c()
  for (i in seq_along(aminodo)) {
    principal <- aminodo[i]  
    amideamigos <- neighbors(red, principal, mode = "out")
    tamañoamideamigo <- length(amideamigos)
    amiamigos <- c(amiamigos, tamañoamideamigo)
  }
  
cat("Nombre elegido:", nodo, "\n",
      "Su número de amigos es:", tamañoamigos, "\n",
      "El número de amigos de sus amigos es: ", amiamigos, "\n")
}
amigos(red)

#trayectoria entre dos nodos arbitrarios
trayectoria<- function(red){
  matrizdist <- distances(red)
  max<- max(distances(red))
  alejados<- which(matrizdist == max, arr.ind = TRUE)
  nodoalejado1 <- rownames(matrizdist)[alejados[1, 1]]
  nodoalejado2 <- colnames(matrizdist)[alejados[1, 2]]
  cat("nodos mas alejados: ", nodoalejado1,  ",",  nodoalejado2)
}

trayectoria(red)

#popularidad
degree<- degree(red)
ordenar<- order(degree, decreasing = TRUE)[1:2]
ordenar<-(V(red)$name[ordenar])
cat("las personas más populares son :", ordenar[1], "y", ordenar[2])

#libreria levadura
install.packages("igraphdata")
library(igraphdata)
data("yeast")

#distribucion de conectividades
degree<- degree(yeast)
degredistri<- degree_distribution(yeast)
tabla<- table (degree)
sum(tabla)
frecrelative<- tabla/sum(tabla)
barplot(frecrelative, xlab = "Grado", ylab = "Frecuencia relativa", 
        main ="Distribución de conectividades",
        col = "purple", 
        ylim = c(0, max(frecrelative+0.1)))

#10 proteinas más conectadas
degree<- degree(yeast)
conectadas<- sort(degree, decreasing = TRUE)[1:10]
conectadas

#diámetro y promedio de las distancias
diameter(yeast)
mean_distance(yeast)

#coeficiente de clusterización
transitivity(yeast)

#Cargue y genere una gráfica en `igraph`.
bnmouse<- read.table("C:/Users/xymen/OneDrive/Documentos/tarea4/bn-mouse_visual-cortex_1.edges ")
View(bnmouse)
grafica <- graph_from_edgelist(as.matrix(bnmouse), directed = FALSE)
plot(grafica)

#medidas de centralidad
centr_degree(grafica) #centralidad de grado
centr_betw(grafica) #centralidad de intermediación
centr_clo(grafica) #centralida de cercanía

#densidad
edge_density(grafica)

#Simula una matriz de expresión con 20 genes en 6 condiciones diferentes.
genes <- matrix(rnorm(20 * 6), nrow = 20, ncol = 6)
rownames(genes) <- paste0("Gen", 1:20)
colnames(genes) <- paste0("Cond", 1:6)
genes

#Calcula la correlación entre todos los pares de genes.
cor(genes)

#coexpresion
cor<- cor(genes)
correlacionyumbral<- cor[abs(cor)>=0.8]
coexpresion<- graph.adjacency(correlacionyumbral,  mode = "directed")
plot(coexpresion)

#red aleatoria tipo Erdos-Rényi.
nodos<- 70
conexion<- 3* nodos
erdos<- sample_gnm(nodos , conexion)
plot(erdos)

nodos<- 70
conexion<- 3
ba<-sample_pa(nodos, conexion, directed = FALSE)
plot(ba)

#Compara su grado promedio, distribución de grados, coeficiente de clusterización y diámetro.
mean(degree(ba))
mean(degree(erdos))

degree_distribution(ba)
degree_distribution(erdos)

transitivity(ba)
transitivity(erdos)

diameter(ba)
diameter(erdos)

#Interpreta las diferencias en el contexto de redes biológicas.

