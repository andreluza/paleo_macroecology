install.packages('chronosphere')
install.packages('divDyn')
library(chronosphere)
###################Parte 1, plotar mapas shapefile da paleogeografia do ambiente####

#Copie o endereço da pasta onde estão seus dados com a barra invertida (esse é o endereço da pasta onde coloquei os dados no meu pc):
setwd("D:/Unipampa/Permian-Triassic/Data Review - Harvard")

####Explorando a ocorrencia dos fósseis no espaço paleogeográfico####

# ler dados de ocorrencia (planilha gerada no tutorial)
##planilha com registros e coordenadas atuais, e intervalo de tempo
dataocc <- read.csv(file = 'pbdb_data_pareiasuariaocc.csv',sep=";", h=T) 

#confira a planilha
head(dataocc)

#verifique os nomes das colunas
colnames(dataocc)

#selecionar os tempos max e min e calcular a média
tempo <- as.data.frame(rowMeans(dataocc[,16:17]))
##adicionar nova coluna
fossilst <- cbind(dataocc, tempo)
##Confira:
head(fossilst)
##agora deu certo, mas essa coluna nao tem nome, então vamos criar um nome:
##
colnames(fossilst)[which(names(fossilst) == "rowMeans(dataocc[, 16:17])")] <- "mid"

#####comando para baixar e instalar o conjunto de dados do PALEOMAP
#usaremos a variável DEM que corresponde um modelo de elevação da paisagem
dems <- fetch(dat="paleomap", var="dem") 
head(dems)

#nova planilha incluindo apenas os compartimentos necessários -tempos (Visualizar (dezenas) para verificar o compartimento de tempo vs. Ma)
interColl <- fossilst[fossilst$mid%in%c(265, 260, 255),
                      c("accepted_name","lng", "lat", "time", "paleolat" ,"paleolng")]
head(interColl)

#Remova duplicatas, adicione idade do mapa aos dados e reconstrua as paleocoordenadas
unique(interColl$time)

#reconstruct paleocoordinates for collections
#converter as paleocoordenadas para o período específico e o modelo usado 
#reconstruir cada linha para a idade correspondente no mapa de vetor ``age``
newCoords <- reconstruct(interColl[, c("lng", "lat")], 
                         age=interColl[, "time"], 
                         enumerate=FALSE, verbose=FALSE) 

colnames(newCoords) <- c("plng", "plat") #rename columns
allColl <- cbind(interColl,newCoords) #combine vectors
head(allColl)


# collections in each bin
coll270 <- allColl[allColl$time==270 ,] 
coll260 <- allColl[allColl$time==260 ,] 
coll255  <- allColl[allColl$time==255 ,]

pr <- fetch(dat="paleomap", var="paleoatlas", res=0.5) #download and install paleomap dataset

#accessing the different elements corresponding to each age separately
pa260 <- pr["260", ] 
pa255 <- pr["255", ]


par(mfrow=c(1,2))
##para diferenciar os registros de uma espécie, p.ex.
prove <- coll260[coll260$accepted_name == "Provelosaurus americanus",]

#Exportar figura
tiff("provelosaurusmap260final.tiff", width = 6, height = 4, units = 'in', res = 600)
mapplot(pa260,  axes = TRUE, box = TRUE, rgb=TRUE) #Map at 260 Ma
points(coll260$plng, coll260$plat, pch = 21, bg = "black", col = "darkgrey", cex=0.7)
points(prove$plng, prove$plat, pch = 21, bg = "red", col = "white", cex=0.8)
dev.off()


###########################Parte 2 plotar mapas com dados raster da paleotopografia (PALEODEM) do ambiente####
#paleomaps
pr <- fetch(dat="paleomap", var="paleoatlas", res=0.5) #download and install paleomap dataset

#accessing the different elements corresponding to each age separately
pa270 <- pr["270", ] 
pa265 <- pr["265", ] 
pa260 <- pr["260", ] 

par(mfrow=c(2,2))
mapplot(pa270, rgb=TRUE, legend.title=TRUE)
points(coll270$plng, coll270$plat,  pch=19, col="black")

mapplot(pa260, rgb=TRUE) #Map at 260 Ma
points(coll260$plng, coll260$plat,  pch=19, col="black")

mapplot(pa270, rgb=TRUE, legend.title=TRUE) #Map at 270 Ma
points(coll265$plng, coll265$plat,  pch=19, col="darkred")

mapplot(pa260, rgb=TRUE) 
points(coll265$plng, coll265$plat,  pch=19, col="darkred")

############# raster paleoDEM
library(raster)
###Baixar dados de Scotese & Wright 2018
#Copie o endereço da pasta onde estão os dados raster com a barra invertida (esse é o endereço da pasta onde coloquei os dados no meu pc):
###Ou na janela à  direita e abaixo clique em Files, navegue até a pasta onde estao seus dados e depois no ícone de configurações nesta janela, clique em Set As Working Directory
setwd("D:/Unipampa/Permian-Triassic/Paleoenvironment/PaleoDEMlayers/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2")

#read raster
Roadian <- brick("Map52_PALEOMAP_1deg_Middle_Permian_270Ma.nc")
plot(Roadian)
pal <- colorRampPalette(c("darkblue","blue4","blue","cyan3","cyan","chartreuse3","chartreuse4","yellow","orange","red"))
Roadian[Roadian  < -2000] <- -2000
Roadian[Roadian  > 2000] <- 2000

plot(Roadian, col=pal(100))
points(coll270$plng, coll270$plat,  pch = 21, bg = "black", col = "white", cex=0.8)


Wordian <- brick("Map51.5_PALEOMAP_1deg_Middle_Permian_265Ma.nc")
Wordian[Wordian  < -2000] <- -2000
Wordian[Wordian  > 2000] <- 2000

plot(Wordian, col=pal(100))
points(coll265$plng, coll265$plat,  pch = 21, bg = "black", col = "white", cex=0.8)

Capitanian <- brick("Map51_PALEOMAP_1deg_Middle_Permian_260Ma.nc")

Capitanian[Capitanian  < -2000] <- -2000
Capitanian[Capitanian > 2000] <- 2000

plot(Capitanian, col=pal(100))

       maxValue(Roadian)
Roadian[Roadian  < -2000] <- -2000
Roadian[Roadian  > 2000] <- 2000

#save figures
tiff("PaleoDEMprov.tiff", width=200, height=200, bg="white", res=600, units="mm")
pal <- colorRampPalette(c("darkblue","blue4","blue","cyan3","cyan","chartreuse3","chartreuse4","yellow","orange","red"))
plot(Capitanian, col=pal(100))
points(prove$plng, prove$plat,  pch = 21, bg = "black", col = "white", cex=0.8)
dev.off()
