### REBECA VEGA ESCAMILLA ###

#2.A partir de un objeto de tipo phyloseq generado de un análisis  de 
#identificación taxonómica  a partir del gen 16S ribosomal elabora un programa 
#que

#primero cargo mis librerias
library(phyloseq)
library(ggplot2)
library(vegan)
library(microbiome)
library(tibble)
library(dada2)
library(ampvis2)
install.packages("remotes")
remotes::install_github("MadsAlbertsen/ampvis2")

# creo que estos análisis se hacen con los datos que pediste


otu_mat<-data("GlobalPatterns")
tax_mat<-data("GlobalPatterns")
samples_df <-data("GlobalPatterns")


#Los objetos Phyloseq deben tener nombres de fila.
#definir los nombres de las filas de la columna otu
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 
#Ídem para las otras dos matrices
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 
#Transformar en matrices otu y tablas de impuestos
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#Transformar a objetos phyloseq
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

gp <- phyloseq(OTU, TAX, samples)
gp

t(otu_table(GlobalPatterns)) -> otu_table(GlobalPatterns)
otutable <- data.frame(OTU = rownames(phyloseq::otu_table(GlobalPatterns)@.Data),
                       phyloseq::otu_table(GlobalPatterns)@.Data,
                       phyloseq::tax_table(GlobalPatterns)@.Data,
                       check.names = FALSE)

metadata <- data.frame(phyloseq::sample_data(GlobalPatterns), 
                       check.names = FALSE)
#1. Calcule distintas medidas de diversidad

plot_richness(GlobalPatterns, measures=c("Chao1", "Shannon"))
#medidas de diversidad alfa
alpha(GlobalPatterns)

#2. Elabore una gráfica de barras de abundancias por muestras.

barplot(GlobalPatterns)

#3. Elabore un análisis de reducción de dimensionalidad
#Creo que eso lo haría con un diagrama de venn

#4.  Muestre el microbioma core de las muestras
av2 <- amp_load(otu_table, metadata)
amp_core(GlobalPatterns)
#5. (Opcional) genere redes de co-abundacia por muestra.
plot_net(GlobalPatterns)
