# Bioinformatics-statistics

En enero del 2023 realizamos, junto a un equipo interdisciplinario, una publicación bioinformática en una reconocida editorial académica:

Fernando Javier Ureta Suelgaray, Viviana Mónica Chiocchio, Federico Ciolfi, Mario Carlos Nazareno Saparrat. Are dark septate endophytes an ancestral ecological state in the evolutionary history of the order Chaetothyriales?. January 2023 Archives of Microbiology 205(2) DOI: 10.1007/s00203-023-03401-6

En este trabajo estuve a cargo del análisis estadístico-bioinformático que detallo a continuación:

Para evaluar la posición filogenética de los DSE dentro del orden de los Chaetothyriales, se analizaron los loci ITS de 115 secuencias correspondientes a especies con distintos hábitos ecológicos. Posteriormente se realizó un alineamiento múltiple de las secuencias utilizando el servidor online de MAFFT (https://mafft.cbrc.jp/alignment/server/) y se optimizaron manualmente utilizando MEGA v10.2.6. Los modelos de evolución de las secuencias de ADN para cada partición de locus, se seleccionaron con jModelTest v.2.1.10 (Darriba et al. 2012) y utilizando el criterio de información de Akaike (AIC, Akaike 1974). El análisis de filogenia se realizó mediante el algoritmo de máxima verosimilitud utilizando el software MEGA v10.2.6. Para ello, se aplicó el modelo GTR con distribución gamma y bootstrapping con 1000 repeticiones. 

Para la reconstrucción de caracteres ancestrales se utilizaron las librerías Phytools (Revell 2012) y APE (Paridis et. al 2019) del lenguaje de programación R. Se estimaron los caracteres ancestrales utilizando el modelo Mk (Lewis, 2001) con ARD (“All rates different”) mediante la función fitMk. Finalmente, el mapeo estocástico fue realizado mediante la función make.simmap (Huelsenbeck 2003, Bollback 2006) aplicando 1000 simulaciones sobre el modelo ARD. Para estimar la diversidad de linajes a través del tiempo (LTT) se utilizó la función LTT y el estadístico gamma se calculó según Pybus y Harvey (2000).

Adjunto en este repositorio el código en R utilizado para realizar el análisis.

-----
In January 2023, together with an interdisciplinary team, we carried out a bioinformatics publication in a renowned academic publisher:

Fernando Javier Ureta Suelgaray, Viviana Mónica Chiocchio, Federico Ciolfi, Mario Carlos Nazareno Saparrat. Are dark septate endophytes an ancestral ecological state in the evolutionary history of the order Chaetothyriales?. January 2023 Archives of Microbiology 205(2) DOI: 10.1007/s00203-023-03401-6

In this work I was in charge of the statistical-bioinformatic analysis that I detail below:

To evaluate the phylogenetic position of the DSEs within the order of Chaetothyriales, the ITS loci of 115 sequences corresponding to species with different ecological habits were analyzed. A multiple alignment of the sequences was subsequently performed using the MAFFT online server (https://mafft.cbrc.jp/alignment/server/) and they were manually optimized using MEGA v10.2.6. The evolutionary models of the DNA sequences for each locus partition were selected with jModelTest v.2.1.10 (Darriba et al. 2012) and using the Akaike information criterion (AIC, Akaike 1974). The phylogeny analysis was performed using the maximum likelihood algorithm using the MEGA v10.2.6 software. For this purpose, the GTR model with gamma distribution and bootstrapping with 1000 repetitions was applied.

For the reconstruction of ancestral characters, the Phytools (Revell 2012) and APE (Paridis et. al. 2019) libraries of the R programming language were used. The ancestral characters were estimated using the Mk model (Lewis, 2001) with ARD (“All rates different”) using the fitMk function. Finally, the stochastic mapping was performed using the make.simmap function (Huelsenbeck 2003, Bollback 2006) applying 1000 simulations on the ARD model. To estimate the diversity of lineages over time (LTT), the LTT function was used and the gamma statistic was calculated according to Pybus and Harvey (2000).

I attach in this repository the R code used to perform the analysis
