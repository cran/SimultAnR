################################################################################
#####               SimultAnR (Simultaneous Analysis) Package              #####
################################################################################
### The package includes six functions:
### CorrAn <- function(data, sr = NA, sc = NA, nd = 2, dp = 2)
### SimAn <- function(data, G, acg, weight = 2, nameg = NA, sr = NA, sc = NA,
###    nd = 2, dp = 2)
### CorrAnSummary <- function(input)
### SimAnSummary <- function(input)
### CorrAnGraph <- function(input, s1 = 1, s2 = 2, screen = TRUE)
### SimAnGraph <- function(input, s1 = 1, s2 = 2, screen = TRUE)



################################################################################
#####               CorrAn function: Correspondence Analysis               #####
################################################################################
CorrAn <- function(data, sr = NA, sc = NA, nd = 2, dp = 2)
{
   ### opciones originales
   datos <- as.matrix(data)
   FS <- sr
   CS <- sc
   r <- dp
   ### Empieza
   ### IFS
   if(is.na(FS[1])) {
      IFS <- 0
      I <- nrow(datos)   }
   else {
      IFS <- length(FS)
      I <- nrow(datos) - IFS   }
   ### JCS
   if(is.na(CS[1])) {
      JCS <- 0
      J <- ncol(datos)   }
   else {
      JCS <- length(CS)
      J <- ncol(datos) - JCS   }
   ### Datos activos
   if(IFS != 0) {
      if(JCS != 0) {
         datosA <- datos[ - FS, - CS]   }
      else {
         datosA <- datos[ - FS, ]   }   }
   else {
      if(JCS != 0) {
         datosA <- datos[, - CS]   }
      else {
         datosA <- datos   }   }
   nomi <- labels(datosA)[[1]]
   nomj <- labels(datosA)[[2]]
   ### Filas suplementarias
   datosFS <- matrix(0, IFS, J)
   if(IFS != 0) {
      if(JCS != 0) {
         datosFS[, ] <- as.matrix(datos[FS, - CS])   }
      else {
         datosFS[, ] <- datos[FS, ]   }   }
   else {
      datosFS[, ] <- 0   }
   nomiFS <- labels(datos)[[1]][FS]
   ### Columnas suplementarias
   datosCS <- matrix(0, I, JCS)
   if(JCS != 0) {
      if(IFS != 0) {
         datosCS[, ] <- datos[ - FS, CS]   }
      else {
         datosCS[, ] <- datos[, CS]   }   }
   else {
      datosCS[, ] <- 0   }
   nomjCS <- labels(datos)[[2]][CS]
   ### Analisis activo
   ### Matriz de frecuencias
   Fd <- datosA/sum(datosA)
   fi <- apply(Fd, 1, sum)
   fj <- apply(Fd, 2, sum)
   I <- length(fi)
   J <- length(fj)
   ### Numero de ejes editados
   ndim <- min(I, J, nd)
   fii <- matrix(0, I, J)
   for(i in 1:J) {
      fii[, i] <- fi   }
   fjj <- matrix(0, I, J)
   for(j in 1:I) {
      fjj[j, ] <- fj   }
   P <- diag(fi, I, I)
   M <- diag(fj, J, J)
   X <- (sqrt(fii)) * ((solve(P) %*% Fd %*% solve(M)) - 1) * (sqrt(fjj))
   dimnames(X) <- list(nomi, nomj)
   ### Diagonalizacion de la matriz de menor orden
   if(I >= J) {
      A <- t(X) %*% X   }
   else {
      A <- (X) %*% t(X)   }
   vls <- eigen(A)$values
   intot <- sum(vls)
   names(intot) <- "Total inertia"
   U1 <- eigen(A)$vectors
   u1n <- abs(diag(t(U1) %*% U1))
   U2 <- t(t(U1)/sqrt(u1n))
   vls <- vls[1:ndim]
   fe <- sqrt(solve(diag(abs(vls))))
   if(I >= J) {
      U <- U2[, 1:ndim]
      V <- X %*% U %*% fe   }
   else {
      V <- U2[, 1:ndim]
      U <- t(X) %*% V %*% fe   }
   ### Fin diagonalizacion
   ### Porcentajes de inercia
   vlstod <- eigen(A)$values
   porin <- vlstod/intot
   poracuin <- cumsum(porin)
   ### Proyecciones
   Fs <- solve(sqrt(P)) %*% X %*% U
   dimnames(Fs) <- list(nomi, paste("F", 1:ndim, sep = ""))
   Gs <- solve(sqrt(M)) %*% t(X) %*% sqrt(P) %*% Fs %*% fe
   dimnames(Gs) <- list(nomj, paste("G", 1:ndim, sep = ""))
   ### Contribuciones absolutas
   ctai <- t(t((Fs * Fs) * fi)/vls)
   ctaj <- t(t((Gs * Gs) * fj)/vls)
   ### Dist de perfiles al origen
   d2i <- apply((((Fd/fii) - fjj)^2) %*% solve(M), 1, sum)
   for(i in 1:I) {
      if(fi[i] == 0) {
         d2i[i] <- 0   }   }
   d2j <- apply(solve(P) %*% (((Fd/fjj) - fii)^2), 2, sum)
   for(j in 1:J) {
      if(fj[j] == 0) {
         d2j[j] <- 0   }   }
   ### Contribuciones relativas
   ctri <- (Fs * Fs)/d2i
   for(i in 1:I) {
      if(fi[i] == 0) {
         ctri[i, ] <- 0   }   }
   ctrj <- (Gs * Gs)/d2j
   for(j in 1:J) {
      if(fj[j] == 0) {
         ctrj[j, ] <- 0   }   }
   ### Resultados
   resin <- cbind(Re(round(vlstod, r + 5)), Re(round(100 * porin, r + 5)),
      Re(round(100 * poracuin, r + 5)))
   dimnames(resin) <- list(paste("s=", 1:length(vlstod), sep = ""),
      c("values", "percentage", "cumulated"))
   resi <- cbind(Re(round(100 * fi, r)), Re(round(d2i, r)),
      Re(round(Fs[, 1:ndim], r)), Re(round(100 * ctai[, 1:ndim], r)),
      Re(round(ctri[, 1:ndim], r)))
   dimnames(resi) <- list(nomi, c("100fi", "d2", paste("F", 1:ndim, sep = ""),
      paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
   resj <- cbind(Re(round(100 * fj, r)), Re(round(d2j, r)),
      Re(round(Gs[, 1:ndim], r)), Re(round(100 * ctaj[, 1:ndim], r)),
      Re(round(ctrj[, 1:ndim], r)))
   dimnames(resj) <- list(nomj, c("100fj", "d2", paste("G", 1:ndim, sep = ""),
      paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
   ########## ########## ########## ########## ########## ########## ###########
   ### Filas suplementarias
   if(IFS == 0) {
      FdFS = fiFS = giFS = giiFS = gjjFS = XFS = nomiFS = d2iFS = resicFS =
      resigFS = FsicFS = FsigFS = todoFsFS = FsFS <- 0
      IFS = 0
      resiFS <- "No supplementary rows"   }
   if(IFS != 0) {
      FdFS <- as.matrix(datosFS/sum(datosA))
      fiFS <- apply(FdFS, 1, sum)
      IFS <- length(fiFS)
      giFS <- fiFS
      giiFS <- matrix(0, IFS, J)
      for(i in 1:J) {
         giiFS[, i] <- giFS      }
      gjjFS <- matrix(0, IFS, J)
      for(j in 1:IFS) {
         gjjFS[j, ] <- fj      }
      ### Matriz X Filas Suplementarias
      XFS <- matrix(0, IFS, J)
      for(i in 1:IFS) {
         for(j in 1:J) {
            XFS[i, j] <- (sqrt(giFS[i]) * (FdFS[i, j]/giFS[i] - 
               fj[j]))/sqrt(fj[j])   }   }
      XFS[!is.finite(XFS)] <- 0
      dimnames(XFS) <- list(nomiFS, nomj)
      ### Proyecciones Filas Suplementarias
      if(length(FS) == 1) {
         FsFS <- solve(sqrt(giFS)) %*% XFS %*% U   }
      else {
         FsFS <- solve(sqrt(diag(giFS))) %*% XFS %*% U   }
      dimnames(FsFS) <- list(nomiFS, paste("F", 1:ndim, sep = ""))
      ### Distancias y Contribuciones FS Filas Suplementarias
      d2iFS <- apply((((FdFS/giiFS) - gjjFS)^2) %*% solve(diag(fj)), 1, sum)
      d2iFS[!is.finite(d2iFS)] <- 0
      names(d2iFS) <- (nomiFS)
      ctriFS <- (FsFS^2)/d2iFS
      dimnames(ctriFS) <- list(nomiFS, NULL)
      ## Resultados FS Filas Suplementarias
      resiFS <- Re(round(cbind(100 * fiFS, as.matrix(d2iFS), FsFS, ctriFS), r))
      dimnames(resiFS) <- list(nomiFS, c("100fi", "d2",
         paste("F", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))   }
   ### fin filas suplementarias
   ########## ########## ########## ########## ########## ########## ###########
   #### CS Columnas Suplementarias
   if(JCS == 0) {
      JCS = GsCS = nomjCS = fjCS = d2jCS <- 0
      resjCS <- "No supplementary columns"   }
   if(JCS != 0) {
      FdCS <- datosCS/sum(datosA)
      fjCS <- apply(FdCS, 2, sum)
      JCS <- length(fjCS)
      gjCS <- fjCS
      giiCS <- matrix(0, I, JCS)
      for(i in 1:JCS) {
         giiCS[, i] <- fi   }
      gjjCS <- matrix(0, I, JCS)
      for(j in 1:I) {
         gjjCS[j, ] <- gjCS   }
      ### Matriz X Columnas Suplementarias
      XCS <- matrix(0, I, JCS)
      for(i in 1:I) {
         for(j in 1:JCS) {
            XCS[i, j] <- (FdCS[i, j]/gjCS[j]) - fi[i]   }   }
      ### Proyecciones Columnas Suplementarias
      GsCS <- t(XCS) %*% Fs %*% fe
      dimnames(GsCS) <- list(nomjCS, paste("G", 1:ndim, sep = ""))
      ### Distancias y Contribuciones Columnas Suplementarias
      d2jCS <- apply(solve(diag(fi)) %*% (((FdCS/gjjCS) - giiCS)^2), 2, sum)
      names(d2jCS) <- (nomjCS)
      ctrjCS <- (GsCS^2)/d2jCS
      dimnames(ctrjCS) <- list(nomjCS, NULL)
      ## Resultados CS columnas suplementarias
      resjCS <- Re(round(cbind(100 * fjCS, as.matrix(d2jCS), GsCS, ctrjCS), r))
      dimnames(resjCS) <- list(nomjCS, c("100fj", "d2",
         paste("G", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))   }
   ### Fin Columnas Suplementarias
   ########## ########## ########## ########## ########## ########## ###########
   return(list(totalin = intot, eig = vls, resin = resin, resi = resi,
      resj = resj, resisr = resiFS, resjsc = resjCS,
      X = X, totalk = sum(datosA),
      I = I, namei = nomi, fi = fi, Fs = Fs, d2i = d2i,
      J = J, namej = nomj, fj = fj, Gs = Gs, d2j = d2j,
      Isr = IFS, nameisr = nomiFS, fisr = fiFS, Fssr = FsFS, d2isr = d2iFS,
      Xsr = XFS,
      Jsc = JCS, namejsc = nomjCS, fjsc = fjCS, Gssc = GsCS, d2jsc = d2jCS))
}
### Fin CorrAn: Correspondence Analysis
################################################################################



################################################################################
#####               SimAn function: Simultaneous Analysis                  #####
################################################################################
SimAn <- function(data, G, acg, weight = 2, nameg = NA, sr = NA, sc = NA,
   nd = 2, dp = 2)
{
   ### Input
   pond <- weight
   r <- dp
   datos <- as.matrix(data)
   FS <- sr
   CAg <- acg
   CS <- sc
   nomg <- nameg
   ### IFS
   if(is.na(FS[1])) {
      IFS <- 0
      I <- nrow(datos)
      nomic <- labels(datos)[[1]]   }
   else {
      IFS <- length(FS)
      I <- nrow(datos) - IFS
      nomic <- labels(datos)[[1]][ - FS]
      nomicFS <- labels(datos)[[1]][FS]   }
   ### Nombres grupos
   if(is.na(nomg[1])) {
      nomg <- paste("G", 1:G, sep = "")   }
   ### Analisis separados
   CAres <- vector("list", G)
   names(CAres) <- paste("Table_", 1:G, sep = "")
   for(g in 1:G) {
      Jg <- length(CAg[g][[1]])
      if(!is.na(CS[1])) {
         tablag <- datos[, c(CAg[g][[1]], CS)]
         CSg <- c((1 + Jg):(ncol(tablag)))   }
      else {
         tablag <- datos[, c(CAg[g][[1]])]
         CSg <- NA   }
      dimnames(tablag)[[1]] <- paste(nomg[g], labels(tablag)[[1]], sep = "")
      CAres[[g]] <- CorrAn(tablag, sr = FS, sc = CSg, nd = nd, dp = r)
      names(CAres[[g]]) <- c("intot", "vls", "resin",
         "resi", "resj", "resiFS", "resjCS", "X", "KtotA",
         "I", "nomi", "fi", "Fs", "d2i",
         "J", "nomj", "fj", "Gs", "d2j",
         "IFS", "nomiFS", "fiFS", "FsFS", "d2iFS", "XFS",
         "JCS", "nomjCS", "fjCS", "GsCS", "d2jCS")
   }
   ### Comienza el Simultaneo
   J <- 0
   for(g in G:1) {
      J <- CAres[[g]]$J + J   }
   ndim <- min(I, J, nd)
   ### Nombres perfiles fila
   nomig <- c()
   for(g in G:1) {
      nomig <- c(CAres[[g]]$nomi, nomig)   }
   ### Nombres perfiles columna
   nomj <- c()
   for(g in G:1) {
      nomj <- c(CAres[[g]]$nomj, nomj)   }
   ### Ponderac tablas: 0=1/1er vls separado; 1=1; 2=1/intot
   alphag <- c()
   if(pond == 1) {
      for(g in G:1) {
         alphag <- c(1, alphag)   }   }
   if(pond == 2) {
      for(g in G:1) {
         alphag <- c(1/CAres[[g]]$vls[[1]], alphag)   }   }
   if(pond == 3) {
      for(g in G:1) {
         alphag <- c(1/CAres[[g]]$intot, alphag)   }   }
   ### Matriz(IxG) raices pesos fila
   rfi <- array(0, c(I, G))
   for(g in G:1) {
      rfi[, g] <- sqrt(CAres[[g]]$fi)   }
   ### Raices pesos compromisos
   srfi <- c(rep(0, I))
   for(g in G:1) {
      srfi <- rfi[, g] + srfi   }
   ### Pesos compromisos (pic); parciales (fiG); columnas (fj)
   pic <- srfi^2
   fiG <- c()
   for(g in G:1) {
      fiG <- c(CAres[[g]]$fi, fiG)   }
   fj <- c()
   for(g in G:1) {
      fj <- c(CAres[[g]]$fj, fj)   }
   ### Matriz X
   X <- c()
   for(g in G:1) {
      X <- cbind(sqrt(alphag[g]) * CAres[[g]]$X, X)   }
   ### Matriz Y para parciales
   cJ <- c()
   for(g in G:1) {
      cJ <- c(CAres[[g]]$J, cJ)   }
   # suma de J(g + 1) para que empiece por 0
   ndimJgm1 <- c(0)
   for(g in G:1) {
      ndimJgm1[g + 1] <- sum(cJ[1:g])   }
   Y <- c()
   Yg <- array(0, c(I, J, G))
   for(g in 1:G) {
      Yg[, (1 + ndimJgm1[g]):ndimJgm1[g + 1], g] <- sqrt(alphag[g]) * 
         CAres[[g]]$X}
   for(g in G:1) {
      Y <- rbind(Yg[, , g], Y)   }
   ### Diagonalizacion de la matriz de menor orden
   if(I >= J) {
      A <- t(X) %*% X   }
   else {
      A <- (X) %*% t(X)   }
   vls <- eigen(A)$values
   intot <- sum(vls)
   U1 <- eigen(A)$vectors
   u1n <- abs(diag(t(U1) %*% U1))
   U2 <- t(t(U1)/sqrt(u1n))
   vls <- vls[1:ndim]
   fe <- sqrt(solve(diag(abs(vls))))
   if(I >= J) {
      U <- U2[, 1:ndim]
      V <- X %*% U %*% fe   }
   else {
      V <- U2[, 1:ndim]
      U <- t(X) %*% V %*% fe   }
   ### Fin diagonalizacion de la matriz de menor orden
   ### Inercias y porcentajes
   vlstod <- eigen(A)$values
   porin <- vlstod/intot
   poracuin <- cumsum(porin)
   ### Proyecciones
   Fsig <- solve(sqrt(diag(fiG, I * G, I * G))) %*% Y %*% U
   dimnames(Fsig) <- list(nomig, paste("F", 1:ndim, sep = ""))
   Fsic <- solve(sqrt(diag(srfi * srfi, I, I))) %*% X %*% U
   dimnames(Fsic) <- list(nomic, paste("F", 1:ndim, sep = ""))
   Gs <- solve(sqrt(diag(fj, J, J))) %*% t(X) %*% V
   dimnames(Gs) <- list(nomj, paste("G", 1:ndim, sep = ""))
   todoFs <- array(0, c(I, ndim, G + 1))
   dimnames(todoFs) <- list(nomic, paste("F", 1:ndim, sep = ""),
      c(paste("Table ", nomg, sep = ""), "overall rows"))
   PYg <- array(0, c(I, J, G))
   for(g in 1:G) {
      PYg[, (1 + ndimJgm1[g]):ndimJgm1[g + 1], g] <- 
         sqrt(alphag[g]) * solve(diag(sqrt(CAres[[g]]$fi))) %*% CAres[[g]]$X
      todoFs[, , g] <- PYg[, , g] %*% U[, 1:ndim]   }
   todoFs[, , G + 1] <- Fsic[, 1:ndim]
   ### Contribuciones absolutas
   ctaic <- t(t((Fsic^2) * (srfi^2))/vls)
   ctag <- c()
   for(g in G:1) {
      Gsg <- ((solve(sqrt(diag(CAres[[g]]$fj))) %*%
         ((sqrt(alphag[g])) * t(CAres[[g]]$X)) %*% V[, 1:ndim]))
      ctajg <- t(t((Gsg * Gsg) * CAres[[g]]$fj)/vls)
      ctag <- rbind(apply(ctajg, 2, sum), ctag)   }
   ctag <- Re(round(ctag, r))
   dimnames(ctag) <- list(paste("Table", nomg[1:G]),
      paste("Axis ", 1:ndim, sep = ""))
   ctaj <- t(t((Gs * Gs) * fj)/vls)
   ### Distancias de puntos al origen
   d2ig <- c()
   for(g in G:1) {
      d2ig <- c((alphag[g]) * CAres[[g]]$d2i, d2ig)   }
   d2ic <- 0
   for(g in G:1) {
      d2ic <- c(((alphag[g]) * (CAres[[g]]$fi/srfi^2) * CAres[[g]]$d2i) +  
         d2ic)   }
   d2j <- c()
   for(g in G:1) {
      d2j <- c((alphag[g]) * matrix(CAres[[g]]$d2j), d2j)   }
   ### Contribuciones relativas
   ctric <- (Fsic^2)/d2ic
   ctrig <- (Fsig^2)/d2ig
   ctrig[!is.finite(ctrig)] <- 0
   ctrj <- (Gs^2)/d2j
   ### Proyecciones de grupos
   Fsg <- t(vls * t(ctag))
   Fsg <- Re(round(Fsg, r))
   # Relaciones compromisos y parciales
   rcp <- matrix(0, G, ndim)
   for(g in 1:G) {
      for(is in 1:ndim) {
         varsg <- t(diag(rfi[, g]) %*% todoFs[, is, g]) %*% (diag(rfi[, g]) %*%
            todoFs[, is, g])
         rcp[g, is] <- t((1/sqrt(varsg[1, 1])) * (diag(rfi[, g]) %*%
            todoFs[, is, g])) %*% ((1/sqrt(vls[is])) * (diag(srfi, I, I) %*%
            todoFs[, is, G + 1]))   }   }
   rcp <- Re(round(rcp, r))
   dimnames(rcp) <- list(paste("Table", nomg[1:G]),
      paste("Axis ", 1:ndim, sep = ""))
   ### Relaciones AC y An Simultaneo
   rcs <- array(0, c(ndim, ndim, G))
   for(g in 1:G) {
      ndimcs <- min((I - 1), ((CAres[[g]]$J) - 1), ndim)
      for(is in 1:ndim) {
         for(sp in 1:ndimcs) {
            if(CAres[[g]]$vls[sp]>0 & vls[is]>0){
               rcs[sp, is, g] <- t((1/sqrt(CAres[[g]]$vls[sp])) * 
                  (diag(rfi[, g]) %*% CAres[[g]]$Fs[, sp])) %*%
                  ((1/sqrt(vls[is])) * (diag(srfi, I, I) %*%
                  todoFs[, is, G + 1]))   }   }   }   }
   dimnames(rcs) <- list(paste("CA Axis", 1:ndim), paste("SA Axis", 1:ndim),
      paste("Table", nomg[1:G]))
   rcs <- Re(round(rcs, r))
   ### Relaciones ejes diferentes AC
   rcc <- array(0, c(ndim, ndim, G * (G - 1) / 2))
   namegh <- c(1:(G * (G - 1) / 2))
   gh <- 1
   for(g in 1:G) {
      ndimis <- min((I - 1), ((CAres[[g]]$J) - 1), ndim)
      for(h in 1:G) {
         ndimsp <- min((I - 1), ((CAres[[h]]$J) - 1), ndim)
         if (g < h){
            for(is in 1:ndimis) {
               for(sp in 1:ndimsp) {
                  if(CAres[[g]]$vls[is]>0 & CAres[[h]]$vls[sp]>0) {
                     rcc[sp, is, gh] <- t((1/sqrt(CAres[[h]]$vls[sp])) * 
                        (diag(rfi[, h]) %*% CAres[[h]]$Fs[, sp])) %*%
                        ((1/sqrt(CAres[[g]]$vls[is])) * (diag(rfi[, g]) %*%
                        CAres[[g]]$Fs[, is]))
                     namegh[gh] <- paste("Table", nomg[g], "Table", nomg[h])
                  }   }   }
            gh <- gh + 1
         }   }   }
   rcc <- Re(round(rcc, r))
   dimnames(rcc) <- list(paste("CA Axis", 1:ndim), paste("CA Axis", 1:ndim),
      namegh)
   ### Para graficas
   maxJg <- 0
   for(g in G:1) {
      maxJg <- max(maxJg, CAres[[g]]$J)   }
   ### Nombres columnas
   Osc <- matrix(0, 1, ndim)
   todoGs <- array(0, c(maxJg, ndim, G))
   for(g in G:1) {
      Gsg <- (solve(sqrt(diag(CAres[[g]]$fj)))) %*% t(sqrt(alphag[g]) * 
         CAres[[g]]$X) %*% V[, 1:ndim]
      todoGs[(1:CAres[[g]]$J), , g] <- Gsg[, 1:ndim]   }
      ### sin nombres porque tiene que ser igual en todas las tablas (1ra
      ### dimension de array)
   ########## ########## ########## ########## ########## ########## ###########
   ### Filas suplementarias
   if(is.na(FS[1])) {
      resicFS = resigFS = FsicFS = FsigFS = todoFsFS = nomicFS <- 
      "No supplementary rows"
      IFS = 0   }
   if(!is.na(FS[1])) {
      IFS <- length(FS)
      fiFS <- c()
      for(g in G:1) {
         fiFS <- c(CAres[[g]]$fiFS, fiFS)   }
      giFS <- fiFS
      nomigFS <- c()
      ### Nombres perfiles Filas Suplementarias
      for(g in G:1) {
         nomigFS <- c(CAres[[g]]$nomiFS, nomigFS)   }
      rfiFS <- array(0, c(IFS, G))
      for(g in G:1) {
         rfiFS[, g] <- sqrt(CAres[[g]]$fiFS)   }
      srfiFS <- c(rep(0, IFS))
      for(g in G:1) {
         srfiFS <- rfiFS[, g] + srfiFS   }
      picFS <- srfiFS^2
      RFS <- diag(srfiFS * srfiFS, IFS, IFS)
      XFS <- c()
      ### Definir en funcion de s
      for(g in G:1) {
         XFS <- cbind(sqrt(alphag[g]) * CAres[[g]]$XFS, XFS)   }
      YgFS <- array(0, c(IFS, J, G))
      for(g in 1:G) {
         YgFS[, (1 + ndimJgm1[g]):ndimJgm1[g + 1], g] <- 
            sqrt(alphag[g]) * CAres[[g]]$XFS   }
      YFS <- c()
      for(g in G:1) {
         YFS <- rbind(YgFS[, , g], YFS)   }
      fiGFS <- c()
      for(g in G:1) {
         fiGFS <- c(CAres[[g]]$fiFS, fiGFS)   }
      QFS <- diag(fiGFS, IFS * G, IFS * G)
      FsicFS <- solve(sqrt(RFS)) %*% XFS %*% U
      dimnames(FsicFS) <- list(nomicFS, paste("F", 1:ndim, sep = ""))
      FsigFS <- solve(sqrt(QFS)) %*% YFS %*% U
      dimnames(FsigFS) <- list(nomigFS, paste("F", 1:ndim, sep = ""))
      ###
      todoFsFS <- array(0, c(IFS, ndim, G + 1))
      dimnames(todoFsFS) <- list(nomicFS, paste("F", 1:ndim, sep = ""),
         c(paste("Table ", nomg, sep = ""), "overall rows"))
      PYgFS <- array(0, c(IFS, J, G))
      for(g in 1:G) {
         if(IFS == 1) {
            PYgFS[, (1 + ndimJgm1[g]):ndimJgm1[g + 1], g] <- sqrt(alphag[g]) * 
               solve(sqrt(CAres[[g]]$fiFS)) %*% CAres[[g]]$XFS   }
         else {
            PYgFS[, (1 + ndimJgm1[g]):ndimJgm1[g + 1], g] <- sqrt(alphag[g]) * 
               solve(diag(sqrt(CAres[[g]]$fiFS))) %*%CAres[[g]]$XFS   }
         todoFsFS[, , g] <- PYgFS[, , g] %*% U[, 1:ndim]   }
      todoFsFS[, , G + 1] <- FsicFS[, 1:ndim]
      ### Distancias
      d2igsFS <- c()
      for(g in G:1) {
         d2igsFS <- c((alphag[[g]]) * CAres[[g]]$d2iFS, d2igsFS)   }
      d2igFS <- matrix(d2igsFS)
      dimnames(d2igFS) <- list(nomigFS, NULL)
      d2icsFS <- 0
      for(g in G:1) {
         d2icsFS <- c(((alphag[[g]]) * (CAres[[g]]$fiFS/srfiFS^2) * 
            CAres[[g]]$d2iFS) + d2icsFS)   }
      d2icFS <- matrix(d2icsFS)
      dimnames(d2icFS) <- list(nomicFS, NULL)
      ctricFS <- (FsicFS^2)/d2icsFS
      dimnames(ctricFS) <- list(nomicFS, NULL)
      ctrigFS <- (FsigFS^2)/d2igsFS
      ctrigFS[!is.finite(ctrigFS)] <- 0
      dimnames(ctrigFS) <- list(nomigFS, NULL)
      ### Resultados FS Filas Suplementarias
      resicFS <- cbind(Re(round(picFS, r)), Re(round(d2icFS, r)),
         Re(round(FsicFS, r)), Re(round(ctricFS, r)))
      dimnames(resicFS) <- list(nomicFS, c("pi", "d2",
         paste("F", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
      resigFS <- cbind(Re(round(100 * giFS, r)), Re(round(d2igFS, r)),
         Re(round(FsigFS, r)), Re(round(ctrigFS, r)))
      dimnames(resigFS) <- list(nomigFS, c("100fig", "d2",
         paste("F", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))   }
   ### Fin Filas Suplementarias
   ########## ########## ########## ########## ########## ########## ###########
   #### Columnas Suplementarias
   if(is.na(CS[1])) {
      todoresjCS = todoGsCS = GsCS = nomjCS <- "No supplementary columns"
      JCS = 0   }
   if(!is.na(CS[1])) {
      ### para todoresjCS, todoGsCS
      JCS <- length(CS)
      datosCS <- matrix(0, I, JCS)
      if(!is.na(FS[1])) {
         datosCS <- as.matrix(datos[ - FS, CS])   }
      else {
         datosCS <- as.matrix(datos[, CS])   }
      nomjCSorig <- labels(datosCS)[[2]]
      todoGsCS <- array(0, c(JCS, ndim, G))
      todoresjCS <- 0
      GsCS <- array(0, c(1, ndim))
      for(g in G:1) {
         nomjCS <- paste(nomg[g], nomjCSorig, sep = "")
         FdCS <- datosCS/CAres[[g]]$KtotA
         fjCS <- apply(FdCS, 2, sum)
         giCS <- CAres[[g]]$fi
         gjCS <- fjCS
         giiCS <- matrix(0, I, JCS)
         for(i in 1:JCS) {
            giiCS[, i] <- giCS   }
         gjjCS <- matrix(0, I, JCS)
         for(j in 1:I) {
            gjjCS[j, ] <- gjCS   }
         XCS <- sqrt(alphag[g]) * (1/(sqrt(giiCS))) * (FdCS/ (gjjCS) - giiCS)
         XCS <- as.matrix(XCS)
         XCS[!is.finite(XCS)] <- 0
         GsgCS <- t(XCS) %*% V[, 1:ndim]
         dimnames(GsgCS) <- list(nomjCS, paste("G", 1:ndim, sep = ""))
         d2jCS <- apply(alphag[g] * solve(diag(giCS)) %*%
            (((FdCS/gjjCS) - giiCS)^2), 2, sum)
         names(d2jCS) <- (nomjCS)
         ctrjCS <- (GsgCS^2)/d2jCS
         dimnames(ctrjCS) <- list(nomjCS, NULL)
         resjCS <- Re(round(cbind(100 * fjCS, as.matrix(d2jCS), GsgCS, ctrjCS),
            r))
         dimnames(resjCS) <- list(nomjCS, c("100fj", "d2",
            paste("G", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
         todoresjCS <- rbind(resjCS, todoresjCS)
         todoGsCS[(1:JCS), , g] <- GsgCS[, 1:ndim]
         GsCS <- rbind(GsgCS, GsCS)
         dimnames(todoGsCS) <- list(CAres[[g]]$nomjCS,
            paste("G", 1:ndim, sep = ""),
            c(paste("Table ", nomg, sep = "")))   }
      todoresjCS <- todoresjCS[1:nrow(todoresjCS) - 1, ]
      GsCS <- GsCS[1:nrow(GsCS) - 1, ]
   }
   ### Fin Columnas Suplementarias
   #######################   resultados GAS   ##################################
   names(CAres) <- paste("Table_", 1:G, sep = "")
   for(g in 1:G) {
      names(CAres[[g]]) <- c("totalin", "eig", "resin", "resi", "resj",
         "resisr", "resjsc",
         "X", "totalk",
         "I", "namei", "fi", "Fs", "d2i",
         "J", "namej", "fj", "Gs", "d2j",
         "Isr", "nameisr", "fisr", "Fssr", "d2isr",
         "Xsr",
         "Jsc", "namejsc", "fjsc", "Gssc", "d2jsc")   }
   resin <- cbind(Re(round(vlstod, r + 5)), Re(round(100 * porin, r + 5)),
      Re(round(100 * poracuin, r + 5)))
   dimnames(resin) <- list(paste("s=", 1:length(vlstod), sep = ""),
      c("values", "percentage", "cumulated"))
   resic <- cbind(Re(round(pic, r)), Re(round(d2ic, r)),
      Re(round(Fsic[, 1:ndim], r)), Re(round(100 * ctaic[, 1:ndim], r)),
      Re(round(ctric[, 1:ndim], r)))
   dimnames(resic) <- list(nomic, c("pi", "d2", paste("F", 1:ndim, sep = ""),
      paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
   resig <- cbind(Re(round(100 * fiG, r)), Re(round(d2ig, r)),
      Re(round(Fsig[, 1:ndim], r)), Re(round(ctrig[, 1:ndim], r)))
   dimnames(resig) <- list(nomig, c("100fig", "d2",
      paste("F", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
   resj <- cbind(Re(round(100 * fj, r)), Re(round(d2j, r)),
      Re(round(Gs[, 1:ndim], r)), Re(round(100 * ctaj[, 1:ndim], r)),
      Re(round(ctrj[, 1:ndim], r)))
   dimnames(resj) <- list(nomj, c("100fjg", "d2", paste("G", 1:ndim, sep = ""),
      paste("ctr", 1:ndim, sep = ""), paste("cor", 1:ndim, sep = "")))
   ########## ########## ########## ########## ########## ########## ###########
   return(list(totalin = intot, resin = resin, resi = resic, resig = resig,
      resj = resj, Fsg = Fsg, ctrg = ctag, riig = rcp, rCACA = rcc, rCASA = rcs,
      Fsi = Fsic, Fsig = Fsig, Gs = Gs, allFs = todoFs, allGs = todoGs,
      I = I, maxJg = maxJg, G = G, namei = nomic, nameg = nomg,
      # namejg = nomjg,
      # d2g = d2g, resg = resg, Lgh = Lgh, d2gh = d2gh,
      resisr = resicFS, resigsr = resigFS, Fsisr = FsicFS,
      Fsigsr = FsigFS, allFssr = todoFsFS,
      Isr = IFS, nameisr = nomicFS,
      resjsc = todoresjCS, Gssc = GsCS, allGssc = todoGsCS, Jsc = JCS,
      namejsc = nomjCS,
      CAres = CAres))
}
### Fin SimAn: Simultaneous Analysis
################################################################################



################################################################################
#####               Summary of CorrAn                                      #####
################################################################################

CorrAnSummary <- function(input)
{
   CAres <- input
   CAR <- list("CA table" = "CA table",
               "Total inertia" = CAres$totalin,
               "Eigenvalues and percentages of inertia" = CAres$resin,
               "Output for rows" = CAres$resi,
               "Output for columns" = CAres$resj,
               "Output for supplementary rows" = CAres$resisr,
               "Output for supplementary columns" = CAres$resjsc)
   return(CAR = CAR)
}
### Fin CorrAnSummary
################################################################################



################################################################################
#####               Summary of SimAn                                       #####
################################################################################
SimAnSummary <- function(input)
{
   ASG <- input
   ### CA
   CAres <- ASG$CAres
   G <- ASG$G
   CAR <- vector("list", G)
   names(CAR) <- paste("CA", ASG$nameg)
   for(g in 1:G) {
      CAR[[g]] <- list(CAres[[g]]$totalin, CAres[[g]]$resin, CAres[[g]]$resi, 
         CAres[[g]]$resj, CAres[[g]]$resisr, CAres[[g]]$resjsc)
      names(CAR[[g]]) <- c("Total inertia",
      "Eigenvalues and percentages of inertia", "Output for rows",
      "Output for columns", "Output for supplementary rows",
      "Output for supplementary columns")   }
   ### SA
    SAR <- list("Total inertia" = ASG$totalin,
         "Eigenvalues and percentages of inertia" = ASG$resin,
         "Output for rows" = ASG$resi,
         "Output for columns" = ASG$resj,
         "Output for partial rows" = ASG$resig,
         "Projections of tables" = ASG$Fsg,
         "Contributions of tables to SA" = ASG$ctrg,
         "Relation between overall and partial rows" = ASG$riig,
         "Relation between factors of separate CA" = ASG$rCACA,
         "Relation between factors of CA and SA" = ASG$rCASA,
         "Output for supplementary rows" = ASG$resisr,
         "Output for supplementary partial rows" = ASG$resigsr,
         "Output for supplementary columns" = ASG$resjsc)
   ### output CA y SA
   return(list(CAR = CAR, SAR = SAR))
}
### Fin SimAnSummary
################################################################################



################################################################################
#####               CorrAnGraph: Graphs for CorrAn                        ######
################################################################################
CorrAnGraph <- function(input, s1 = 1, s2 = 2, screen = TRUE)
{

   SAPplot <- function(x, y)
   {
         plot(x, y, type = "n", xlab = paste("  "), ylab = paste("  "), col = 1)
   }
   SAPtitle <- function(titulo, titsub, s1, s2, resin)
   {
      title(titulo, col = 1) 
      title("", 
         xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)   "), 
         ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)   "),
         adj = 1, col = 1)
      mtext(titsub, side = 3, line = .4, outer = FALSE)
   }
   SAPplotpuntos <- function(g, puntos, gra, s1, s2, etiq,
      colg = c(rep(2:8, 10)), pchg = c(rep(c(15, 17, 18, 19), 10)))
      {
          points(c(puntos[gra, s1]), c(puntos[gra, s2]), 
             pch = pchg[g], col = colg[g])
          text(c(puntos[gra, s1]), c(puntos[gra, s2]),
             c(paste("  ", etiq[gra], sep = "")), col = colg[g], adj = 0)
      }
   ### Opciones graficos
   par(cex = 1, cex.axis = .6, cex.lab = 1, cex.main = 2, cex.sub = 1)
   colg <- c(rep(2:8, 10))
   pchg <- c(rep(c(15, 17, 18, 19), 10))
   pchgs <- c(rep(c(0, 2, 5, 1), 10))
   ### Input
   G <- 1
   CAres <- vector("list", G)
   names(CAres) <- paste("Table_", 1:G, sep = "")
   CAres[[1]] <- input
   names(CAres[[1]]) <- c("totalin", "eig", "resin", "resi", "resj",
      "resisr", "resjsc",
      "X", "totalk",
      "I", "namei", "fi", "Fs", "d2i",
      "J", "namej", "fj", "Gs", "d2j",
      "Isr", "nameisr", "fisr", "Fssr", "d2isr",
      "Xsr",
      "Jsc", "namejsc", "fjsc", "Gssc", "d2jsc")
   titulinparc <- "CA"
   ### Dimensiones
   I <- CAres[[1]]$I
   ndim <- ncol(CAres[[1]]$Gs)
   IFS <- CAres[[1]]$Isr
   JCS <- CAres[[1]]$Jsc
   J <- 0
   for(g in G:1) {
      J <- CAres[[g]]$J + J}
   ### Nombres
   ###
   nomig <- c(0)
   for(g in G:1) {
      nomig <- c(paste(g, CAres[[g]]$namei, sep = ""), nomig)   }
   nomig <- nomig[1:length(nomig) - 1]
   ###
   nomj <- c(0)
   for(g in G:1) {
      nomj <- c(paste(g, CAres[[g]]$namej, sep = ""), nomj)   }
   nomj <- nomj[1:length(nomj) - 1]
   ###
   nomjg <- c(0)
   for(g in G:1) {
      nomjg <- c(CAres[[g]]$namej, nomjg)   }
   ###
   g <- 1
   nomic <- CAres[[g]]$namei
   nomicFS <- CAres[[g]]$nameisr
   ### Datos para ACS
   Fsgrupos <- array(0, c(I, ndim, G))
   Jg <- CAres[[g]]$J
   Gsgrupos <- array(0, c(Jg, ndim, G))
   for(g in G:1) {
      Fsgrupos[, , g] <- CAres[[g]]$Fs
      Gsgrupos[1:CAres[[g]]$J, , g] <- CAres[[g]]$Gs   }
   if(IFS != 0) {
      FsgruposFS <- array(0, c(IFS, ndim, G))
      for(g in G:1) {
         FsgruposFS[, , g] <- CAres[[g]]$Fssr   }   }
   if(JCS != 0) {
      GsgruposCS <- array(0, c(JCS, ndim, G))
      for(g in G:1) {
         GsgruposCS[1:CAres[[g]]$Jsc, , g] <- CAres[[g]]$Gssc}   }
   #######################   AC activo   #######################################
   for(g in 1:G) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulinparc, sep = "")
      titsub <- "Active elements"
      nomi <- CAres[[g]]$namei
      puntosplot <- rbind(Fsgrupos[, , g], Gsgrupos[1:CAres[[g]]$J, , g])
      resin <- round(CAres[[g]]$resin, 2)
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      SAPplotpuntos(g, puntos = Fsgrupos[, , g], gra = c(1:I), s1, s2,
         etiq = CAres[[g]]$namei, colg = rep(1, G + 1), pchg = rep(24, G + 1))
      SAPplotpuntos(g, puntos = Gsgrupos[1:CAres[[g]]$J, , g],
         gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej)   }
   #######################   AC activo + suplementario   #######################
   if(IFS != 0 || JCS != 0) {
      for(g in 1:G) {
         if(screen == TRUE) {dev.new()}
         titulo <- paste(titulinparc, sep = "")
         titsub <- "Active and supplementary elements"
         resin <- round(CAres[[g]]$resin, 2)
         puntosplot <- rbind(Fsgrupos[, , g], Gsgrupos[1:CAres[[g]]$J, , g])
         nomi <- CAres[[g]]$namei
         ### Filas suplementarias
         if(IFS != 0) {
            puntosplot <- rbind(puntosplot, FsgruposFS[ , , g])   }
         ### Columnas suplementarias
         if(JCS != 0) {
            puntosplot <- rbind(puntosplot, 
               GsgruposCS[1:CAres[[g]]$Jsc, , g])   }
         ### Grafico
         SAPplot(puntosplot[, s1], puntosplot[, s2])
         SAPtitle(titulo, titsub, s1, s2, resin)
         SAPplotpuntos(g, puntos = Fsgrupos[, , g], gra = c(1:I), s1, s2, 
            etiq = CAres[[g]]$namei, colg = rep(1, G + 1), 
            pchg = rep(24, G + 1))
         SAPplotpuntos(g, puntos = Gsgrupos[1:CAres[[g]]$J, , g],
            gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej)
         if(IFS != 0) {
            puntos <- matrix(0, IFS, ndim)
            puntos[1:IFS, ] <- FsgruposFS[1:IFS, , g]
            SAPplotpuntos(g, puntos = puntos, gra = c(1:IFS), s1, s2, 
               etiq = CAres[[g]]$nameisr, colg = rep(1, G + 1),
               pchg = rep(4, G + 1))   }
         if(JCS != 0) {
            puntos <- matrix(0, CAres[[g]]$Jsc, ndim)
            puntos[1:CAres[[g]]$Jsc, ] <- GsgruposCS[1:CAres[[g]]$Jsc, , g]
            SAPplotpuntos(g, puntos = puntos, 
               gra = c(1:CAres[[g]]$Jsc), s1, s2,
               etiq = CAres[[g]]$namejsc, pchg = pchgs)   }   }   }
}
### Fin CorrAnGraph
################################################################################



################################################################################
#####               SimAnGraph: Graphs for SimAn                          ######
################################################################################
SimAnGraph <- function(input, s1 = 1, s2 = 2, screen = TRUE)
{
   SAPplot <- function(x, y)
   {
         plot(x, y, type = "n", xlab = paste("  "), ylab = paste("  "), col = 1)
   }
   SAPtitle <- function(titulo, titsub, s1, s2, resin)
   {
      title(titulo, col = 1)
      title("", 
         xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)   "), 
         ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)   "),
         adj = 1, col = 1)
      mtext(titsub, side = 3, line = 0.4, outer = FALSE)
   }
   SAPplotpuntos <- function(g, puntos, gra, s1, s2, etiq,
      colg = c(rep(2:8, 10)), pchg = c(rep(c(15, 17, 18, 19), 10)))
      {
          points(c(puntos[gra, s1]), c(puntos[gra, s2]), 
             pch = pchg[g], col = colg[g])
          text(c(puntos[gra, s1]), c(puntos[gra, s2]),
             c(paste("  ", etiq[gra], sep = "")), col = colg[g],
             adj = 0)
      }
   ### Opciones graficos
   par(cex = 1, cex.axis = .6, cex.lab = 1, cex.main = 1.2, cex.sub = 1)
   colg <- c(rep(2:8, 10))
   pchg <- c(rep(c(15, 17, 18, 19), 10))
   pchgs <- c(rep(c(0, 2, 5, 1), 10))
   ###
   titulin <- "SA"
   titulinparc <- "CA"
   ASG <- input
   ### Dimensiones
   I <- ASG$I
   Jg <- ASG$maxJg
   G <- ASG$G
   nomgrup <- ASG$nameg
   ndim <- ncol(ASG$Gs)
   IFS <- ASG$Isr
   JCS <- ASG$Jsc
   ### CA
   CAres <- ASG$CAres
   J <- 0
   for(g in G:1) {
      J <- CAres[[g]]$J + J}
   ### Nombres
   ###
   nomig <- c(0)
   for(g in G:1) {
      nomig <- c(paste(g, CAres[[g]]$namei, sep = ""), nomig)   }
   nomig <- nomig[1:length(nomig) - 1]
   ###
   nomj <- c(0)
   for(g in G:1) {
      nomj <- c(paste(g, CAres[[g]]$namej, sep = ""), nomj)   }
   nomj <- nomj[1:length(nomj) - 1]
   ###
   nomjg <- c(0)
   for(g in G:1) {
      nomjg <- c(CAres[[g]]$namej, nomjg)   }
   ###
   nomic <- ASG$namei
   nomicFS <- ASG$nameisr
   ### Datos para CA
   Fsgrupos <- array(0, c(I, ndim, G))
   Gsgrupos <- array(0, c(Jg, ndim, G))
   for(g in G:1) {
      Fsgrupos[, , g] <- CAres[[g]]$Fs
      Gsgrupos[1:CAres[[g]]$J, , g] <- CAres[[g]]$Gs   }
   if(IFS != 0) {
      FsgruposFS <- array(0, c(IFS, ndim, G))
      for(g in G:1) {
         FsgruposFS[, , g] <- CAres[[g]]$Fssr   }   }
   if(JCS != 0) {
      GsgruposCS <- array(0, c(JCS, ndim, G))
      for(g in G:1) {
         GsgruposCS[1:CAres[[g]]$Jsc, , g] <- CAres[[g]]$Gssc}   }
   ### Datos SA
   Fsic <- ASG$Fsi
   Fsig <- ASG$Fsig
   Gs <- ASG$Gs
   todoFs <- ASG$allFs
   todoGs <- ASG$allGs
   ### Resultados SA FS (Filas Suplementarias)
   if(IFS != 0) {
      FsicFS <- ASG$Fsisr
      FsigFS <- ASG$Fsigsr
      todoFsFS <- ASG$allFssr   }
   ### Resultados SA CS (Columnas Suplementarias)
   if(JCS != 0) {
      GsCS <- ASG$Gssc
      todoGsCS <- ASG$allGssc   }
   ### Datos SA para grupos (cuadrado) y correlacios CA y SA (circulo)
   vls <- ASG$resin[, 1]
   ctag <- ASG$ctrg
   Fsg <- ASG$Fsg
   fact <- ASG$rCASA
   ###################   CA activo   ###########################################
   for(g in 1:G) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulinparc, ": Table ", nomgrup[g], sep = "")
      titsub <- "Active elements"
      nomi <- CAres[[g]]$namei
      puntosplot <- rbind(Fsgrupos[, , g], Gsgrupos[1:CAres[[g]]$J, , g])
      resin <- round(CAres[[g]]$resin, 2)
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      SAPplotpuntos(g, puntos = Fsgrupos[, , g], gra = c(1:I), s1, s2,
         etiq = CAres[[g]]$namei, colg = rep(1, G + 1), pchg = rep(24, G + 1))
      SAPplotpuntos(g, puntos = Gsgrupos[1:CAres[[g]]$J, , g],
         gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej)   }
   ###################   CA activo + suplementario   ###########################
   if(IFS != 0 || JCS != 0) {
      for(g in 1:G) {
         if(screen == TRUE) {dev.new()}
         titulo <- paste(titulinparc, ": Table ", nomgrup[g], sep = "")
         titsub <- "Active and supplementary elements"
         resin <- round(CAres[[g]]$resin, 2)
         puntosplot <- rbind(Fsgrupos[, , g], Gsgrupos[1:CAres[[g]]$J, , g])
         nomi <- CAres[[g]]$namei
         ### Filas suplementarias
         if(IFS != 0) {
            puntosplot <- rbind(puntosplot, FsgruposFS[ , , g])   }
         ### Columnas suplementarias
         if(JCS != 0) {
            puntosplot <- rbind(puntosplot,
               GsgruposCS[1:CAres[[g]]$Jsc, , g])   }
         ### Grafico
         SAPplot(puntosplot[, s1], puntosplot[, s2])
         SAPtitle(titulo, titsub, s1, s2, resin)
         SAPplotpuntos(g, puntos = Fsgrupos[, , g], gra = c(1:I), s1, s2, 
            etiq = CAres[[g]]$namei, colg = rep(1, G + 1),
            pchg = rep(24, G + 1))
         SAPplotpuntos(g, puntos = Gsgrupos[1:CAres[[g]]$J, , g],
            gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej)
         if(IFS != 0) {
            puntos <- matrix(0, IFS, ndim)
            puntos[1:IFS, ] <- FsgruposFS[1:IFS, , g]
            SAPplotpuntos(g, puntos = puntos, gra = c(1:IFS), s1, s2, 
               etiq = CAres[[g]]$nameisr, colg = rep(1, G + 1),
               pchg = rep(4, G + 1))   }
         if(JCS != 0) {
            puntos <- matrix(0, CAres[[g]]$Jsc, ndim)
            puntos[1:CAres[[g]]$Jsc, ] <- GsgruposCS[1:CAres[[g]]$Jsc, , g]
            SAPplotpuntos(g, puntos = puntos, gra = c(1:CAres[[g]]$Jsc), 
               s1, s2, etiq = CAres[[g]]$namejsc, pchg = pchgs)   }   }   }
   ###################   SA   ##################################################
   resin <- round(ASG$resin, 2)
   ###################   SA Fsic   #############################################
   if(screen == TRUE) {dev.new()}
   titulo <- paste(titulin, ": Overall rows", sep = "")
   titsub <- "Active elements"
   puntosplot <- todoFs[, , (G + 1)]
   SAPplot(puntosplot[, s1], puntosplot[, s2])
   SAPtitle(titulo, titsub, s1, s2, resin)
   SAPplotpuntos(g = G + 1, puntos = todoFs[, , (G + 1)], gra = c(1: I),
      s1, s2, etiq = nomic)
   ###################   SA Gsj   ##############################################
   if(screen == TRUE) {dev.new()}
   titulo <- paste(titulin, ": Columns", sep = "")
   titsub <- "Active elements"
   puntosplot <- Gs
   SAPplot(puntosplot[, s1], puntosplot[, s2])
   SAPtitle(titulo, titsub, s1, s2, resin)
   for(g in 1:G) {
      SAPplotpuntos(g, puntos = todoGs[1:CAres[[g]]$J, , g],
      gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej)   }
   ###################   SA Fsic + Gsj   #######################################
   if(screen == TRUE) {dev.new()}
   titulo <- paste(titulin, ": Overall rows and columns", sep = "")
   titsub <- "Active elements"
   puntosplot <- rbind(todoFs[, , (G + 1)], Gs)
   SAPplot(puntosplot[, s1], puntosplot[, s2])
   SAPtitle(titulo, titsub, s1, s2, resin)
   SAPplotpuntos(g = G + 1, puntos = todoFs[, , (G + 1)], gra = c(1:I), 
      s1, s2, etiq = nomic)
   for(g in 1:G) {
      SAPplotpuntos(g, puntos = todoGs[1:CAres[[g]]$J, , g],
      gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej)   }
   ###################   SA Fsic + Fsig   ######################################
   if(screen == TRUE) {dev.new()}
   titulo <- paste(titulin, ": Overall and partial rows", sep = "")
   titsub <- "Active elements"
   puntosplot <- todoFs
   SAPplot(puntosplot[, s1, ], puntosplot[, s2, ])
   SAPtitle(titulo, titsub, s1, s2, resin)
   SAPplotpuntos(G + 1, puntos = todoFs[, , (G + 1)], gra = c(1:I), s1, s2, 
      etiq = nomic)
   for(g in 1:G) { 
      SAPplotpuntos(g, puntos = todoFs[, , g], gra = c(1:I), s1, s2, 
         etiq = c(paste(nomgrup[g], nomic, sep = "")))   }
   ###################   fin SA activo   #######################################
   ###################   SA Fsic + FsicFS   ####################################
   if(IFS != 0) {
      if(screen == TRUE) {dev.new()}
      ### Fsic
      titulo <- paste(titulin, ": Overall rows", sep = "")
      titsub <- "Active and supplementary elements"
      g <- G + 1
      puntosplot <- rbind(todoFs[, , (G + 1)], todoFsFS[, , (G + 1)])
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      ### Fsic
      puntos <- todoFs[, , (G + 1)]
      etiq <- c(nomic)
      points(c(puntos[, s1]), c(puntos[, s2]), pch = pchg[g], col = colg[g])
      text(c(puntos[, s1]), c(puntos[, s2]), c(paste("  ", etiq, sep = "")),
         col = colg[g], adj = 0)
      ### FsicFS
      puntos <- matrix(0, IFS, ndim)
      puntos[1:IFS, ] <- todoFsFS[1:IFS, , (G + 1)]
      etiq <- c(nomicFS)
      points(c(puntos[, s1]), c(puntos[, s2]), pch = pchgs[g], col = colg[g])
      text(c(puntos[, s1]), c(puntos[, s2]), c(paste("  ", etiq, sep = "")),
         col = colg[g], adj = 0)   }
   else{   }
   ###################   SA Gsj + GsjCS   ######################################
   if(JCS != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Columns", sep = "")
      titsub <- "Active and supplementary elements"
      puntosplot <- rbind(Gs, GsCS)
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      for(g in 1:G) {
         SAPplotpuntos(g, puntos = todoGs[1:CAres[[g]]$J, , g],
            gra = c(1:CAres[[g]]$J), s1, s2, etiq = CAres[[g]]$namej)   }
      for(g in 1:G) {
         puntos <- matrix(0, JCS, ndim)
         puntos[1:JCS, ] <- todoGsCS[1:JCS, , g]
         SAPplotpuntos(g, puntos = puntos, gra = c(1:JCS), s1, s2, 
            etiq = paste(nomgrup[g], CAres[[g]]$namejsc, sep = ""),
            pchg = pchgs)   }   }
   else{   }
   ###################   SA Fsic + FsicFS + Gsj + GsjCS   ######################
   if(IFS != 0 || JCS != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Overall rows and columns", sep = "")
      titsub <- "Active and supplementary elements"
      puntosplot <- rbind(todoFs[, , (G + 1)], Gs)
      if(IFS != 0) {
         puntosplot <- rbind(puntosplot, todoFsFS[, , (G + 1)])   }
      if(JCS != 0) {
         puntosplot <- rbind(puntosplot, GsCS)   }
      ### Grafico
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      SAPplotpuntos(g = G + 1, puntos = todoFs[, , (G + 1)], gra = c(1:I),
         s1, s2, etiq = nomic)
      for(g in 1:G) {
         SAPplotpuntos(g, puntos = todoGs[1:CAres[[g]]$J, , g],
            gra = c(1:CAres[[g]]$J), s1, s2,
            etiq = CAres[[g]]$namej)   }
      if(IFS != 0) {
         puntos <- matrix(0, IFS, ndim)
         puntos[1:IFS, ] <- todoFsFS[1:IFS, , (G + 1)]
         SAPplotpuntos(g = G + 1, puntos = puntos,
            gra = c(1:IFS), s1, s2, etiq = nomicFS, pchg = pchgs)   }
      if(JCS != 0) {
         for(g in 1:G) {
            puntos <- matrix(0, JCS, ndim)
            puntos[1:JCS, ] <- todoGsCS[1:JCS, , g]
            SAPplotpuntos(g, puntos = puntos, gra = c(1:JCS), s1, s2, 
               etiq = paste(nomgrup[g], CAres[[g]]$namejsc,
               sep = ""), pchg = pchgs)   }   }   }
   else{   }
   ###################   SA Fsic + Fsig + FsicFS + FsigFS   ####################
   if(IFS != 0) {
      if(screen == TRUE) {dev.new()}
      titulo <- paste(titulin, ": Overall and partial rows", sep = "")
      titsub <- "Active and supplementary elements"
      ### En puntosplot deberia ir rbind(todoFs, todoFsFS) que no funciona
      puntosplot <- rbind(Fsic, Fsig, FsicFS, FsigFS)
      SAPplot(puntosplot[, s1], puntosplot[, s2])
      SAPtitle(titulo, titsub, s1, s2, resin)
      ### Activo
      puntos <- todoFs
      etiq <- c(nomic)
      SAPplotpuntos(G + 1, puntos = todoFs[, , (G + 1)], gra = c(1:I), 
         s1, s2, etiq = nomic)
      for(g in 1:G) {
         SAPplotpuntos(g, puntos = todoFs[, , g], gra = c(1:I), s1, s2,
            etiq = c(paste(nomgrup[g], nomic, sep = "")))   }
      ### Suplementario
      puntos <- todoFsFS
      etiq <- c(nomicFS)
      ### FcFS
      points(c(puntos[, s1, (G + 1)]), c(puntos[, s2, (G + 1)]),
         pch = pchgs[G + 1], col = colg[G + 1])
      text(c(puntos[, s1, (G + 1)]), c(puntos[, s2, (G + 1)]),
         c(paste("  ", etiq, sep = "")), col = colg[G + 1], adj = 0)
      ### FgFS
      for(g in 1:G) {
         points(c(puntos[, s1, g]), c(puntos[, s2, g]), pch = pchgs[g],
            col = colg[g])
         text(c(puntos[, s1, g]), c(puntos[, s2, g]), 
            c(paste("  ", nomgrup[g], etiq, sep = "")), 
            col = colg[g], adj = 0)   }   }
   else{   }
   ###################   fin SA activo + suplementario   #######################
   ###################   Fsg y SA + CA (cuadrado y circulo)   ##################
   ###################   SA Fsg (cuadrado)   ###################################
   if(screen == TRUE) {dev.new()}
   titulo <- paste(titulin, ": Tables", sep = "")
   puntosplot <- Fsg
   puntos <- Fsg
   etiq <- nomgrup
   oldpar <- par(pin = c(5, 5), las = 1)
   plot(puntosplot[, s1], puntosplot[, s2], type = "n",
      xlab = paste("  "), ylab = paste("  "),
      col = 1, xlim = c(0, 1), ylim = c(0, 1))
   title(titulo, col = 1) 
   title("", 
      xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)    "), 
      ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)    "),
      adj = 1, col = 1)
   for(g in 1:G) {
      points(c(puntos[g, s1]), c(puntos[g, s2]), pch = pchg[g], col = colg[g])
      text(c(puntos[g, s1]), c(puntos[g, s2]), 
         c(paste("  ", etiq[g], sep = "")), col = colg[g], adj = 0)   }
   ### Volver a pin originales (no grafico cuadrado)
   par(oldpar)
   ###################   SA + CA (circulo)   ###################################
   if(screen == TRUE) {dev.new()}
   titulo <- "Relation between factors of CA and SA"
   puntosplot <- fact
   puntos <- fact
   etiq <- nomgrup
   ### Circulo
   z <- seq(0, 2 * pi, length = 1000)
   x <- sin(z)
   y <- cos(z)
   oldpar <- par(pin = c(5, 5), las = 1)
   plot(x, y, type = "l", axes = FALSE, xlab = "", ylab = "")
   axis(1, pos = 0, labels = FALSE)
   axis(2, pos = 0, labels = FALSE)
   title(titulo, col = 1) 
   title("", 
      xlab = paste("Axis ", s1, "   (", resin[s1, 2], "%)    "), 
      ylab = paste("Axis ", s2, "   (", resin[s2, 2], "%)    "),
      adj = 1, col = 1)
   for(sac in 1:2) {
      for(g in 1:G) {
         points(c(puntos[sac, s1, g]), c(puntos[sac, s2, g]), pch = pchg[g],
            col = colg[g])
         text(c(puntos[sac, s1, g]), c(puntos[sac, s2, g]),
            c(paste("  F", sac, "(", etiq[g], ")", sep = "")), col = colg[g],
            adj = 0)   }   }
   ### Volver a pin originales (no grafico redondo)
   par(oldpar)
   ###################   fin Fsg y SA + CA (cuadrado y circulo)   ##############
}
### fin SimAnGraph
################################################################################



################################################################################
#####               End of SimultAnR (Simultaneous Analysis) Package       #####
################################################################################

