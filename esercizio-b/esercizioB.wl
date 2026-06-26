(* ::Package:: *)

(* =====================================================================
   Esercizio B - Sistema dinamico LTI a tempo discreto (SISO, proprio)
   Fondamenti di Automatica - Traccia 28 - Bruni Christian, matr. 240008

   Sistema:  x(k+1) = A x(k) + B u(k)
             y(k)   = C x(k)

   Differenza chiave rispetto all'Esercizio A: qui A NON e'
   diagonalizzabile. L'autovalore -1/5 ha molteplicita' algebrica 3 ma
   geometrica 1, quindi JordanDecomposition restituisce un unico blocco
   di Jordan 3x3. I modi naturali a tempo discreto hanno la forma
   binomiale  C(k,j) lambda^(k-j),  da cui i termini Binomial[k,j].

   Punti della traccia:
     1.  modi naturali
     2.  risposta libera per x0 assegnato
     3.  stati iniziali che attivano/disattivano i modi
     4.  funzione di trasferimento, poli e zeri
     5.  risposta al gradino (Heaviside generalizzata per il polo triplo)
     6.  modello ARMA e risposta all'ingresso a finestra (1 per 0<=k<=10)
     7.  condizioni iniziali ARMA che annullano il transitorio
   ===================================================================== *)



(* ===== Analisi Sistema LTI-TD ===== *)
ClearAll["Global`*"]

A={{-92/125, -107/125, -43/50}, {217/125, -18/125, -7/50}, 
 {-184/125, 36/125, 7/25}};

B={{-1/2}, {1/2}, {-1}};

C1= {{9/4, -7/4, -2}};


(* ===== Calcolo gli autovalori ===== *)
\[Lambda]= Eigenvalues[A]

(* A non e' diagonalizzabile: JordanDecomposition restituisce T (cambio
   base) e Lambda con l'unico blocco di Jordan 3x3 associato a -1/5.
   I tre DiscretePlot sotto mostrano i modi binomiali del blocco:
   lambda^k, C(k,1) lambda^(k-1), C(k,2) lambda^(k-2). *)
{T,\[CapitalLambda]}=JordanDecomposition[A]

\[CapitalLambda] // MatrixForm

MatrixPower[\[CapitalLambda],k] // MatrixForm

Subscript[x, 0]={{0},{0},{-3}}

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

DiscretePlot[(-(1)/(5))^(k),{k,0,4},PlotRange->All]

DiscretePlot[Binomial[k,1](-(1)/(5))^(k-1),{k,0,7},PlotRange->All]

DiscretePlot[Binomial[k,2](-(1)/(5))^(k-2),{k,0,7},PlotRange->All]


(* ===== Calcolo la risposta libera ===== *)
Subscript[x, l][k_]:=Simplify[T.MatrixPower[\[CapitalLambda],k].Subscript[z, 0]]

Subscript[x, l][k]

Subscript[y, l][k_]:=Simplify[C1.Subscript[x, l][k]]

Subscript[y, l][k]

DiscretePlot[Evaluate[Subscript[y, l][k]],{k,0,10},PlotRange->All]

DiscretePlot[{Evaluate[Subscript[y, l][k]], Evaluate[Subscript[x, l][k]]},{k,0,8},PlotRange->All, PlotLegends-> {"componente 1","componente 2","componente 3", "componente 4"}]


(* ===== \[OpenCurlyDoubleQuote]Accendiamo\[CloseCurlyDoubleQuote] i modi naturali ===== *)
Subscript[x, 0]=Transpose[T][[1]]

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

DiscretePlot[{Subscript[y, l][t],Subscript[x, l][t]},{t,0,4},PlotRange->All, PlotLegends->{"(*SubscriptBox[(y), (l)])[t]","(*SubscriptBox[(x), (l)])[t]"}]

Subscript[x, 0]=Transpose[T][[2]]

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

DiscretePlot[{Subscript[y, l][t],Subscript[x, l][t]},{t,0,4},PlotRange->All, PlotLegends->{"(*SubscriptBox[(y), (l)])[t]","(*SubscriptBox[(x), (l)])[t]"}]

Subscript[x, 0]=Transpose[T][[3]]

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

DiscretePlot[{Subscript[y, l][t],Subscript[x, l][t]},{t,0,4},PlotRange->All, PlotLegends->{"(*SubscriptBox[(y), (l)])[t]","(*SubscriptBox[(x), (l)])[t]"}]


(* ===== Primi due modi: ===== *)
T // MatrixForm

Subscript[x, 0]= T[[All,2]]1+T[[All,3]]1

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

DiscretePlot[{Subscript[y, l][t],Subscript[x, l][t]},{t,0,4},PlotRange->All, PlotLegends->{"(*SubscriptBox[(y), (l)])[t]","(*SubscriptBox[(x), (l)])[t]"}]

Subscript[x, 0]= T[[All,1]]1+T[[All,3]]1

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

DiscretePlot[{Subscript[y, l][t],Subscript[x, l][t]},{t,0,4},PlotRange->All, PlotLegends->{"(*SubscriptBox[(y), (l)])[t]","(*SubscriptBox[(x), (l)])[t]"}]


(* ===== calcolo la FdT a partire dalla definizione Cell[ TemplateBox[<|boxesG(z)\[LongEqual]C(z I-A)^(-1)B,errors -> {},input -> G(z)=C \left( z \, I - A \right)^{-1} B,state -> Boxes|>,TeXAssistantTemplate],ExpressionUUID-> c334d1bd-ce97-0f49-9e39-0abb9ce6088b] ===== *)
G[z_]:=Simplify[C1.Inverse[z IdentityMatrix[3]-A].B][[1]][[1]]

G[z]

Solve[Denominator[G[z]]==0,z]

Solve[Numerator[G[z]]==0,z]

\[CapitalSigma]=StateSpaceModel[{A,B,C1},SamplingPeriod->1]

TransferFunctionModel[\[CapitalSigma]]

TransferFunctionPoles[\[CapitalSigma]]

TransferFunctionZeros[\[CapitalSigma]]

{{{-(1)/(8)}}}


(* ===== Risposta al gradino unitario ===== *)
Subscript[Y, f][z_]:=G[z]((z)/(z-1))

Factor[(Subscript[Y, f][z])/(z)]

Subscript[C, 1]((1)/(z-1))+Subscript[C, 21]((1)/(z+(1)/(5)))+Subscript[C, 22]((1)/((z+(1)/(5))^(2)))+Subscript[C, 23]((1)/((z+(1)/(5))^(3)))

Subscript[C, 1]=Underscript[\[Limit], z\[Rule]1](z-1)((Subscript[Y, f][z])/(z))


(* ===== Heaviside generalizzata ===== *)
Subscript[C, 23]=Underscript[\[Limit], z\[Rule](-1)/(5)](z+(1)/(5))^(3)((Subscript[Y, f][z])/(z))

Subscript[C, 22]=Underscript[\[Limit], z\[Rule](-1)/(5)]D[(z+(1)/(5))^(3)((Subscript[Y, f][z])/(z)),z]

Subscript[C, 21]=((1)/(2))Underscript[\[Limit], z\[Rule](-1)/(5)]D[D[(z+(1)/(5))^(3)((Subscript[Y, f][z])/(z)),z],z]

Subscript[C, 1]((z)/(z-1))+Subscript[C, 21]((z)/(z+(1)/(5)))+Subscript[C, 22]((z)/((z+(1)/(5))^(2)))+Subscript[C, 23]((z)/((z+(1)/(5))^(3)))


(* ===== Risposta forzata nel dominio del tempo ===== *)
Subscript[y, f][k_]:=Subscript[C, 1]UnitStep[k]+Subscript[C, 21](-(1)/(5))^(k)UnitStep[k]+Subscript[C, 22]Binomial[k,1](-(1)/(5))^(k-1)UnitStep[k]+Subscript[C, 23]Binomial[k,2](-(1)/(5))^(k-2)UnitStep[k]

Subscript[y, ss][k_]:=Subscript[C, 1]UnitStep[k]

Subscript[y, tr][k_]:=Subscript[C, 21](-(1)/(5))^(k)UnitStep[k]+Subscript[C, 22]Binomial[k,1](-(1)/(5))^(k-1)UnitStep[k]+Subscript[C, 23]Binomial[k,2](-(1)/(5))^(k-2)UnitStep[k]

DiscretePlot[{Evaluate[Subscript[y, f][k]],Evaluate[Subscript[y, ss][k]],Evaluate[Subscript[y, tr][k]]},{k,0,8},PlotRange->All,PlotLegends->{"(*SubscriptBox[(y), (f)])[k]","(*SubscriptBox[(y), (ss)])[k]","(*SubscriptBox[(y), (tr)])[k]"}]


(* ===== Modello ARMA ===== *)
Subscript[x, 0]={{0},{-2},{-3}}

Subscript[y, 0]={{Subscript[y, 1]},{Subscript[y, 2]},{Subscript[y, 3]}}

Subscript[y, 1]= FullSimplify[C1.Subscript[x, 0]][[1]][[1]]

Subscript[y, 2]=FullSimplify[C1.A.Subscript[x, 0]][[1]][[1]]

Subscript[y, 3]=FullSimplify[C1.A.A.Subscript[x, 0]][[1]][[1]]

Subscript[y, 0]

(* Modello ARMA: eguaglio den(G) Y(z) = num(G) U(z), forma implicita I/U. *)
Intedita= Expand[Denominator[G[z]]*Y[z]]== Expand[Numerator[G[z]]U[z]]

(* Antitrasformo z^n Y(z) -> y(k+n) per ottenere l'equazione alle
   differenze (ricorrenza) nel dominio del tempo discreto. *)
Ricorrenza=Intedita/.{Y[z]\[Rule]y[k],z Y[z]\[Rule]y[k+1],z^(2) Y[z]\[Rule]y[k+2],z^(3)Y[z]\[Rule]y[k+3],U[z]\[Rule]u[k],z U[z]\[Rule]u[k+1]}

Ricorrenza=(ZTransform[Ricorrenza,k,z])/.{ZTransform[y[k],k,z]\[Rule]Y[z],ZTransform[u[k],k,z]\[Rule]U[z],u[0]\[Rule]0};

(*Adesso vado a risovlere l'equazione differenziale, ottenendo la risposta del sistema nel dominio della trasformata Z*)

Risposta=Collect[Solve[Ricorrenza, Y[z]][[1,1]][[2]], U[z]]

RispostaLibera = (1)/(4 (1+5 z)^(3))5 (12 z y[0]+60 z^(2) y[0]+100 z^(3) y[0]+60 z y[1]+100 z^(2) y[1]+100 z y[2])

(*Determino la risposta libera,sfruttando i valori appena calcolati*)

RispostaLibera = ((1)/(4 (1+5 z)^(3))5 (12 z y[0]+60 z^(2) y[0]+100 z^(3) y[0]+60 z y[1]+100 z^(2) y[1]+100 z y[2]))/.{y[0] -> Subscript[y, 0][[1,1]], y[1]-> Subscript[y, 0][[2,1]], y[2]->Subscript[y, 0][[3,1]]}

RispostaLiberaConfronto =Simplify[ C1.Inverse[z IdentityMatrix[3]-A].Subscript[x, 0]][[1,1]]

(*Come si pu\[OGrave] vedere i due valori coincidono, allora ho trovato un modello ARMA equivalente. l'equivalenza tra le due rappresentazioni conferma la correttezza del calcolo delle condizioni iniziali, e di conseguenza della risposta libera del sistema. *)


(* ===== Risposta all\[CloseCurlyQuote]ingresso u(k) ===== *)
(* definisco il mio ingresso*)

Subscript[U, pw][z]= Underoverscript[\[Sum], k=0, 10]z^(-k)

Subscript[Y, pw][z]= Risposta /. {U[z]->Subscript[U, pw][z], y[0]->Subscript[y, 1],y[1]->Subscript[y, 2],y[2]->Subscript[y, 3]}

(*effettuo l'antitrasformata*)

Subscript[y, pw][k_]:= InverseZTransform[Subscript[Y, pw][z],z,k]

Subscript[y, pw][k]

DiscretePlot[Evaluate[Subscript[y, pw][k]],{k,-3,20},PlotRange->All]

(*risposta forzata*)

Subscript[y, fpw][k_]:= InverseZTransform[(1)/(4 (1+5 z)^(3))5 (1+(1)/(z^(10))+(1)/(z^(9))+(1)/(z^(8))+(1)/(z^(7))+(1)/(z^(6))+(1)/(z^(5))+(1)/(z^(4))+(1)/(z^(3))+(1)/(z^(2))+(1)/(z)) (25+200 z),z,k]

DiscretePlot[Evaluate[Subscript[y, fpw][k]],{k,-3,20},PlotRange->All]


(* ===== Condizioni Iniziali ARMA ===== *)
y0={{Subscript[y, 11]},{Subscript[y, 22]},{Subscript[y, 33]}}

x0 = {{Subscript[x, 11]},{Subscript[x, 22]},{Subscript[x, 33]}}

U[z]= ZTransform[UnitStep[k],k,z]

RispostaCI= ((5 (25+200 z) U[z])/(4 (1+5 z)^(3))+(1)/(4 (1+5 z)^(3))5 (12 z y[0]+60 z^(2) y[0]+100 z^(3) y[0]+60 z y[1]+100 z^(2) y[1]+100 z y[2]))/.{y[0]->Subscript[y, 11],y[1]->Subscript[y, 22],y[2]->Subscript[y, 33]}

(* divido le componenti *)Rlibera= (5 (12 z Subscript[y, 11]+60 z^(2) Subscript[y, 11]+100 z^(3) Subscript[y, 11]+60 z Subscript[y, 22]+100 z^(2) Subscript[y, 22]+100 z Subscript[y, 33]))/(4 (1+5 z)^(3))Rforzata= (5 z (25+200 z))/(4 (-1+z) (1+5 z)^(3))Rregime=  (G[1]z)/(z-1)Rtransitoria = Factor[Rforzata-Rregime]

(* azzero libera+transitoria*)

Numerator[Simplify[Expand[Rlibera+Rtransitoria]]]

CoefficientList[Numerator[Simplify[Expand[Rlibera+Rtransitoria]]],z]

Solve[CoefficientList[Numerator[Simplify[Expand[Rlibera+Rtransitoria]]],z]=={0,0,0,0},{Subscript[y, 11],Subscript[y, 22],Subscript[y, 33]}]

Ob = {{C1}[[1,1]],{C1.A}[[1,1]],{C1.A.A}[[1,1]]}

y0= {{(125)/(96)},{(125)/(96)},{-(67)/(96)}}

x0= Inverse[Ob].y0

