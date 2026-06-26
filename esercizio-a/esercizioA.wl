(* ::Package:: *)

(* =====================================================================
   Esercizio A - Sistema dinamico LTI a tempo continuo (SISO, proprio)
   Fondamenti di Automatica - Traccia 28 - Bruni Christian, matr. 240008

   Sistema:  x'(t) = A x(t) + B u(t)
             y(t)  = C x(t)              (D = 0  ->  sistema proprio)

   Lo script ripercorre, nell'ordine, i punti della traccia:
     1.  modi naturali (autovalori di A)
     2.  risposta libera per un x0 assegnato
     3.  stati iniziali che "accendono"/"spengono" singoli modi
     4.  funzione di trasferimento, poli e zeri
     5.  risposta al gradino (componente di regime + transitoria)
     6.  risposta al segnale sinusoidale e fattori di distorsione
     7.  modello ARMA equivalente e risposta alla rampa
     8.  x0 che annulla la componente transitoria nella risposta al gradino
     9.  risposta all'ingresso 1(-t)

   Nota: i blocchi 2x2 nella forma canonica derivano dalla coppia di
   autovalori complessi coniugati; A resta diagonalizzabile perché
   per ogni autovalore m.a. = m.g.
   ===================================================================== *)



(* ===== Analisi Sistema LTI-TC ===== *)
ClearAll["Global`*"]

A = {{-179/4, -331/4, 205/8, -409/8}, {183/4, 327/4, -213/8, 409/8},
      {-179/2, -331/2, 205/4, -413/4}, {-179/2, -327/2, 209/4, -409/4}};

B= {{-1/2}, {1/2}, {-1}, {-1}};

C1= {{-10, 0, 0, 5}};

\[Lambda]= Eigenvalues[A]


(* ===== Grafico dei modi naturali: ===== *)
Plot[Exp[\[Lambda][[1]]t],{t,0,4},PlotRange->All]

Plot[Exp[\[Lambda][[2]]t],{t,0,4},PlotRange->All]

Plot[Exp[Re[\[Lambda][[3]]]t]Cos[Im[\[Lambda][[3]]]t],{t,0,4},PlotRange->All]

Plot[Exp[Re[\[Lambda][[4]]]t]Sin[Im[\[Lambda][[4]]]t],{t,0,4},PlotRange->All]

T0=Transpose[Eigenvectors[A]]

(* Matrice di cambio base "ibrida": per la coppia complessa coniugata
   uso parte reale e immaginaria di UN solo autovettore, così T resta
   reale e Lambda assume la forma di Jordan reale (blocco 2x2). *)
T=Transpose[{T0[[All,1]],T0[[All,2]],Re[T0[[All,3]]],Im[T0[[All,4]]]}]


(* ===== Vediamo se la Forma Canonica \[EGrave] diagonale a blocchi, ci aspettiamo la presenza di un blocco di Jordan 2x2 data la presenza di autovalori complessi e coniugati ===== *)
\[CapitalLambda]=Simplify[Inverse[T].A.T]

\[CapitalLambda] // MatrixForm


(* ===== Calcolo della Risposta Libera ===== *)

(* ===== Inizializzo il vettore delle condizioni iniziali ===== *)
Subscript[x, 0]={{-1},{0},{-3},{0}};


(* ===== Proietto lo stato iniziale lungo le colonne della matrice \[OpenCurlyDoubleQuote]ibrida\[CloseCurlyDoubleQuote] T ===== *)
Subscript[z, 0]=Inverse[T].Subscript[x, 0]

\[Sigma]= Re[\[Lambda][[3]]]

\[Omega]= Im[\[Lambda][[3]]]


(* ===== Calcolo la risposta libera sfruttando la decomposizione modale ===== *)
(*Subscript[x, l][t_]:=FullSimplify[T.(\[NoBreak]{{Exp[\[Lambda][[1]]t], 0, 0, 0}, {0, Exp[\[Lambda][[2]]t], 0, 0}, {0, 0, Exp[\[Sigma] t]Cos[\[Omega] t], Exp[\[Sigma] t]Sin[\[Omega] t]}, {0, 0, -Exp[\[Sigma] t]Sin[\[Omega] t], Exp[\[Sigma] t]Cos[\[Omega] t]}}\[NoBreak]).Subscript[z, 0]];*)

Subscript[x, l][t_]:=FullSimplify[T.MatrixExp[\[CapitalLambda] t].Subscript[z, 0]]

Subscript[y, l][t_]:=C1.Subscript[x, l][t];

Plot[{Subscript[x, l][t],Subscript[y, l][t]},{t,0,4},PlotRange->All,PlotLegends->"Expressions"]


(* ===== Vediamo cosa accade quando cambiano le condizioni iniziali ===== *)
T //MatrixForm


(* ===== Tramite l\[CloseCurlyQuote]analisi della combinazione lineare delle colonne associate agli autovalori reali, ottengo ===== *)
(* Idea della selezione dei modi: la risposta libera è combinazione
   lineare delle colonne di T pesate da z0. Scegliendo x0 come somma
   delle SOLE colonne dei modi che voglio "accendere", le restanti
   componenti di z0 si annullano e quei modi non compaiono. *)
Subscript[x, 0]= T[[All,1]](-6663)/(344)+T[[All,2]](6663)/(344)

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

Subscript[x, l][t_]:=Simplify[T.MatrixExp[\[CapitalLambda] t].Subscript[z, 0]]

Subscript[x, l][t]

Plot[{Subscript[x, l][t],Subscript[y, l][t]},{t,0,4},PlotRange->All]


(* ===== Tramite l\[CloseCurlyQuote]analisi della combinazione lineare delle colonne associate agli autovalori complessi e coniugati, ottengo ===== *)
Subscript[x, 0]= T[[All,3]]1+T[[All,4]]1

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

Subscript[x, l][t_]:=Simplify[T.MatrixExp[\[CapitalLambda] t].Subscript[z, 0]]

Plot[{Subscript[x, l][t],Subscript[y, l][t]},{t,0,4},PlotRange->All]


(* ===== \[OpenCurlyDoubleQuote]Accendo\[CloseCurlyDoubleQuote] il primo modo naturale ===== *)
Subscript[x, 0]=Transpose[T][[1]]

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

Plot[{Subscript[x, l][t],Subscript[y, l][t]},{t,0,1},PlotRange->All]


(* ===== \[OpenCurlyDoubleQuote]Accendo\[CloseCurlyDoubleQuote] il secondo modo naturale ===== *)
Subscript[x, 0]=Transpose[T][[2]]

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

Plot[{Subscript[x, l][t],Subscript[y, l][t]},{t,0,2},PlotRange->All]


(* ===== \[OpenCurlyDoubleQuote]Accendo\[CloseCurlyDoubleQuote] il terzo modo naturale ===== *)
Subscript[x, 0]=Transpose[T][[3]]

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

Plot[{Subscript[x, l][t],Subscript[y, l][t]},{t,0,2},PlotRange->All]


(* ===== \[OpenCurlyDoubleQuote]Accendo\[CloseCurlyDoubleQuote] il quarto modo naturale ===== *)
Subscript[x, 0]=Transpose[T][[4]]

Subscript[z, 0]=Inverse[T].Subscript[x, 0]

Plot[{Subscript[x, l][t],Subscript[y, l][t]},{t,0,4},PlotRange->All]


(* ===== Calcolo della Funzione di Trasferimento Cell[ TemplateBox[<|boxesG(s)\[LongEqual]C (s I-A)^(-1) B,errors -> {},input -> $G(s)=C \, \left( s \, I - A\right)^{-1} \, B$,state -> Boxes|>,TeXAssistantTemplate],ExpressionUUID-> 08042abc-0c50-2e4a-95ff-e03c4c1a4e31] ===== *)
G[s_]:=Simplify[(C1.Inverse[s IdentityMatrix[4]-A].B)]

G[s]


(* ===== normalizzo la funzione di trasferimento ===== *)
FdT[s_]:= (Numerator[G[s]]/4)/(s^(4)+14s^(3)+(261)/(4)s^(2)+(253)/(2)s+(357)/(4))

FdT[s]

\[CapitalSigma]= StateSpaceModel[{A,B,C1}];

TransferFunctionModel[\[CapitalSigma]]


(* ===== Calcolo i Poli della Funzione di Trasferimento ===== *)
Solve[Denominator[FdT[s][[1]]]==0,s]

{{s\[Rule]-7},{s\[Rule]-3},{s\[Rule]-2-(\[ImaginaryI])/(2)},{s\[Rule]-2+(\[ImaginaryI])/(2)}}

TransferFunctionPoles[\[CapitalSigma]]


(* ===== Calcolo gli Zeri della Funzione di Trasferimento ===== *)
Solve[Numerator[FdT[s][[1]]]==0,s]

TransferFunctionZeros[\[CapitalSigma]]


(* ===== Risposta al Gradino Unitario ===== *)
Factor[FdT[s]]


(* ===== La risposta al Gradino unitario e\[CloseCurlyQuote] pari a Cell[ TemplateBox[<|boxes(FdT(s))/(s),errors -> {},input -> \frac{FdT(s)}{s},state -> Boxes|>,TeXAssistantTemplate],ExpressionUUID-> a3d65515-f587-1642-9046-f3cdb633eecd] ===== *)
U[s_]:=LaplaceTransform[UnitStep[t],t,s]

Y[s_]:= Factor[FdT[s] U[s]]

Y[s]

Apart[Y[s]]


(* ===== Scrivo in maniera verbosa i fratti semplici ===== *)
(Subscript[C, 1])/(s)+(Subscript[C, 1])/(s+3)+(Subscript[C, 3])/(s+7)+(Subscript[C, 4])/(s+2-(\[ImaginaryI])/(2))+(Subscript[C, 5])/(s+2+(\[ImaginaryI])/(2))


(* ===== Trovo il valore delle C, tramite la formula di Heaviside ===== *)
(* Ogni residuo si ottiene col limite (s - p_i) Y(s) per s->p_i.
   C1 è legato al polo s=0 introdotto dal gradino: è la componente
   di REGIME. I residui sui poli del sistema danno il TRANSITORIO.
   Per la coppia complessa basta C4: C5 = Conjugate[C4]. *)
Subscript[C, 1]=Underscript[\[Limit], s\[Rule]0](s) Y[s][[1,1]]

Subscript[C, 2]=Underscript[\[Limit], s\[Rule]-3](s+3) Y[s][[1,1]]

Subscript[C, 3]=Underscript[\[Limit], s\[Rule]-7](s+7) Y[s][[1,1]]

Subscript[C, 4]=Underscript[\[Limit], s\[Rule]-2+(\[ImaginaryI])/(2)](s+2-(1)/(2)\[ImaginaryI]) *Y[s][[1,1]]

Subscript[C, 5]=Conjugate[Subscript[C, 4]]


(* ===== Scrivo ora la componente di regime della risposta forzata al segnale periodico elementare 1(t) ===== *)
Subscript[y, f][t_]:=Subscript[C, 1] UnitStep[t]+Subscript[C, 2] Exp[-3 t] UnitStep[t]+Subscript[C, 3] Exp[-7 t] UnitStep[t]+2 Exp[-2 t]*ComplexExpand[Re[Subscript[C, 4] Exp[(-\[ImaginaryI])/(2) t]]]\[NonBreakingSpace]UnitStep[t]

Subscript[y, ss][t_]:=Subscript[C, 1] UnitStep[t]

Subscript[y, tr][t_]:=Subscript[C, 2] Exp[-3 t] UnitStep[t]+Subscript[C, 3] Exp[-7 t] UnitStep[t]+2 Exp[-2 t]*ComplexExpand[Re[Subscript[C, 4] Exp[(-\[ImaginaryI])/(2) t]]]\[NonBreakingSpace]UnitStep[t]

Plot[{Subscript[y, f][t],Subscript[y, ss][t],Subscript[y, tr][t]},{t,0,4},PlotRange->All,PlotLegends->"Expressions"]


(* ===== Forma Amplitude Phase ===== *)
(*condisero lap di sin(1)1(t)*)

U[s_]:= LaplaceTransform[Sin[ t]UnitStep[t],t,s]

Y[s_]:= FdT[s][[1,1]]U[s]

Apart[Y[s]]

(Subscript[D, 1])/(s+3)+(Subscript[D, 2])/(s+7)+(Subscript[D, 3])/(s+\[ImaginaryI])+(Subscript[D, 4])/(s-\[ImaginaryI])+(Subscript[D, 5])/(s+2-(\[ImaginaryI])/(2))+(Subscript[D, 6])/(s+2+(\[ImaginaryI])/(2))

Subscript[D, 1]=Underscript[\[Limit], s\[Rule]-3](s+3) Y[s]

Subscript[D, 2]=Underscript[\[Limit], s\[Rule]-7](s+7)Y[s]

Subscript[D, 3]=Underscript[\[Limit], s\[Rule]-\[ImaginaryI]](s+\[ImaginaryI])Y[s]

Subscript[D, 4]=Underscript[\[Limit], s\[Rule]\[ImaginaryI]](s-\[ImaginaryI])Y[s]

Subscript[D, 5]=Underscript[\[Limit], s\[Rule]-2+(\[ImaginaryI])/(2)](s+2-(1)/(2)\[ImaginaryI]) Y[s]

Subscript[D, 6]=Conjugate[Subscript[D, 5]]

Subscript[y, ss][t_]:= 2ComplexExpand[Re[Subscript[D, 4]Exp[I t]]]

Subscript[y, ss][t]

(* determino il fattore di distorsione, fra uscita a regime d'ingresso, dobbiamo scrivere la regime nella forma Ampiezza-fase, la pulsazione \[EGrave] identica nelle due formme, la distorsione non avviene sulla pulsazione*)

Solve[{(2 ((11Cos[t])/(425)+(7 Sin[t])/(425))==X Sin[t+\[Psi]])/.{t->0},(D[2 (+(11Cos[t])/(425)+(7 Sin[t])/(425)),t]==D[X Sin[ t+\[Psi]],t])/.{t->0},X>0},{X,\[Psi]}]

(*il sistema ha un'ampiezza minore dell'unit\[AGrave], il segnale \[EGrave] attenuato, il sistema \[EGrave] in anticipo di psi*)N[(2 Sqrt[(2)/(85)])/(5)]

N[-2 ArcTan[(-170+7 Sqrt[170])/(11 Sqrt[170])]]((180)/(\[Pi]))

Subscript[y, trs][t_]:=Subscript[D, 1]Exp[-3 t]+Subscript[D, 2]Exp[-7 t]+2 Exp[-2 t]*ComplexExpand[Re[Subscript[D, 5] Exp[(-\[ImaginaryI])/(2) t]]]\[NonBreakingSpace]UnitStep[t]

(* siccome la pulsazione non varia fra una forma e l'altra, il fattore di alterazione \[EGrave] dato da ampiezza e fase *)

Plot[{Subscript[y, trs][t],Subscript[y, ss][t],Subscript[y, ss][t]+Subscript[y, trs][t]},{t,0,10},PlotRange->All,PlotLegends->"Expressions"]

X= 0.06135719910778964`;\[Omega]=1;\[Psi]=57.528807709151515`;

Plot[{X Sin[\[Omega] t+ \[Psi]], Subscript[y, ss][t]},{t,0,10},PlotRange->All,PlotLegends->"Expressions"]


(* ===== Risposta al segnale periodico ===== *)
X= 2;

\[Omega]=2;

\[Psi]=0;

U[s_]:= LaplaceTransform[X Sin[ \[Omega] t + \[Psi]] UnitStep[t],t,s]

Y[s_]:= FdT[s][[1]][[1]]U[s]

Y[s]

Apart[Y[s]]

(Subscript[Q, 1])/(s+3)+(Subscript[Q, 2])/(s+7)+(Subscript[Q, 3])/(s+2\[ImaginaryI])+(Subscript[Q, 4])/(s-2\[ImaginaryI])+(Subscript[Q, 5])/(s+2-(\[ImaginaryI])/(2))+(Subscript[Q, 6])/(s+2+(\[ImaginaryI])/(2));

Subscript[Q, 1]=Underscript[\[Limit], s\[Rule]-3](s+3) Y[s];

Subscript[Q, 2]=Underscript[\[Limit], s\[Rule]-7](s+7)Y[s];

Subscript[Q, 3]=Underscript[\[Limit], s\[Rule]-2\[ImaginaryI]](s+2\[ImaginaryI])Y[s];

Subscript[Q, 4]=Underscript[\[Limit], s\[Rule]2\[ImaginaryI]](s-2\[ImaginaryI])Y[s];

Subscript[Q, 5]=Underscript[\[Limit], s\[Rule]-2+(\[ImaginaryI])/(2)](s+2-(1)/(2)\[ImaginaryI]) Y[s];

Subscript[Q, 6]=Conjugate[Subscript[Q, 5]];

Subscript[y, ssin][t_]:= 2ComplexExpand[Re[Subscript[Q, 4]Exp[2*I t]]]

Subscript[y, trsin][t_]:=Subscript[Q, 1]Exp[-3 t]+Subscript[Q, 2]Exp[-7 t]+2 Exp[-2 t]*ComplexExpand[Re[Subscript[Q, 5] Exp[(-\[ImaginaryI])/(2) t]]]\[NonBreakingSpace]UnitStep[t]

Plot[{Subscript[y, trsin][t],Subscript[y, ssin][t],Subscript[y, ssin][t]+Subscript[y, trsin][t]},{t,0,15},PlotRange->All,PlotLegends->"Expressions"]


(* ===== Modello Arma ===== *)
Subscript[x, 0]={{-3},{2},{-3},{3}};

(* Le condizioni iniziali del modello ARMA (y, y', y'', y''' in 0) sono
   legate allo stato iniziale dalla matrice di osservabilità:
   [y(0); y'(0); ...] = Ob . x0, con Ob = [C; C A; C A^2; C A^3]. *)
Ob={C1[[1]],C1[[1]].A,C1[[1]].A.A,C1[[1]].A.A.A}

vettoreCondizioniIniziali= {{Subscript[y, 0]},{Subscript[y, 1]},{Subscript[y, 2]},{Subscript[y, 3]}}

vettoreCondizioniIniziali= Ob.Subscript[x, 0]

FdT[s][[1]][[1]]

modelloARMA= Expand[Denominator[FdT[s][[1]][[1]]]*Subscript[Y, ARMA][s]]== Expand[Numerator[FdT[s][[1]][[1]]]Subscript[U, ARMA][s]]

(*Abbiamo una rappresentazione quivalente I/U, (modello ARMA). E' noto come "Equivalente", in quanto, per qualsiasi ingresso, la risposta forzata sar\[AGrave] cos\[IGrave], ossia equivalente in ambo le rappresentazioni*)

InverseLaplaceTransform[(357 Subscript[Y, ARMA][s])/(4)+(253)/(2) s Subscript[Y, ARMA][s]+(261)/(4) s^(2) Subscript[Y, ARMA][s]+14 s^(3) Subscript[Y, ARMA][s]+s^(4) Subscript[Y, ARMA][s]\[Equal]-5 Subscript[U, ARMA][s]+5 s Subscript[U, ARMA][s],s,t]/.{Subscript[Y, ARMA][s]->LaplaceTransform[y[t],t,s], Subscript[U, ARMA][s]->LaplaceTransform[u[t],t,s]} /.{y[0]->0, y'[0]->0, y''[0]->0, y'''[0]->0, u[0]->0}

EqDifferenziale= (357 y[t])/(4)+(253 y^(\[Prime])[t])/(2)+(261 y^(\[Prime]\[Prime])[t])/(4)+14 y^((3))[t]+y^((4))[t]\[Equal]-5 u[t]+5 u^(\[Prime])[t]

(* porto l'eqauzione del dominio di Laplace*)

eqDiffLap=Simplify[LaplaceTransform[EqDifferenziale,t, s]/.{LaplaceTransform[y[t],t,s]->Subscript[Y, ARMA][s],LaplaceTransform[u[t],t,s]->Subscript[U, ARMA][s],u[0]->0}];

eqDiffLap

Ydif[s_]:=Collect[Solve[eqDiffLap,Subscript[Y, ARMA][s]][[1,1]][[2]], Subscript[U, ARMA][s]]

Ydif[s]

RispostaLiberaARMA= ((506 y[0]+261 s y[0]+56 s^(2) y[0]+4 s^(3) y[0]+261 y^(\[Prime])[0]+56 s y^(\[Prime])[0]+4 s^(2) y^(\[Prime])[0]+56 y^(\[Prime]\[Prime])[0]+4 s y^(\[Prime]\[Prime])[0]+4 y^((3))[0])/(357+506 s+261 s^(2)+56 s^(3)+4 s^(4))) /. {y[0] ->vettoreCondizioniIniziali[[1,1]],y^(\[Prime])[0]->vettoreCondizioniIniziali[[2,1]], y^(\[Prime]\[Prime])[0]->vettoreCondizioniIniziali[[3,1]],y^((3))[0]-> vettoreCondizioniIniziali[[4,1]]}

(* per assicurarci che sia corretto, ottengo la risposta libera anche dalla definizione*)

RispostaLiberaConfronto =Simplify[C1.Inverse[s IdentityMatrix[4] - A].Subscript[x, 0]][[1,1]]


(* ===== Risposta alla rampa unitaria ===== *)
EqDifferenzialeRampa= (357 y[t])/(4)+(253 y^(\[Prime])[t])/(2)+(261 y^(\[Prime]\[Prime])[t])/(4)+14 y^((3))[t]+y^((4))[t]\[Equal]-5 u[t]+5 u^(\[Prime])[t]

(* Uso il modello ARMA per calcolare la risposta forzata*)armaRampa=Simplify[LaplaceTransform[EqDifferenzialeRampa,t,s]/.{LaplaceTransform[y[t],t,s]->Subscript[Y, ARMA][s],LaplaceTransform[u[t],t,s]->Subscript[U, ARMA][s],u[0]->0}]

Yarm[s]=Collect[Solve[armaRampa,Subscript[Y, ARMA][s]][[1,1]][[2]],Subscript[U, ARMA][s]]

(*sostiuisoo Subscript[U, ARMA][s] con il valore della Trasformata di Laplace della rampa e i valori di y'[0] con quelli calcolati in precedenza*)

Subscript[y, RA][t_]:=FullSimplify[InverseLaplaceTransform[Yarm[s]/.{Subscript[U, ARMA][s]->LaplaceTransform[t UnitStep[t],t,s],y[0] ->vettoreCondizioniIniziali[[1,1]],y^(\[Prime])[0]->vettoreCondizioniIniziali[[2,1]], y^(\[Prime]\[Prime])[0]->vettoreCondizioniIniziali[[3,1]],y^((3))[0]-> vettoreCondizioniIniziali[[4,1]]},s,t]]

Subscript[y, RA][t]

(*Separo i valori della risposta di regime e della risposta transitoria*)

Subscript[y, rampaSS][t]=(17260)/(127449)-(20 t)/(357)

Subscript[y, rampaTR][t]=-(150390 \[ExponentialE]^(-7 t))/(4949)+(6791 \[ExponentialE]^(-3 t))/(9)-(8 \[ExponentialE]^(-2 t) (2478522 Cos[(t)/(2)]-5089109 Sin[(t)/(2)]))/(29189)


(* ===== Determinare lo stato iniziale x0 ===== *)
(* il blu \[EGrave] un asintoto boliquo*)xIncognite= {{Subscript[x, 1]},{Subscript[x, 2]},{Subscript[x, 3]},{Subscript[x, 4]}}

libera=Simplify[C1.Inverse[s IdentityMatrix[4]-A].xIncognite][[1]][[1]]

(*numeratore di 3\[Degree] grado, che dipende dai coefficenti del vettore delle x incognite, il denominatore \[EGrave] il polinomio carateristico*)forzata= Simplify[(FdT[s])/(s)][[1]][[1]]

(* calcoliamo la risposta a regime *)regime= (FdT[0])/(s)[[1]][[1]]

(* calcolo transitorio come forzata-regime, se il termine in s al deniminatore scompare, il processo sar\[AGrave] andato a buon fine*)transitorio= Factor[forzata-regime]

(* sommo libera e transitoria, ne estraggo il denominatore, per verificare che sia identicamente nullo*)Numerator[Simplify[Expand[libera+transitorio]]]

(* il polinomio \[EGrave] di 3\[Degree] grado e ho 4 incognite (le x)*)CoefficientList[Numerator[Simplify[Expand[libera+transitorio]]],s]

(* eguaglio a 0 il sistema per trovare i cooefficenti*)Solve[CoefficientList[Numerator[Simplify[Expand[libera+transitorio]]],s]=={0,0,0,0},{Subscript[x, 1],Subscript[x, 2],Subscript[x, 3],Subscript[x, 4]}]

response= FullSimplify[OutputResponse[{\[CapitalSigma],{(2)/(357),(-2)/(357),(4)/(357),0}},1,t]]

FdT[0]


(* ===== Calcolo della risposta per un sistema LTI-TC all\[CloseCurlyQuote]ingresso 1(-t) ===== *)
(*1. Valutazione della risposta per t<0.*)yNeg= FdT[0][[1]][[1]]

(*2. Valutazione della risposta per t>0.*)

Subscript[x, 0]=Inverse[Ob].{{FdT[0][[1]][[1]]},{0},{0},{0}}

yliberas=Simplify[(C1.Inverse[s IdentityMatrix[4]-A].Subscript[x, 0])[[1]][[1]]]

Apart[yliberas]

Subscript[H, 1]=Underscript[\[Limit], s\[Rule]-3](s+3) yliberas

Subscript[H, 2]=Underscript[\[Limit], s\[Rule]-7](s+7)yliberas

Subscript[H, 3]=Underscript[\[Limit], s\[Rule]-2+(\[ImaginaryI])/(2)]((2-(\[ImaginaryI])/(2))+s)yliberas

Subscript[H, 4]= Conjugate[Subscript[H, 3]]

yliberat=Subscript[H, 1]Exp[-3 t]+Subscript[H, 2]Exp[-7 t]+2 Exp[-2 t]*ComplexExpand[Re[Subscript[H, 4] Exp[(-\[ImaginaryI])/(2) t]]]\[NonBreakingSpace]UnitStep[t]

y[t_]:={{\[Piecewise], {{yNeg, t<0}, {yliberat, t>=0}}}}

Plot[y[t],{t,-5,5},PlotRange->All]
