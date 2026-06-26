# Progetto di Fondamenti di Automatica — Traccia 28

Progetto d'esame per il corso di Fondamenti di Automatica.
Autore: **Bruni Christian** — matricola 240008.

Il lavoro analizza tre sistemi assegnati nella traccia n° 28, svolti tra
Mathematica, MATLAB e Simulink, con una relazione che commenta tutti i passaggi.

## Contenuto

| Cartella | Esercizio | Cosa contiene |
|----------|-----------|---------------|
| `esercizio-a/` | Sistema LTI a tempo continuo (SISO, proprio) | `EsercizioA.nb` (notebook originale) e `esercizioA.wl` (sorgente Wolfram estratto e commentato) |
| `esercizio-b/` | Sistema LTI a tempo discreto (SISO, proprio) | `EsercizioB.nb` e `esercizioB.wl` |
| `esercizio-c/` | Catena di Markov a tempo discreto | `ParametriMarkov.mlx` (Live Script), `parametriMarkov.m` (codice commentato) e `Simulink2a.slx` (modello della ricorsione) |
| `relazione/` | Relazione | `bruni_christian_240008.pdf` (versione consegnata) e `latex/` (riscrittura in LaTeX) |
| `traccia/` | Testo del progetto | `traccia_23_24_28.pdf` |

## Gli esercizi in breve

- **A — tempo continuo.** Sistema del 4° ordine: modi naturali, risposta libera,
  attivazione selettiva dei modi tramite le condizioni iniziali, funzione di
  trasferimento (poli e zeri), risposta al gradino e a un segnale sinusoidale,
  modello ARMA equivalente, risposta alla rampa e al segnale `1(-t)`.
- **B — tempo discreto.** Sistema del 3° ordine con matrice **non** diagonalizzabile
  (blocco di Jordan 3×3, modi binomiali): stessi punti del caso continuo, con
  Heaviside generalizzata per il polo triplo e risposta a un ingresso a finestra.
- **C — catena di Markov.** Grafo di transizione, stato stazionario calcolato sia
  numericamente (ricorsione simulata in Simulink) sia in forma chiusa, e individuazione
  di uno spanning tree.

## Software

- **Wolfram Mathematica 14** — apre i notebook `.nb`. I file `.wl` sono lo stesso
  codice in formato testo (leggibile anche senza Mathematica) con commenti aggiuntivi.
- **MATLAB + Simulink** — `ParametriMarkov.mlx` definisce la matrice di transizione e
  lo stato iniziale; `Simulink2a.slx` itera `x(k+1) = A x(k)` fino alla convergenza.

## Relazione in LaTeX

La cartella `relazione/latex/` contiene la relazione riscritta in LaTeX: la matematica
è ri-tipografata nativamente, mentre i grafici e gli output di Mathematica/Simulink
sono le immagini estratte dalla relazione originale (cartella `img/`).

Per compilare il PDF:

```bash
cd relazione/latex
latexmk -pdf main.tex
```

## Note

I sorgenti `.wl` e `.m` sono stati ricavati dai notebook originali per renderli
leggibili e diffabili direttamente su GitHub; i file `.nb`/`.mlx`/`.slx` originali
restano nel repository. La cache di compilazione di Simulink (`slprj/`, `*.slxc`) è
esclusa tramite `.gitignore`.
