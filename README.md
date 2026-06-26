# Progetto di Fondamenti di Automatica — Traccia 28

Progetto d'esame per il corso di Fondamenti di Automatica, anno accademico 2023/2024.
Autore: **Bruni Christian** - matricola 240008.

Il lavoro analizza tre sistemi assegnati nella traccia n° 28, svolti tra
Mathematica, MATLAB e Simulink.

## Contenuto

| Cartella | Esercizio | Cosa contiene |
|----------|-----------|---------------|
| `esercizio-a/` | Sistema LTI a tempo continuo (SISO) | `EsercizioA.nb` (notebook originale) e `esercizioA.wl` (sorgente Wolfram estratto e commentato) |
| `esercizio-b/` | Sistema LTI a tempo discreto (SISO) | `EsercizioB.nb` e `esercizioB.wl` |
| `esercizio-c/` | Catena di Markov a tempo discreto | `ParametriMarkov.mlx` (Live Script), `parametriMarkov.m` (codice commentato) e `Simulink2a.slx` (modello della ricorsione) |
| `relazione/` | Relazione | `latex/main.pdf`
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

- **Wolfram Mathematica 14** - apre i notebook `.nb`. I file `.wl` sono lo stesso
  codice in formato testo (leggibile anche senza Mathematica) con commenti aggiuntivi.
- **MATLAB + Simulink** - `ParametriMarkov.mlx` definisce la matrice di transizione e
  lo stato iniziale; `Simulink2a.slx` itera `x(k+1) = A x(k)` fino alla convergenza.
