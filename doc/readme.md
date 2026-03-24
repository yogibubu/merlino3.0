This document describes the internal architecture of Merlino 3.0.

# Merlino 3.0

Merlino 3.0 è un framework minimale e modulare per la preparazione,
visualizzazione e manipolazione di strutture molecolari,
basato su un **unico canale di comunicazione**: il file `xyzin`.

Il progetto è pensato per sviluppo incrementale, robustezza
e separazione rigorosa delle responsabilità.

---

## 📁 Struttura delle directory

merlino3.0/
├─ app.py
├─ manager.py
├─ working/
│ └─ xyzin
├─ gui/
├─ geometry/
├─ merlino_fit/
├─ projects/
│ ├─ se_library/
│ ├─ pcs2_library/
│ └─ hpcs2_library/
└─ doc/


`projects/` qui non è un contenitore di "progetti" nel senso organizzativo.
È una libreria dati locale usata dai workflow di similarità e frammentazione.

Le aree operative principali sono:
- `gui/` per l'interfaccia e i workflow utente
- `geometry/` per pipeline rotazionali, vibrazionali e termochimiche
- `merlino_fit/` per similarity, fragment pipeline e delta correction
- `working/` come workspace runtime

---

## 🔑 Concetti chiave

### `working/`
- Workspace **effimero**
- Esiste sempre
- Viene pulito **a fine sessione**
- In modalità `DEBUG` non viene pulito

Nota pratica: dal punto di vista architetturale `working/` è runtime-only.
Se qualche file dentro `working/` risulta tracciato da git, quello è un problema
di stato del repository, non una scelta di design.

### `xyzin`
- **Single source of truth**
- Contiene:
  - blocco XYZ
  - sezioni (`#BASIC`, `#SMILES`, …)
- Nessun altro canale di comunicazione è permesso

---

## 🔁 Ciclo di vita

1. `manager.py` crea `working/xyzin` con H₂ e `#BASIC`
2. La GUI modifica `xyzin` tramite i readers
3. Il viewer legge **solo** `xyzin`
4. Alla chiusura:
   - `working/` viene pulita
   - se `DEBUG=True`, nulla viene rimosso

Nota: la GUI scrive un log unico in `working/gui.log`. Se un eventuale
`working/fchkin` ha un numero di atomi diverso da `xyzin`, la parte
vibrazionale viene saltata con un warning.

---

## 🧪 Workflow rovib (CLI Gaussian)

Quando usi `cli_gaussian.py`, oltre al workflow standard viene eseguito
anche il passo `rovib` se disponibile. Se manca `dos_vib.dat`, il DOS
vibrazionale viene generato automaticamente con i seguenti default:

- vmax = 6 (tutti i modi)
- emax = 8000 cm^-1
- bin = 50 cm^-1
- ncap = 10

Questo consente di ottenere in automatico:
- `dos_vib.dat`
- `dos_rovib.dat`
- `rovib_qt.dat` (Q(T) rovibrazionale)

## 🧩 Workflow attivi aggiuntivi

### Fragment pipeline e delta correction
- Fragment-level similarity su librerie `SE`, `PCS2`, `HPCS2`
- Preparazione bundle di correzione locale
- Applicazione di correzioni geometriche locali su frammenti

### Bridge DeltaVib da matrice vibro-rotazionale `alpha`
- Lettura di `alpha` da Gaussian log
- Somma selettiva per modo
- inversione opzionale del segno per frequenze immaginarie
- scrittura di `ΔVib` nella sezione rotazionale di `xyzin`

Nota di architettura:
- questa non è più la linea scientifica principale per il problema vibro-rotazionale
- il lavoro principale vive nella nuova linea `CeDiTT + alpha_resonances`
- in `Merlino 3.0` il blocco `DeltaVib/alpha` va mantenuto come ponte di compatibilità applicativa

---

## 🧩 Readers

Tutti i readers rispettano il contratto:

NPUT → (file, SMILES, ecc.)
OUTPUT → working/xyzin
RETURN → None


Reader attivi:
- `smiles_reader.py` → aggiorna XYZ + `#SMILES`
- `xyz_reader.py` → sostituisce solo XYZ

La gestione del formato `xyzin` è centralizzata in:
- `xyzin_utils.py`

---

## 🖼️ Visualizzazione

- Viewer 2D basato su RDKit (`viewer2d.py`)
- Fallback automatico al logo se la molecola non è visualizzabile

---

## 🧪 Avogadro

- Avogadro (v1) è usato come **default viewer** (macOS via `open`)
- Avogadro2 è disponibile su richiesta
- Avogadro lavora **sempre su una copia temporanea**
- L’import in `xyzin` è **esplicito e controllato**

---

## 🐞 DEBUG

Per preservare `working/` dopo la chiusura:

```python
run(debug=True)

Utile per:

ispezione di xyzin

debug dei readers

test manuali

🎯 Filosofia

minimalismo

zero stato nascosto

nessuna logica duplicata

ogni modulo fa una cosa sola

semplicità > feature premature

Merlino 3.0 è una base solida, non un prodotto finito.
