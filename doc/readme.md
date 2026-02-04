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
│ ├─ main_window.py
│ ├─ input_panel.py
│ ├─ viewer_panel.py
│ ├─ readers.py
│ ├─ smiles_reader.py
│ ├─ xyz_reader.py
│ ├─ xyzin_utils.py
│ ├─ viewer2d.py
│ └─ merlino_logo.png


Ogni directory è **autoconsistente** e può essere sviluppata e testata
indipendentemente dalle altre.

---

## 🔑 Concetti chiave

### `working/`
- Workspace **effimero**
- Esiste sempre
- Viene pulito **a fine sessione**
- In modalità `DEBUG` non viene pulito

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
