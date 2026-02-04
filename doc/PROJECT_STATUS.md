PROJECT STATUS — Merlino 3.0
============================

Stato attuale (checkpoint)
--------------------------
Data: 2026-01-28

Obiettivo
---------
Riorganizzare la GUI e la pipeline per:
- distinguere geometria vs proprietà (Gaussian log separato)
- calcolare DOS/Q(T) vib e rovib con Emin/Emax e T singola
- riconoscere automaticamente minimo vs TS (frequenza negativa)
- generare anche numero di stati per TS

Implementato
------------
GUI e pipeline:
- Campo separato “Gaussian properties” (log/.out) senza sovrascrivere geometria
- FCHK/FCH esclusi da questa pipeline (tooltip e blocco)
- DOS/Q(T) settings dialog con validazione live + reset default
- Persistenza settings in working/gui_settings.json
- Status panel con link ai file + warning + path
- gui.log con timestamp
- summary.txt generato automaticamente
- Toolbar: separatori + “Open working folder”

DOS/Q(T):
- Calcolo DOS vib (harm/anharm) da freq/chi
- Convoluzione con rotazionale (rovib)
- TS: esclusione modo immaginario e output N(E) cumulativo
  - n_vib_ts.dat
  - n_rovib_ts.dat
- Warning TS se Emin > 0

Documentazione:
- doc/gui.md aggiornato con nuovo flusso + output
- doc/rovib_dos_manual.md con nota metodologica e test minimi

File principali modificati
--------------------------
- gui/main_window.py
- gui/input_panel.py
- gui/gaussian.py
- gui/project_manager.py
- cli_modules.py
- cli_gaussian.py
- geometry/vib_anh.py
- geometry/rovib_pipeline.py
- doc/gui.md
- doc/rovib_dos_manual.md
- doc/readme.md
- working/vib_anh.py (wrapper)

Come riprendere
--------------
1) Se serve, attivare venv: `source ~/.bashrc && merlino-set`
2) Avviare GUI: `python3 /Users/vincenzobarone/merlino3.0/manager.py`
3) Per test TS:
   - usare un xyzin con frequenza negativa
   - verificare n_vib_ts.dat e n_rovib_ts.dat in working/

Da fare (opzionale)
-------------------
- Test TS reale con file di esempio
- RRKM/TST (N^‡(E), k(E), k(T))
- Hindered rotors per modi < 100 cm^-1
