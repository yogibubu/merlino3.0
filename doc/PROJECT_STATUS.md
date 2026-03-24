PROJECT STATUS — Merlino 3.0
============================

Stato attuale (checkpoint)
--------------------------
Data: 2026-03-24

Obiettivo attuale
-----------------
Rimettere sotto controllo il repository mentre prosegue lo sviluppo attivo su tre linee:
- GUI e pipeline Gaussian/DOS/Q(T)
- fragment pipeline + delta correction PCS2/HPCS2
- compatibilità futura con la nuova linea vibro-rotazionale `CeDiTT + alpha_resonances`

Situazione corrente
-------------------
Il repository non è in stato "sporco da cleanup", ma in stato "sporco da sviluppo attivo":
- vari file tracciati sono modificati lungo le tre linee sopra
- alcuni file sorgente nuovi erano presenti ma non ancora tracciati
- `projects/` non è un registry di progetti: è una libreria dati (`se_library`, `pcs2_library`, `hpcs2_library`)
- `working/` dovrebbe restare effimera, ma `working/xyzin` è ancora tracciato e va trattato in una pulizia git dedicata, non mescolato allo sviluppo

Gerarchia scientifica corretta:
- la linea principale su `DeltaVib`, `alpha` e trattamento vibro-rotazionale non è più `merlino3.0`
- la sede principale è ora la linea `CeDiTT + alpha_resonances`
- in `merlino3.0` questa parte va mantenuta solo come ponte operativo e futura integrazione compatibile

Implementato e verificato
-------------------------
GUI e pipeline:
- campo separato “Gaussian properties” (log/.out) senza sovrascrivere geometria
- DOS/Q(T) settings dialog con validazione live + reset default
- status panel con link ai file + warning + path
- `gui.log` con timestamp
- `summary.txt` generato automaticamente
- toolbar con accessi rapidi a working, symmetry, similarity e fragment pipeline

DOS/Q(T):
- calcolo DOS vib (harm/anharm) da freq/chi
- convoluzione con rotazionale (rovib)
- TS: esclusione modo immaginario e output `N(E)` cumulativo
- warning TS se `Emin > 0`

Fragment delta correction:
- presente un modulo dedicato `survibfit.fragment_delta_correction`
- presenti i workflow `prepare`, `prepare-hpcs2`, `apply`
- presenti test mirati sul blocco delta correction

Bridge DeltaVib / alpha:
- presente una finestra GUI dedicata `gui/deltavib_alpha_dialog.py`
- integrata in `gui/input_panel.py`
- supporta inversione del segno `alpha` per frequenze immaginarie e scrittura di `ΔVib` in `xyzin`
- questo blocco va inteso come bridge temporaneo verso la linea principale `CeDiTT + alpha_resonances`

Verifiche minime fatte
---------------------
- `python -m pytest tests/test_fragment_delta_correction.py` -> `4 passed`
- `python -m py_compile` sui moduli GUI collegati al nuovo dialog `ΔVib/alpha` -> ok

File attivi principali
----------------------
- `gui/main_window.py`
- `gui/input_panel.py`
- `gui/fragment_pipeline_window.py`
- `gui/deltavib_alpha_dialog.py`
- `gui/gaussian.py`
- `geometry/rotational_pipeline.py`
- `merlino_fit/survibfit/fragment_pipeline.py`
- `merlino_fit/survibfit/fragment_delta_correction.py`
- `merlino_fit/tests/test_fragment_delta_correction.py`

Documenti di triage
------------------
- `doc/REPO_TRIAGE_2026-03-24.md`
- `doc/WORKTREE_TRIAGE_2026-03-24.md`
- `doc/ROVIB_COMPATIBILITY_INTERFACE_2026-03-24.md`

Come riprendere
---------------
1) Se serve, attivare l'ambiente: `source ~/.bashrc && merlino-set`
2) Avviare GUI: `python3 /Users/vincenzobarone/merlino3.0/manager.py`
3) Per il blocco delta correction:
   - generare `fragment_pipeline.json`
   - usare `prepare-hpcs2` o `prepare`
   - applicare la correzione con `apply`
4) Per il bridge `ΔVib/alpha`:
   - caricare un Gaussian log con matrice vibro-rotazionale `alpha`
   - usare il dialog dedicato dal pannello input solo come integrazione locale provvisoria
   - considerare `CeDiTT + alpha_resonances` come sorgente scientifica primaria del metodo

Da sistemare dopo
-----------------
- decidere come gestire git per `working/xyzin`
- classificare e committare in modo logico i file attivi oggi modificati
- fare una pulizia separata degli artefatti/documenti non attivi
- definire l'interfaccia finale di compatibilità tra Merlino e la nuova linea `CeDiTT + alpha_resonances`
