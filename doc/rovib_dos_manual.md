Manuale uso: DOS rovibrazionale e Q(T) da DOS
=============================================

Scopo
-----
Questo modulo calcola il DOS rovibrazionale e Q(T) rovibrazionale a partire da:
1) DOS vibrazionale (da vib_anh)
2) livelli rotazionali (modello quantistico coerente con thermo_rot.py)
3) convoluzione in log-space

File coinvolti
--------------
- Input: `xyzin` nel working directory
  - usa #BASIC (T_K)
  - usa #ROTATIONAL (rotor_type, A/B/C in MHz, Symm. Number)
- Input: DOS vibrazionale (file testo con colonne `E_cm1  log_g`)
- Output:
  - `dos_rovib.dat` (DOS rovibrazionale)
  - `rovib_qt.dat` (Q(T) rovibrazionale, una riga)
  - opzionale: `dos_rot.dat` (DOS rotazionale)

Formato DOS
-----------
File testo, una riga per bin:
```
# format: E_cm1 log_g
123.456789  5.4321000000
...
```

Passo 1: DOS vibrazionale (da vib_anh)
--------------------------------------
Esempio (somma diretta):
```
python merlino3.0/working/vib_anh.py xyzin direct \
  --vmax 6 --emax 8000 --bin 50 --ncap 10 --out dos_vib.dat
```

Esempio (WL + stitching, opzionale):
```
python merlino3.0/working/vib_anh.py xyzin wl \
  --vmax 12 --emin 8000 --emax 14000 --bin 50 --ncap 10 --out dos_w1.dat

python merlino3.0/working/vib_anh.py xyzin stitch \
  --windows dos_vib.dat,dos_w1.dat --anchor 0 --out dos_vib_full.dat
```

Passo 2: Pipeline principale (CLI)
----------------------------------
Esecuzione minima (usa default nel working dir):
```
python merlino3.0/cli_modules.py <working_dir> --rovib
```

Default:
- vib DOS: `dos_vib.dat`
- rovib DOS: `dos_rovib.dat`
- Q(T) rovib: `rovib_qt.dat`
- emax_rot: usa `Emax_vib` se non specifichi `--rovib-emax-rot` o `--rovib-jmax`

Override principali:
```
python merlino3.0/cli_modules.py <working_dir> --rovib \
  --rovib-vib-dos dos_vib_full.dat \
  --rovib-emax-rot 2000 \
  --rovib-rot-out dos_rot.dat \
  --rovib-out dos_rovib.dat \
  --rovib-qout rovib_qt.dat \
  --rovib-t 298.15
```

Nota: se vuoi forzare costanti o rotor_type:
```
--rovib-rotor-type linear|spherical|symmetric
--rovib-A <MHz> --rovib-B <MHz> --rovib-C <MHz>
--rovib-sigma <int>
```

Workflow automatico con cli_gaussian.py
--------------------------------------
Se usi `cli_gaussian.py`, il DOS vibrazionale viene generato automaticamente
quando manca `dos_vib.dat`, con i seguenti parametri di default:
- vmax = 6 (tutti i modi)
- emax = 8000 cm^-1
- bin = 50 cm^-1
- ncap = 10

Questo permette di eseguire sempre il passo `rovib` senza intervento manuale.

Output Q(T)
-----------
File `rovib_qt.dat`:
```
# T_K Q_rovib
298.150000  1.234567890123e+08
```

Uso diretto (senza pipeline)
----------------------------
Se preferisci usare direttamente il comando di convoluzione:
```
python merlino3.0/working/vib_anh.py xyzin convrot \
  --vib-dos dos_vib.dat --emax-rot 2000 --out dos_rovib.dat
```
e poi:
```
python merlino3.0/working/vib_anh.py xyzin q --dos dos_rovib.dat --t 298.15
```

Note e limiti
-------------
- Per rotori non lineari si usa il modello symmetric top con Beff = max(B,C).
- La precisione alle alte energie dipende da:
  - scelta di ncap (saturazione vib)
  - emax_rot o jmax
  - bin del DOS
- Se il DOS vibrazionale non copre abbastanza energia, Q(T) puo essere sottostimato.

Nota metodologica (breve)
-------------------------
- Il DOS vibrazionale viene calcolato in base a frequenze (armonico) o a ω+χ
  (anarmonico con saturazione erf sui numeri quantici).
- Il DOS rotazionale usa un modello quantistico di symmetric top:
  per rotori non lineari si usa Beff = max(B,C).
- La convoluzione vib‑rot fornisce un DOS rovibrazionale statistico
  (non include accoppiamenti Coriolis o effetti rovibrazionali fini).

Validazione minima consigliata (2–3 test rapidi)
------------------------------------------------
1) H2O (test rapido)
   - Calcola DOS vib armonico vs anarmonico (se χ presente).
   - Confronta Q_vib a 298 K e 1000 K; controlla la sensibilita di Emax.
   - Verifica che Q_rovib > Q_vib e che cresca con T.
   - Pass/fail minimo:
     - Q_vib(1000 K) > Q_vib(298 K)
     - Q_rovib(298 K) > Q_vib(298 K)
     - Q_rovib(1000 K) > Q_rovib(298 K)

2) Molecola lineare (es. CO2 o HCN)
   - Verifica che il rotore lineare dia Q_rovib coerente con Q_rot classico
     a temperature medio‑alte (ordine di grandezza).
   - Pass/fail minimo:
     - Q_rot (da thermo_rot) e Q_rot (da DOS rot) differiscono < 1 ordine di grandezza.

3) Molecola asimmetrica semplice (es. H2O o H2CO)
   - Confronta Q_rot (thermo_rot) con Q_rot estratto dal DOS rotazionale
     integrato; differenze moderate sono attese per il modello Beff.
   - Pass/fail minimo:
     - Rapporto Q_rot(DOS)/Q_rot(thermo) in [0.3, 3] per T medio‑alte.

Tabella rapida test → input → comando
-------------------------------------
Test              | File richiesti             | Comando essenziale
H2O               | xyzin con freq (+chi)      | vib_anh.py direct + convrot + q
Lineare (CO2/HCN) | xyzin con #ROTATIONAL      | vib_anh.py convrot + q
Asimmetrica       | xyzin con #ROTATIONAL      | vib_anh.py convrot + q

Esempio H2O (DOS vib + rovib, T=298.15):
```
python merlino3.0/working/vib_anh.py xyzin direct \
  --vmax 6 --emax 8000 --bin 50 --ncap 10 --out dos_vib.dat
python merlino3.0/working/vib_anh.py xyzin convrot \
  --vib-dos dos_vib.dat --emax-rot 2000 --out dos_rovib.dat
python merlino3.0/working/vib_anh.py xyzin q --dos dos_vib.dat --t 298.15
python merlino3.0/working/vib_anh.py xyzin q --dos dos_rovib.dat --t 298.15
```
