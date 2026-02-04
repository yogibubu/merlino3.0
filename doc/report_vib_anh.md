Report sintetico: modulo DOS/Q(T) anarmonico con ibrido diretto + WL

Obiettivo
- Calcolare DOS e/o Q(T) da ω e χ (VPT2/HDCPT2) con stabilizzazione alle alte eccitazioni usando una saturazione continua dei numeri quantici.
- Strategia ibrida: somma diretta per bassi numeri quantici + WL su finestre per energia alta, con stitching e normalizzazione su overlap.

Modello energetico
- Energia vibrazionale (cm⁻¹) con saturazione sui numeri quantici:
  - n_i = v_i + 1/2
  - n_i → ñ_i = n_c,i * erf(n_i / n_c,i)
  - E = Σ_i ω_i ñ_i + Σ_{i≤j} χ_ij ñ_i ñ_j
- Per n_i piccoli: recupera i termini lineari e quadratici (VPT2).
- Per n_i grandi: ñ_i satura → energia non “si ripiega” verso valori non fisici.

Metodo DOS/Q(T)
- Low‑energy: enumerazione deterministica (vmax e/o cutoff energetico).
- High‑energy: Wang–Landau in finestre con overlap, facilmente parallelizzabile.
- Stitching: si sceglie un “anchor” (tipicamente la parte diretta) e si riallineano le finestre WL sui bin comuni.
- Q(T): derivata da DOS con Q(β) = Σ g(E) e^{-βE}.

Implementazione
- File: vib_anh.py legge xyzin (freq_cm1, chi_cm1, ZPEh_KJmol, ZPEa_KJmol), calcola DOS e Q(T).
- Comandi:
  - direct: DOS da somma diretta (low‑v)
  - wl: DOS WL su una finestra
  - stitch: normalizza e unisce finestre
  - q: Q(T) da DOS
- Saturazione impostabile con --ncap; vmax per controllare lo spazio degli stati.

Uso tipico (ibrido)
- Direct:
  python vib_anh.py xyzin direct --vmax 6 --emax 8000 --bin 50 --ncap 10 --out dos_low.dat
- WL su finestre (parallelo):
  python vib_anh.py xyzin wl --vmax 12 --emin 8000 --emax 14000 --bin 50 --ncap 10 --out dos_w1.dat
- Stitch ancorato al direct:
  python vib_anh.py xyzin stitch --windows dos_low.dat,dos_w1.dat --anchor 0 --out dos_full.dat
- Q(T):
  python vib_anh.py xyzin q --dos dos_full.dat --t 298.15

Rischi e limiti
- La saturazione erf è un’ipotesi di regolarizzazione; va tarata (n_c,i) e validata.
- WL richiede vmax: senza limiti lo spazio è infinito.
- La normalizzazione tra finestre dipende dalla qualità dell’overlap.

Criteri pratici per scegliere n_c,i
- Se hai D0 (o un limite energetico fisico): scegli n_c,i in modo che la saturazione del modo i contribuisca in modo coerente alla soglia energetica (es. frazione di D0 proporzionale a ω_i).
- Se non hai D0: usa il turning point della quadratica unaria come stima:
  n_tp,i ≈ ω_i / (2 |χ_ii|) − 1/2
  poi scegli n_c,i ≈ 0.8–0.9 * n_tp,i per mantenere monotonia e prevenire curvature non fisiche.
- In alternativa: calibra n_c,i per far coincidere un piccolo set di livelli alti con dati di riferimento (se disponibili).

Test numerico minimo consigliato (H2O come prova rapida)
- Calcola DOS con quadratico puro vs erf‑saturato.
- Confronta Q(T) a 298 K e ad alte T (es. 1000 K) per misurare l’impatto della regolarizzazione.
- Controlla la monotonia dell’energia e l’assenza di livelli negativi ad alti v.

Originalità vs banalità
- “Banalità” nel senso tecnico: usare VPT2/Dunham per DOS e Q(T) e WL su finestre con stitching è standard.
- “Originalità” nel tuo progetto: l’uso di una saturazione morbida sui numeri quantici (erf) per prevenire la non‑fisicità della forma quadratica, mantenendo i primi termini VPT2, è una scelta non standard ma plausibile. Non è rivoluzionaria, però è elegante e difendibile se motivata e validata.

Schema algoritmico (ibrido)
1) Leggi ω, χ, ZPE da xyzin.
2) Definisci n_c,i e vmax.
3) Calcola DOS low‑E con enumerazione diretta fino a E_split (o vmax_low).
4) Suddividi [E_split, E_max] in finestre con overlap.
5) Esegui WL per ciascuna finestra in parallelo.
6) Stitching: ancora la DOS WL alla DOS diretta usando l’overlap.
7) Calcola Q(T) via somma su g(E) se serve.

Scelta numero finestre WL
- Regola pratica: ogni finestra dovrebbe contenere abbastanza bin per avere istogramma “piatto”.
- Se ΔE_tot è il range, usa finestre da 2–4k cm⁻¹ con overlap 10–20%.
- Troppo piccole: overhead di stitching; troppo grandi: WL lento e meno uniforme.
- Verifica che l’overlap includa regioni con g(E) non troppo rumorosa.

Parametri WL consigliati (indicativi)
- Piccolo (≤6 modi): bin 25–50 cm⁻¹, finestra 2–3k cm⁻¹, overlap 15–25%, sweeps 2–5k, logf_final 1e‑6.
- Medio (7–15 modi): bin 50–100 cm⁻¹, finestra 3–5k cm⁻¹, overlap 10–20%, sweeps 5–10k, logf_final 1e‑6.
- Grande (≥16 modi): bin 100–200 cm⁻¹, finestra 5–8k cm⁻¹, overlap 10–15%, sweeps 1–2e4, logf_final 1e‑7.

Stima pratica di E_max
- Se hai D0: usa E_max ≈ D0 − E0 (energia disponibile sopra lo ZPE).
- Se non hai D0: usa un limite interno basato sulla monotonia delle energie:
  - Per ciascun modo i: n_tp,i ≈ ω_i / (2|χ_ii|) − 1/2.
  - Definisci n_cap,i ≲ 0.8–0.9 n_tp,i e valuta E al set (v_i ≈ n_cap,i).
  - Usa E_max ≈ E(v_cap) − E0 come limite conservativo.
- Per DOS cinetico: può essere utile estendere E_max del 10–20% per garantire convergenza delle integrali in TST/RRKM, poi verificare la stabilità dei risultati.
