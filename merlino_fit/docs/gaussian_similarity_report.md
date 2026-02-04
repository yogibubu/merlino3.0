# Report: Similarita Gaussiana in Spazio dei Sintoni

## Obiettivo
Definire una metrica continua di similarita molecolare coerente con il layer descrittivo dei sintoni, evitando classi atomiche rigide.

## Descrittori usati
Per ogni atomo viene costruito il vettore:
- `charge` (`q_i`)
- `covalency` (`C_i`)
- `delocalization` (`D_i`)
- `strain` (`S_i`)
- `zeff` (`Z_i^eff`)

Questi 5 descrittori formano la matrice `X` (`natoms x 5`) della molecola.

In aggiunta, per ogni ciclo rilevato viene costruito un vettore ring-aware:
- `ring_size`
- `ring_planarity`
- `ring_radial_pos` (posizione del centro ciclo rispetto al centro molecolare)
- `ring_fused_degree` (connettivita tra cicli)
- `ring_nn_centroid_dist` (distanza al centroide del ciclo piu vicino)
- `ring_radius` (raggio RMS del ciclo)

Questi descrittori formano la matrice `R` (`nrings x 6`).

## Modello Gaussiano molecolare
Ogni molecola e rappresentata con un modello gaussiano:
- media `mu = mean(X)`
- covarianza `Sigma = cov(X) + lambda*I`

Modalita disponibili:
- `full`: covarianza completa
- `diag`: solo diagonale della covarianza

Regolarizzazione:
- `lambda > 0` (default `5e-2`) per stabilita numerica e invertibilita.

## Distanza e similarita
Date due molecole `A` e `B`, con modelli gaussiani `(mu_A, Sigma_A)` e `(mu_B, Sigma_B)`:
- `Sigma = 0.5*(Sigma_A + Sigma_B)`
- `Delta = mu_A - mu_B`

Distanza di Bhattacharyya:

`D_B = 1/8 * Delta^T * pinv(Sigma) * Delta + 1/2 * [logdet(Sigma) - 1/2*(logdet(Sigma_A)+logdet(Sigma_B))]`

Similarita finale:

`Sim = exp(-D_B)` con `Sim in (0, 1]`.

Per i cicli si applica lo stesso schema gaussiano ottenendo `Sim_ring`.
Il punteggio finale usato per il confronto e:

`Sim_comb = (1 - w_ring) * Sim + w_ring * Sim_ring`

dove `w_ring` e configurabile (`--ring-weight`, default `0.25`).
Se entrambe le molecole non hanno cicli, il contributo ring viene disattivato automaticamente.

Interpretazione:
- `Sim ~ 1`: ambienti locali molto simili
- `Sim` bassa: distribuzioni sinthoniche diverse

## Standardizzazione
Di default viene fatta standardizzazione globale prima del fit:
- si concatena `X_A` e `X_B` (o query+library)
- si normalizza ogni feature con media/deviazione standard globali

Questo evita che una feature domini la distanza solo per scala numerica.

## Uso operativo
Confronto a coppie:

```bash
python -m survibfit.synthon_similarity --xyz-a mol1.xyz --xyz-b mol2.xyz
```

Ranking contro libreria:

```bash
python -m survibfit.synthon_similarity \
  --query-xyz query.xyz \
  --library-dir ./library_xyz \
  --library-glob "*.xyz" \
  --top-k 10
```

Opzioni chiave:
- `--covariance-mode full|diag`
- `--regularization <float>`
- `--no-standardize`
- `--ring-weight <float>`
- `--no-ring-comparison`
- `--json-out report.json`

## Validazioni presenti
Test automatici gia presenti:
- similarita alta per molecola identica (`Sim > 0.999`)
- similarita minore per molecole diverse
- ranking libreria ordinato per `Sim` decrescente

File test: `merlino_fit/tests/test_synthon_similarity.py`.

## Limiti e note
- La qualita dipende dalla robustezza dei descrittori topologici/sinthonici.
- Molecole con pochi atomi hanno stima covarianza piu fragile (mitigata da `lambda`).
- Se una covarianza non e definita positiva, il codice segnala errore (controllo su `slogdet`).

## Conclusione
La similarita gaussiana fornisce una misura continua, interpretabile e coerente con il framework dei sintoni. E adatta sia al confronto pairwise sia al retrieval su librerie molecolari.
