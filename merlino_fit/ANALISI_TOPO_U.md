# Analisi: coordinate topologiche e uso della matrice G

## Perché la base topologica (soprattutto sugli anelli)
- Le combinazioni cicliche (breathing, CVB, torsioni cicliche) sono stabili e continue perché dipendono solo dalla connettività e dall’ordine canonico dell’anello.
- La riduzione via G su una singola geometria può cambiare ordine degli autovettori, generare mixing tra coordinate equivalenti e introdurre discontinuità quando la geometria cambia.

## Topologia ≠ indipendenza lineare completa
La topologia determina l’equivalenza locale, ma non garantisce che tutte le coordinate siano indipendenti.

Esempi:
- Angoli attorno a un centro (CH3): 3 angoli H–C–H → solo 2 indipendenti.
- Dihedrali attorno a un legame: molte combinazioni possibili, ma spesso 1–2 indipendenti.
- OOP: per centri con molti vicini ci sono dipendenze note.
- Anelli fusi: serve la regola M−1 per cluster.

Quindi la strategia migliore è:
1) costruire U “topologico”, deterministico e continuo
2) usare G solo come filtro finale per eliminare ridondanze residue

## Proposta tecnica (concisa)

### 1) Classi atomiche
- Discrete: AtomicSynthons.canonical_signature (Z, NED, D, aromaticità)
- Continua/priority: Zeff o numero quantico principale (solo come tie‑breaker)

### 2) Signature delle primitive
- bond: (kind, class(i), class(j), ring_tag); per X–H aggiungi l’indice dell’atomo pesante (vincolo geminale)
- angle: (kind, class(center), multiset(class(i),class(k)), ring_tag)
- dihedral: (kind, class(j), class(k), class(i), class(l), ring_tag)
- out_of_plane/linear_bend analoghi
- regola: X–H separati da non‑H; geminali sì, vicinali no

### 3) Base locale per gruppo equivalente
- per gruppo di n primitive:
  - v0 = (1,1,…,1)/sqrt(n) (modo simmetrico)
  - v1..v(n−1) da differenze e_i − e_{i+1}, ortonormalizzate
- base deterministica e continua

### 4) Anelli
- breathing su bond dell’anello
- CVB su angoli ciclici
- torsioni cicliche su dihedrali dell’anello
- cluster fusi: applicare M−1 (topologico)

### 5) Filtro G opzionale
- calcolare G = B B^T (o massa pesata se serve) solo alla fine
- SVD/eig e scartare valori < tol → rimuove ridondanze residue
- l’ordine rimane stabile perché U_topo è stabile

## Rischi / casi limite
- Se la firma discreta cambia (es. NED/D cambia per rounding) la base può cambiare. Raro se non ci sono transizioni topologiche.
- Per famiglie quasi degeneri (es. legami quasi lineari) il rango può cambiare → qui il filtro G è necessario.

## Differenza rispetto allo stato attuale
- Oggi: valence_angle_u, ring_*_u, out_of_plane_u usano G‑blocks per scegliere i modi.
- Proposto: costruire prima una base topologica locale, poi eventualmente filtrare con G.
