# Rovibrational Compatibility Interface

Date: 2026-03-24

## Purpose

Define the downstream interface that `merlino3.0` should expect from the new
primary rovibrational line:

- `CeDiTT`
- `alpha_resonances`

Merlino is not the primary scientific home of `DeltaVib` / `alpha`.
Merlino should only consume a stable compatibility payload.

## Architectural Rule

Primary methodology lives outside this repository:

- theory and filtering logic: `CeDiTT + alpha_resonances`
- downstream consumer and GUI integration: `merlino3.0`

Therefore:

- Merlino must not hard-code scientific assumptions beyond the minimum
  compatibility contract
- Merlino should accept externally prepared `DeltaVib` values and metadata
- unknown compatibility keys in `#ROTATIONAL` must remain harmless

## Minimum Contract

The minimum stable contract for Merlino is the presence of these lines in
`#ROTATIONAL`:

- `DVibA_MHz= ...`
- `DVibB_MHz= ...`
- `DVibC_MHz= ...`

These are the only values that Merlino strictly needs for downstream use.

## Recommended Extended Contract

When the values come from `CeDiTT + alpha_resonances`, the producer should also
be allowed to write optional provenance keys in `#ROTATIONAL`, for example:

- `DVibSource= CeDiTT+alpha_resonances`
- `DVibMethod= alpha_sum`
- `DVibImagHandling= invert`
- `DVibExternalReport= /abs/path/to/report.json`
- `DVibVersion= 1`

Merlino should ignore unknown keys unless a later feature explicitly uses them.

## File-Level Interchange

Preferred interchange for now:

1. external tool computes `DeltaVib`
2. external tool or bridge writes the three `DVib*_MHz` values into `xyzin`
3. Merlino reads them from `#ROTATIONAL`

Optional richer interchange:

- a sidecar JSON report may exist externally
- Merlino may store only a path/reference to that report in `#ROTATIONAL`
- Merlino should not need the full scientific payload to remain functional

## Bridge Behavior Inside Merlino

Current local `DeltaVib/alpha` code in Merlino should be treated as:

- temporary bridge
- local convenience
- fallback/manual integration path

It should not define the canonical scientific meaning of the rovibrational data.

## Implementation Consequence

The `gui.xyzin_utils` layer should support generic writes into `#ROTATIONAL`
without rejecting unknown keys. This allows future importers to attach
provenance from `CeDiTT + alpha_resonances` while keeping current Merlino
parsers stable.

## Next Step

If a direct bridge is implemented later, it should expose a function or import
workflow equivalent to:

- input: external `DeltaVib` payload
- output: update of `#ROTATIONAL` in `xyzin`

with the three mandatory `DVib*_MHz` keys and optional metadata.
