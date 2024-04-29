# PointSpreadFunction calculation

```@docs
PointSpreadFunctions.psf
PointSpreadFunctions.apsf
```

## Excitation Modifications

These excitation modifications will be applied to incorporate the photophysics of the excitation process.
This is particularly relevant for non-linear excitation effects such as exploited in two-photon excitation `modify_square`
or stimulated emission depletion (STED) `modify_STED_exp`, `modify_STED_hyper`.

```@docs
PointSpreadFunctions.modify_ident
PointSpreadFunctions.modify_square
PointSpreadFunctions.modify_STED_exp
PointSpreadFunctions.modify_STED_hyper
```
