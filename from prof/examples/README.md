# Examples

Quick overview of the runnable demos in this folder. Run them from the repo root after `python -m pip install -e .`.

- `example_cylinder_flow.py` — Analytical cylinder solution: Cp distribution on the surface plus flow-field streamlines/heatmap, optional circulation.
- `example_panel_method.py` — NACA 4412 panel-method solution: Cp distribution and flow field computed from the solver.
- `example_circle_cp_comparison.py` — Convergence study: panel-method Cp on a circle vs analytical, includes max-error vs panel-count log–log fit.
- `example_panel_geometry_plot.py` — Panel geometry for NACA 4412: segments, tangents, and normals visualized for a coarse discretization.
- `example_singularity_flow_fields.py` — Flow-field gallery: freestream at two angles, point source/vortex at the origin, and uniform panel source/vortex strengths on a short line.
- `example_utils_contours.py` — Utility contours: circle and NACA 2412 generation visualized with key reference points.
- `example_utils_naca4_parametric.py` — Parametric NACA 4-series family: vary camber, camber location, or thickness and plot resulting shapes.
