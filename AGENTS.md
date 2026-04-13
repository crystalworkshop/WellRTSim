# Project Instructions

## MATLAB Code Style

- Do not write production-defensive MATLAB code by default.
- Avoid excessive guard clauses, fallback branches, compatibility scaffolding, and silent returns.
- Prefer direct code that assumes valid project inputs and fails loudly when assumptions are violated.
- Add checks only for real invariants that are known to break in this codebase.
- Optimize for readability and scientific workflow, not generic library robustness.
- When a simpler implementation is possible, prefer it over a more defensive one.
