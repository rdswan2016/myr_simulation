---
name: feedback-numerical-debugging
description: "User wants numerical/ODE bugs root-caused and confirmed via diagnostics before proceeding, not patched blindly"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 99165efa-e454-4d46-833e-8127ba9c7a8c
---

When a mechanistic ODE model produces implausible results (e.g. mass-balance violations, values exceeding physical
bounds), stop and diagnose the actual numerical mechanism (check solver step points, raw un-interpolated state
trajectories, not just the dense-output/interpolated view) before applying a fix. Two fix attempts were tried and
rejected in sequence on the DS-Pool reactor model (hard clamp → absorbing negative state; sign-preserving
continuation → runaway blowup from large rate constants) before the smooth-floor regularization worked — each
failure was diagnosed with concrete printed evidence (solver step sizes collapsing, exact moles trajectories)
rather than guessed at.

**Why:** the user paused mid-fix (rejected a tool call, said "check again") when a third numerical patch was about
to be applied without first surfacing the tradeoff. They wanted to be asked whether to keep debugging the
fractional-order kinetics numerically or fall back to a simpler, numerically well-behaved kinetic model (n=m=1
baseline), rather than have me keep unilaterally patching. See [[project_dspool_reactor_model]].

**How to apply:** for reactor/ODE modeling work in this project (and similar first-principles simulation work),
when a numerical fix isn't obviously correct on the first attempt, pause and offer the user a choice (keep
debugging vs. simplify the model) rather than silently iterating through multiple patches. Always validate fixes
with a concrete check (e.g. exact mass balance, no blow-up) before treating them as done.
