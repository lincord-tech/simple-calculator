# Section 23. Explicit Hamiltonian Density in the Gauge-Fixed Theory

## 23.1 Setup

We work in the gauge-fixed theory after harmonic (de Donder) gauge,
with the two physical sectors completely decoupled (as established in Sect. 18–19):

**Tensor sector** (TT graviton):
$$\mathcal{L}_{TT} = \frac{M_{Pl}^2}{8}\, \dot{h}^{TT}_{\mu\nu}\dot{h}^{TT\,\mu\nu} - \frac{M_{Pl}^2}{8}(\nabla h^{TT}_{\mu\nu})(\nabla h^{TT\,\mu\nu})$$

**Scalar sector** (physical chameleon scalar Φ, canonically normalized after diagonalization):
$$\mathcal{L}_\Phi = \frac{1}{2}\dot{\Phi}^2 - \frac{1}{2}(\nabla\Phi)^2 - \frac{1}{2}m^2(\rho)\Phi^2$$

where $m^2(\rho) = \sigma\rho - \alpha \geq 0$ (positive in dense matter, tuned near zero in vacuum).

These are two decoupled sectors. We compute $\mathcal{H}$ for each via the Legendre transform.

---

## 23.2 Tensor sector Hamiltonian

Define conjugate momentum:
$$\pi^{\mu\nu}_{TT} = \frac{\partial \mathcal{L}_{TT}}{\partial \dot{h}^{TT}_{\mu\nu}} = \frac{M_{Pl}^2}{4}\,\dot{h}^{TT\,\mu\nu}$$

Legendre transform:
$$\mathcal{H}_{TT} = \pi^{\mu\nu}_{TT}\,\dot{h}^{TT}_{\mu\nu} - \mathcal{L}_{TT}$$

Substituting $\dot{h}^{TT}_{\mu\nu} = \frac{4}{M_{Pl}^2}\pi_{\mu\nu}^{TT}$:

$$\mathcal{H}_{TT} = \frac{4}{M_{Pl}^2}\pi^{\mu\nu}_{TT}\pi_{\mu\nu}^{TT} - \frac{M_{Pl}^2}{8}\left(\frac{4}{M_{Pl}^2}\right)^2\pi^{\mu\nu}_{TT}\pi_{\mu\nu}^{TT} + \frac{M_{Pl}^2}{8}(\nabla h^{TT}_{\mu\nu})^2$$

$$\mathcal{H}_{TT} = \frac{4}{M_{Pl}^2}\pi^{\mu\nu}\pi_{\mu\nu} - \frac{2}{M_{Pl}^2}\pi^{\mu\nu}\pi_{\mu\nu} + \frac{M_{Pl}^2}{8}(\nabla h^{TT}_{\mu\nu})^2$$

$$\boxed{\mathcal{H}_{TT} = \frac{2}{M_{Pl}^2}\,\pi^{\mu\nu}_{TT}\pi_{\mu\nu}^{TT} + \frac{M_{Pl}^2}{8}(\nabla h^{TT}_{\mu\nu})(\nabla h^{TT\,\mu\nu})}$$

Both terms manifestly positive (squares). $\mathcal{H}_{TT} \geq 0$. ✓

---

## 23.3 Scalar sector Hamiltonian

Conjugate momentum:
$$\Pi = \frac{\partial \mathcal{L}_\Phi}{\partial \dot{\Phi}} = \dot{\Phi}$$

Legendre transform:
$$\mathcal{H}_\Phi = \Pi\,\dot{\Phi} - \mathcal{L}_\Phi = \Pi^2 - \left[\frac{1}{2}\Pi^2 - \frac{1}{2}(\nabla\Phi)^2 - \frac{1}{2}m^2\Phi^2\right]$$

$$\boxed{\mathcal{H}_\Phi = \frac{1}{2}\Pi^2 + \frac{1}{2}(\nabla\Phi)^2 + \frac{1}{2}m^2(\rho)\Phi^2}$$

All three terms are manifestly non-negative:
- $\frac{1}{2}\Pi^2 \geq 0$ — kinetic energy (positive definite)
- $\frac{1}{2}(\nabla\Phi)^2 \geq 0$ — gradient energy (positive definite)
- $\frac{1}{2}m^2(\rho)\Phi^2 \geq 0$ — mass term, since $m^2(\rho) = \sigma\rho - \alpha \geq 0$ in dense matter

$\mathcal{H}_\Phi \geq 0$ without any conditions on $\Phi$ or $\Pi$. ✓

---

## 23.4 Total Hamiltonian density

$$\boxed{\mathcal{H} = \underbrace{\frac{2}{M_{Pl}^2}\pi^{\mu\nu}_{TT}\pi_{\mu\nu}^{TT} + \frac{M_{Pl}^2}{8}(\nabla h^{TT}_{\mu\nu})^2}_{\mathcal{H}_{TT}\,\geq\, 0} + \underbrace{\frac{1}{2}\Pi^2 + \frac{1}{2}(\nabla\Phi)^2 + \frac{1}{2}m^2(\rho)\Phi^2}_{\mathcal{H}_\Phi\,\geq\, 0}}$$

**$\mathcal{H} \geq 0$ — explicitly, term by term, without reference to eigenvalues.**

The theory is bounded below. No ghost, no tachyon (for $m^2 \geq 0$), no runaway.

---

## 23.5 Note on the vacuum ($\rho \approx \rho_{vac}$, $m \approx 0$)

In the cosmological vacuum $m^2(\rho_{vac}) = \sigma\rho_{vac} - \alpha \approx 0$ by tuning (Sect. 11).
The mass term vanishes and $\mathcal{H}_\Phi = \frac{1}{2}\Pi^2 + \frac{1}{2}(\nabla\Phi)^2 \geq 0$.
The scalar is massless (or ultra-light) in vacuum — cosmologically active — but still bounded below.

In matter ($\rho \gg \alpha/\sigma$), $m^2 \approx \sigma\rho > 0$: mass term adds to positivity.
The Hamiltonian is strictly positive and grows with density. This is the kinetic screening: the scalar is energetically suppressed inside dense objects.

---

## 23.6 Summary table

| Sector | $\mathcal{H} \geq 0$? | Condition |
|--------|----------------------|-----------|
| TT graviton | ✓ manifest | none |
| Scalar (matter) | ✓ manifest | $m^2 = \sigma\rho - \alpha > 0$ |
| Scalar (vacuum) | ✓ manifest | $m^2 \approx 0$, still bounded below |
| Ghost | ✗ absent | all kinetic coefficients positive |
| Tachyon | ✗ absent | $m^2 \geq 0$ enforced by parameter tuning |

**Stability at the quadratic level is established completely and explicitly.**
