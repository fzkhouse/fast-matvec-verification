# ==================================================
# w is from discrete Gaussian distribution over Z^d with center 0
# C: matrix sampling from the distribution of GHL21 over {-1,0,1}
# Simulate the probability distribution of ||Cw|| using the chi‑distribution,
# and compare the PDF/CDF/CCDF of the true distribution with the theoretical chi‑distribution.
# ==================================================

import numpy as np
from scipy.stats import chi
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


def discrete_gaussian_rejection_sampling(sigma, d):

    w = np.zeros(d, dtype=int)
    for i in range(d):
        while True:
            k = np.random.randint(-int(3 * sigma), int(3 * sigma) + 1)
            rho = np.exp(-np.pi * (k ** 2) / (sigma ** 2))
            if np.random.uniform(0, 1) < rho:
                w[i] = k
                break
    return w

def estimate_Cw_norm_dist(w, lam):
    # a : the "Scaling Factor" of the theoretical chi‑square distribution
    a = np.linalg.norm(w) / np.sqrt(2)
    k = 2 * lam
    return a, k

def generate_C(lam, d):
    """ C ∈ {-1,0,1}^{2λ×d}, satisfying Pr[0]=1/2, Pr[±1]=1/4"""
    return np.random.choice([-1, 0, 1], size=(2 * lam, d), p=[1 / 4, 1 / 2, 1 / 4])


np.random.seed(42)
d = 256            # vector dimension
lam = 128          # security parameter
sigma_gauss = 30   # Gaussian parameter (corresponding to \sigma in the paper)


w = discrete_gaussian_rejection_sampling(sigma_gauss, d)


a, k = estimate_Cw_norm_dist(w, lam)


n_trials = 1000000
norm_vals = np.zeros(n_trials)

print("Monte‑Carlo simulation is in progress, please wait...")
for i in range(n_trials):
    C = generate_C(lam, d)
    Z = np.matmul(C, w)
    norm_vals[i] = np.linalg.norm(Z)

# ===================== Theoretical Chi distribution =====================

x = np.linspace(0, max(norm_vals) * 1.05, 1000)

# 1.  PDF
pdf_theory = (1 / a) * chi.pdf(x / a, k)
# 2. CDF
cdf_theory = chi.cdf(x / a, k)
# 3.  CCDF (Pr[X > x])
ccdf_theory = chi.sf(x / a, k)


plt.figure(figsize=(18, 5))

# ----------------- First figure：PDF -----------------
plt.subplot(131)
counts, bins, patches = plt.hist(norm_vals, bins=50, density=True, alpha=0.5, color='steelblue', label='Empirical Histogram')
bin_centers = 0.5 * (bins[:-1] + bins[1:])
smooth_counts = gaussian_filter1d(counts, sigma=1.2)
plt.plot(bin_centers, smooth_counts, color='darkblue', linewidth=2, label='Smoothed Empirical PDF')
plt.plot(x, pdf_theory, 'r-', linewidth=2, label='Analytical Chi PDF')
plt.xlabel(r'$\mathbb{Y} = \|\mathbf{C}\boldsymbol{\omega}\|_2$')
plt.ylabel('Probability Density')
plt.legend()
plt.title(rf'PDF Comparison')

# ----------------- Second figure：CDF -----------------
plt.subplot(132)
plt.hist(norm_vals, bins=50, density=True, cumulative=True, alpha=0.5, color='steelblue', label='Empirical CDF')
plt.plot(x, cdf_theory, 'r-', linewidth=2, label='Analytical Chi CDF')
plt.xlabel(r'$\mathbb{Y} = \|\mathbf{C}\boldsymbol{\omega}\|_2$')
plt.ylabel('Cumulative Probability')
plt.legend()
plt.title('CDF Comparison')

# ----------------- Third figure：CCDF  -----------------
plt.subplot(133)

sorted_norms = np.sort(norm_vals)

y_ccdf_empirical = 1.0 - np.arange(1, len(sorted_norms) + 1) / len(sorted_norms)

plt.plot(sorted_norms[:-1], y_ccdf_empirical[:-1], color='darkblue', linewidth=2, label='Empirical CCDF')

plt.plot(x, ccdf_theory, 'r-', linewidth=2, label='Analytical Chi CCDF')

plt.yscale('log')
plt.xlabel(r'$\mathbb{Y} = \|\mathbf{C}\boldsymbol{\omega}\|_2$')
plt.ylabel(r'Tail Probability: $\Pr[\mathbb{Y} > y]$ (Log Scale)')


plt.ylim(bottom=1e-6, top=1.5)
plt.legend()
plt.title('Log-Scale Tail (CCDF) Comparison')



plt.savefig("norm_distribution_analysis.pdf", format="pdf", bbox_inches="tight")

# ==========================================

plt.show(block=False) # 这行必须在 savefig 的后面

# ===================== 输出结果 =====================
print("=" * 60)
print(f"📌 w ：discrete Gaussian， center 0，sigma = {sigma_gauss}")
print(f"✅ ||w|| = {np.linalg.norm(w):.2f}")
print(f"✅ Chi distribution: scaling factor a = {a:.2f}, freedom k = {k}")
print(f"✅ Theoretical  E[||Cw||] = {a * chi.mean(k):.2f}")
print(f"✅ Experimental E[||Cw||] = {np.mean(norm_vals):.2f}")
print(f"✅ Theoretical Var[||Cw||] = {a ** 2 * chi.var(k):.2f}")
print(f"✅ Experimental Var[||Cw||] = {np.var(norm_vals):.2f}")
print("=" * 60)

plt.pause(0.01)
input("Press Enter to exit...")
plt.close('all')