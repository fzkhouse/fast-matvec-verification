import numpy as np
from scipy.stats import chi
from scipy.integrate import quad
import matplotlib.pyplot as plt
import time


def chi_integral(B2, k, M):

    sqrt_60 = np.sqrt(60)
    sqrt_674 = np.sqrt(674)
    sqrt_30_337 = np.sqrt(30 / 337)
    sqrt_337_30 = np.sqrt(337 / 30)
    jacobi = 1 / M

    t_lower = sqrt_30_337 * B2
    t_upper = sqrt_337_30 * B2

    def integrand(t):
        a = sqrt_60 * B2 / t
        b = sqrt_674 * B2 / t
        prob_y = chi.cdf(b, 256) - chi.cdf(a, 256)
        prob_y = np.clip(prob_y, 0.0, 1.0)
        pdf_t = chi.pdf(t / M, k) * jacobi
        integrand_val = prob_y * pdf_t
        integrand_val = np.clip(integrand_val, 0.0, 1.0)
        return integrand_val


    try:
        integral_value, error = quad(
            func=integrand,
            a=t_lower,
            b=t_upper,
            epsabs=1e-4,
            epsrel=1e-4,
            limit=50,
            maxp1=50
        )[0:2]
    except Exception as e:
        print(f"⚠️ B2={B2} Integral Error：{e}，Return 0")
        integral_value, error = 0.0, 1.0

    return integral_value, error


def find_first_positive_b2(B2_list, integral_results, epsilon=1e-8):

    for idx, (b2, integral) in enumerate(zip(B2_list, integral_results)):
        if integral > epsilon:
            return round(b2, 1)
    return None



if __name__ == "__main__":
    k_list = [256, 512, 1024]
    M_list = np.arange(20, 21, 1)


    plot_data = {}
    step = 1
    zero_integral_results = []
    first_positive_b2_results = {}
    start_time = time.time()


    for k in k_list:
        for M in M_list:
            print(f"\n{'=' * 70}")
            print(f"Current Parameters: k = {k}  |  M = {M}")
            print(f"{'=' * 70}")


            B2_start = 0.0001 * M * np.sqrt(k)
            B2_end = M * np.sqrt(k)
            B2_list = np.arange(np.ceil(B2_start), np.floor(B2_end) + step, step)
            if len(B2_list) == 0:
                B2_list = np.array([B2_start, B2_end])

            print(f"📊 B2 Range：{B2_list[0]:.0f} ~ {B2_list[-1]:.0f}，{len(B2_list)} points in total")
            integral_results = []

            total_B2 = len(B2_list)
            for idx, B2 in enumerate(B2_list):
                if idx % 10 == 0:
                    elapsed = time.time() - start_time
                    print(f"⏳ Progress：{idx}/{total_B2} B2 Values. Time: {elapsed:.1f} seconds")

                res, err = chi_integral(B2, k, M)
                integral_results.append(res)
                print(f"B2 = {B2:4.1f} | Integral Result = {res:.8f} | Numerical Error = {err:.2e}")


                res_rounded = round(res, 8)
                if res_rounded == 0.0:

                    C_value = B2 / M

                    zero_integral_results.append({
                        "k": k,
                        "M": M,
                        "B2": round(B2, 1),
                        "C": C_value,
                        "integral_value": res_rounded
                    })

            first_positive_b2 = find_first_positive_b2(B2_list, integral_results)
            first_positive_b2_results[(k, M)] = first_positive_b2


            plot_data[(k, M)] = {
                "B2": [round(b, 1) for b in B2_list],
                "integral": integral_results,
                "first_positive_b2": first_positive_b2
            }

            max_integral = np.max(integral_results)
            max_B2_index = np.argmax(integral_results)
            max_B2 = B2_list[max_B2_index]
            print(f"\n✅ Calculation completed! Max Integral = {max_integral:.8f}, B2 = {max_B2:.1f}")


    print(f"\n\n{'=' * 90}")
    print("📌 Summary of B2 and C Values with Integral Result = 0.00000000 📌")
    print(f"{'=' * 90}")
    if len(zero_integral_results) == 0:
        print("⚠️  No B2 values with integral result = 0.00000000 found!")
    else:
        for idx, item in enumerate(zero_integral_results, 1):
            print(
                f"{idx}. k={item['k']}, M={item['M']} → B2 = {item['B2']:.1f} | C = {item['C']:.4f} | Integral Value = {item['integral_value']:.8f}")

    print(f"\n\n{'=' * 80}")
    print("📌 Summary of First Positive Integral B2 Values 📌")
    print(f"{'=' * 80}")
    for (k, M), b2_val in first_positive_b2_results.items():
        if b2_val is not None:
            C_critical = b2_val / M
            print(f"k={k}, M={M} → First Positive Integral B2 Values   = {b2_val:.1f} (C = {C_critical:.4f})")
        else:
            print(f"k={k}, M={M} → No Positive Integral B2 Values found!")
    print(f"{'=' * 80}")


    print("\n\n🎨 Start plotting the relationship between B2 and integral value...")
    try:
        fig, axes = plt.subplots(len(plot_data), 1, figsize=(12, 6 * len(plot_data)))
        if len(plot_data) == 1:
            axes = [axes]

        for idx, ((k, M), data) in enumerate(plot_data.items()):
            ax = axes[idx]
            B2_vals = data["B2"]
            integral_vals = data["integral"]
            first_positive_b2 = data["first_positive_b2"]

            ax.plot(B2_vals, integral_vals,
                    color='#1f77b4', linewidth=2, marker='o', markersize=4,
                    label=f'k={k}, M={M} Integral Value')

            if first_positive_b2 is not None:
                b2_idx = np.where(np.array(B2_vals) == first_positive_b2)[0][0]
                integral_at_first = integral_vals[b2_idx]

                ax.axvline(x=first_positive_b2, color='red', linestyle='--', linewidth=2, alpha=0.8,
                           label=f'First Positive B2 = {first_positive_b2}')
                ax.scatter(first_positive_b2, integral_at_first, color='red', s=100, zorder=5)
                ax.annotate(f'First Positive\nB2 = {first_positive_b2}',
                            xy=(first_positive_b2, integral_at_first),
                            xytext=(10, 20), textcoords='offset points',
                            fontsize=10, color='red', fontweight='bold',
                            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='red', alpha=0.8))

            ax.set_title(f'Relationship between B2 and Integral Value (k={k}, M={M})', fontsize=14, pad=10)
            ax.set_xlabel('B2 Value', fontsize=12)
            ax.set_ylabel('Integral Result (E(abort))', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best', fontsize=10)
            ax.set_xlim(min(B2_vals) - 5, max(B2_vals) + 5)

        plt.tight_layout()
        plt.savefig('B2_Integral_Relationship.png', dpi=300, bbox_inches='tight')
        print("✅ Image saved as 'B2_Integral_Relationship.png' (no need to wait for window close)")

        plt.ion()
        plt.show(block=False)
        plt.pause(5)
        plt.close('all')

    except Exception as e:
        print(f"⚠️ Plotting error: {e}")

    total_elapsed = time.time() - start_time
    print(f"\n🎉 All tasks completed! Total elapsed time: {total_elapsed:.2f} seconds")
    print("✅ Program terminated successfully!")