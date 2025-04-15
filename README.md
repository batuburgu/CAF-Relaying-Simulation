# CAF-Relaying-Simulation

This repository provides MATLAB simulations for analyzing the performance of **Cooperative Amplify-and-Forward (CAF) relaying** systems under **Nakagami-m fading channels**. The simulations cover various scenarios including multi-branch relaying, SER under unequal channel parameters, and relay diversity techniques.

---

## 📁 File Overview and References

### 1. `different_channels_caf_ser.m`

Simulates the **Symbol Error Rate (SER)** of a CAF system under **Nakagami-m fading** where the channels between source-relay and relay-destination have **different m parameters**.

📄 **Reference**:  
> A. Agustin and J. Vidal,  
> *Symbol error probability analysis of cooperative diversity protocols under Nakagami-m fading*,  
> IEEE Communications Letters, vol. 9, no. 1, pp. 20–22, Jan. 2005.  
> [DOI: 10.1109/LCOMM.2004.839952](https://doi.org/10.1109/LCOMM.2004.839952)

📊 **Reproduced Content**:  
- SER vs. SNR plots for asymmetric fading conditions

---

### 2. `multibranch_caf_ber_approx.m`

Analyzes the **BER performance of multi-branch CAF systems** using an **analytical approximation** method under Nakagami-\( m \) fading conditions.

📄 **Reference**:  
> S. S. Ikki and M. H. Ahmed,  
> *Performance analysis of cooperative diversity wireless networks over Nakagami-m fading channel*,  
> IEEE Communications Letters, vol. 11, no. 4, pp. 334–336, Apr. 2007.  
> [DOI: 10.1109/LCOMM.2007.061420](https://doi.org/10.1109/LCOMM.2007.061420)

📊 **Reproduced Content**:  
- BER approximation curves for various number of relays

---

### 3. `multirelay_diversity.m`

Explores the **diversity gain** in CAF systems by evaluating the effect of increasing the number of relays on BER performance.

📄 **Reference**:  
> M. O. Hasna and M. S. Alouini,  
> *Harmonic mean and end-to-end performance of transmission systems with relays*,  
> IEEE Transactions on Communications, vol. 52, no. 1, pp. 130–135, Jan. 2004.  
> [DOI: 10.1109/TCOMM.2003.822168](https://doi.org/10.1109/TCOMM.2003.822168)

📊 **Reproduced Content**:  
- BER vs. SNR curves showing improved diversity with increasing relays

