# CAF-Relaying-Simulation

This repository provides MATLAB simulations for analyzing the performance of **Cooperative Amplify-and-Forward (CAF) relaying** systems under **Nakagami-m fading channels**. The simulations cover various scenarios including multi-branch relaying, SER under unequal channel parameters, and relay diversity techniques.

---

## 📁 File Overview and References

### 1. `different_channels_caf_ser.m`

Simulates the **Symbol Error Rate (SER)** of a CAF system under **Nakagami-m fading** where the channels between source-relay and relay-destination have **different m parameters**.

📄 **Reference**:  
> Yang, L.-L., & Chen, H.-H.
> *Error probability of digital communications using relay diversity over Nakagami-m fading channels*,  
> IEEE Transactions on Wireless Communications, 7(5), 1806–1811, 2008.  

📊 **Reproduced Content**:  
- Figure 3. (only BPSK)

---

### 2. `multibranch_caf_ber_approx.m`

Analyzes the **BER performance of multi-branch CAF systems** using an **analytical approximation** method under Nakagami-\( m \) fading conditions.

📄 **Reference**:  
> S. S. Ikki and M. H. Ahmed,  
> *Performance analysis of cooperative diversity wireless networks over Nakagami-m fading channel*,  
> IEEE Communications Letters, vol. 11, no. 4, pp. 334–336, Apr. 2007.  

📊 **Reproduced Content**:  
- Figure 2.

---

### 3. `multirelay_diversity.m`

Explores the **diversity gain** in CAF systems by evaluating the effect of increasing the number of relays on BER performance.

📄 **Reference**:  
> Anghel, P. A., & Kaveh, M.,  
> *Exact symbol error probability of a cooperative network in a Rayleigh-fading environment*,  
> IEEE Transactions on Wireless Communications, 3(5), 1416–1421, 2004.  

📊 **Reproduced Content**:  
- Figure 2.

