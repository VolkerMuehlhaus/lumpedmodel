* WIDEBAND RFIC CENTER-TAPPED INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py 1.2 from final_inductor3_N2_do82.29_w4_s4.s3p
* Fit band 0.010 to 60.000 GHz, 2 coil segment(s) per half coil
* k is the effective coupling including the center tap lead inductance
* Total element values:
*   Rs_h1 = 0.720612
*   Ls_h1 = 1.50842e-10
*   Rskin1_h1 = 0.758223
*   Lskin1_h1 = 1.10383e-11
*   Rskin2_h1 = 7.63857
*   Lskin2_h1 = 1.07664e-11
*   Rs_h2 = 0.998512
*   Ls_h2 = 1.50842e-10
*   Rskin1_h2 = 1.13045
*   Lskin1_h2 = 1.10383e-11
*   Rskin2_h2 = 11.4561
*   Lskin2_h2 = 1.07664e-11
*   k = 0.414297
*   Cs = 4.84023e-15
*   Rct = 0.0957629
*   Cox1 = 6.47696e-15
*   Rsi1 = 526.49
*   Csi1 = 5.58014e-14
*   Cox2 = 7.1815e-15
*   Rsi2 = 561.325
*   Csi2 = 5.38498e-14
*   Coxct = 1.90625e-14
*   Rsict = 8355.23
*   Csict = 2.91883e-15
*
.SUBCKT ct_inductor_model p1 p2 ct
* half coil 1, segment 1
Rs_h1_1 p1 a_h1_1 0.360306
Ls_h1_1 a_h1_1 b_h1_1 7.5421e-11
Rskin1_h1_1 b_h1_1 d_h1_1 0.379111
Lskin1_h1_1 b_h1_1 d_h1_1 5.51915e-12
Rskin2_h1_1 d_h1_1 c1 3.81928
Lskin2_h1_1 d_h1_1 c1 5.38322e-12
* half coil 1, segment 2
Rs_h1_2 c1 a_h1_2 0.360306
Ls_h1_2 a_h1_2 b_h1_2 7.5421e-11
Rskin1_h1_2 b_h1_2 d_h1_2 0.379111
Lskin1_h1_2 b_h1_2 d_h1_2 5.51915e-12
Rskin2_h1_2 d_h1_2 m 3.81928
Lskin2_h1_2 d_h1_2 m 5.38322e-12
* half coil 2, segment 1
Rs_h2_1 m a_h2_1 0.499256
Ls_h2_1 a_h2_1 b_h2_1 7.5421e-11
Rskin1_h2_1 b_h2_1 d_h2_1 0.565224
Lskin1_h2_1 b_h2_1 d_h2_1 5.51915e-12
Rskin2_h2_1 d_h2_1 c3 5.72805
Lskin2_h2_1 d_h2_1 c3 5.38322e-12
* half coil 2, segment 2
Rs_h2_2 c3 a_h2_2 0.499256
Ls_h2_2 a_h2_2 b_h2_2 7.5421e-11
Rskin1_h2_2 b_h2_2 d_h2_2 0.565224
Lskin1_h2_2 b_h2_2 d_h2_2 5.51915e-12
Rskin2_h2_2 d_h2_2 p2 5.72805
Lskin2_h2_2 d_h2_2 p2 5.38322e-12
* magnetic coupling of the half coils (dot at the first node, aiding for differential current)
K_h1_1_h2_1 Ls_h1_1 Ls_h2_1 0.207148
K_h1_1_h2_2 Ls_h1_1 Ls_h2_2 0.207148
K_h1_2_h2_1 Ls_h1_2 Ls_h2_1 0.207148
K_h1_2_h2_2 Ls_h1_2 Ls_h2_2 0.207148
Cs p1 p2 4.84023e-15
Rct m ct 0.0957629
* substrate network at p1
Cox_n0 p1 s0 3.23848e-15
Rsi_n0 s0 0 1052.98
Csi_n0 s0 0 2.79007e-14
* substrate network at c1
Cox_n1 c1 s1 8.00411e-15
Rsi_n1 s1 0 1020.82
Csi_n1 s1 0 2.86304e-14
* substrate network at m
Cox_n2 m s2 9.53127e-15
Rsi_n2 s2 0 16710.5
Csi_n2 s2 0 1.45942e-15
* substrate network at c3
Cox_n3 c3 s3 8.35638e-15
Rsi_n3 s3 0 1086.16
Csi_n3 s3 0 2.76546e-14
* substrate network at p2
Cox_n4 p2 s4 3.59075e-15
Rsi_n4 s4 0 1122.65
Csi_n4 s4 0 2.69249e-14
.ENDS
