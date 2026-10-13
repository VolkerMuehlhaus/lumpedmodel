* WIDEBAND RFIC CENTER-TAPPED INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py 1.2 from final_inductor3_N2_do82.29_w4_s4.s3p
* Fit band 0.010 to 60.000 GHz, 2 coil segment(s) per half coil
* k is the effective coupling including the center tap lead inductance
* Total element values:
*   Rs_h1 = 0.719309
*   Ls_h1 = 1.52407e-10
*   Rskin1_h1 = 0.682082
*   Lskin1_h1 = 1.08001e-11
*   Rskin2_h1 = 5.25196
*   Lskin2_h1 = 9.68558e-12
*   Rs_h2 = 0.998181
*   Ls_h2 = 1.52407e-10
*   Rskin1_h2 = 1.0548
*   Lskin1_h2 = 1.08001e-11
*   Rskin2_h2 = 8.31206
*   Lskin2_h2 = 9.68558e-12
*   k = 0.409973
*   Cs = 4.84887e-15
*   Rct = 0.0956391
*   Cox1 = 6.48262e-15
*   Rsi1 = 526.78
*   Csi1 = 5.57673e-14
*   Cox2 = 7.1869e-15
*   Rsi2 = 561.786
*   Csi2 = 5.38027e-14
*   Coxct = 1.90342e-14
*   Rsict = 8343.33
*   Csict = 2.92926e-15
*
.SUBCKT ct_inductor_model p1 p2 ct
* half coil 1, segment 1
Rs_h1_1 p1 a_h1_1 0.359655
Ls_h1_1 a_h1_1 b_h1_1 7.62034e-11
Rskin1_h1_1 b_h1_1 d_h1_1 0.341041
Lskin1_h1_1 b_h1_1 d_h1_1 5.40005e-12
Rskin2_h1_1 d_h1_1 c1 2.62598
Lskin2_h1_1 d_h1_1 c1 4.84279e-12
* half coil 1, segment 2
Rs_h1_2 c1 a_h1_2 0.359655
Ls_h1_2 a_h1_2 b_h1_2 7.62034e-11
Rskin1_h1_2 b_h1_2 d_h1_2 0.341041
Lskin1_h1_2 b_h1_2 d_h1_2 5.40005e-12
Rskin2_h1_2 d_h1_2 m 2.62598
Lskin2_h1_2 d_h1_2 m 4.84279e-12
* half coil 2, segment 1
Rs_h2_1 m a_h2_1 0.499091
Ls_h2_1 a_h2_1 b_h2_1 7.62034e-11
Rskin1_h2_1 b_h2_1 d_h2_1 0.527398
Lskin1_h2_1 b_h2_1 d_h2_1 5.40005e-12
Rskin2_h2_1 d_h2_1 c3 4.15603
Lskin2_h2_1 d_h2_1 c3 4.84279e-12
* half coil 2, segment 2
Rs_h2_2 c3 a_h2_2 0.499091
Ls_h2_2 a_h2_2 b_h2_2 7.62034e-11
Rskin1_h2_2 b_h2_2 d_h2_2 0.527398
Lskin1_h2_2 b_h2_2 d_h2_2 5.40005e-12
Rskin2_h2_2 d_h2_2 p2 4.15603
Lskin2_h2_2 d_h2_2 p2 4.84279e-12
* magnetic coupling of the half coils (dot at the first node, aiding for differential current)
K_h1_1_h2_1 Ls_h1_1 Ls_h2_1 0.204986
K_h1_1_h2_2 Ls_h1_1 Ls_h2_2 0.204986
K_h1_2_h2_1 Ls_h1_2 Ls_h2_1 0.204986
K_h1_2_h2_2 Ls_h1_2 Ls_h2_2 0.204986
Cs p1 p2 4.84887e-15
Rct m ct 0.0956391
* substrate network at p1
Cox_n0 p1 s0 3.24131e-15
Rsi_n0 s0 0 1053.56
Csi_n0 s0 0 2.78837e-14
* substrate network at c1
Cox_n1 c1 s1 7.99987e-15
Rsi_n1 s1 0 1021.32
Csi_n1 s1 0 2.8616e-14
* substrate network at m
Cox_n2 m s2 9.51712e-15
Rsi_n2 s2 0 16686.7
Csi_n2 s2 0 1.46463e-15
* substrate network at c3
Cox_n3 c3 s3 8.35201e-15
Rsi_n3 s3 0 1086.98
Csi_n3 s3 0 2.76337e-14
* substrate network at p2
Cox_n4 p2 s4 3.59345e-15
Rsi_n4 s4 0 1123.57
Csi_n4 s4 0 2.69014e-14
.ENDS
