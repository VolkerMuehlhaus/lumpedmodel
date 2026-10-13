* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py 1.2 from inductor_500p.s2p
* Fit band 0.010 to 78.000 GHz, 3 coil segment(s)
* Total element values:
*   Rs = 1.5926
*   Ls = 3.6692e-10
*   Rskin1 = 1.6564
*   Lskin1 = 2.58086e-11
*   Rskin2 = 127.767
*   Lskin2 = 4.49805e-11
*   Cs = 4.8369e-15
*   Cox1 = 1.49691e-14
*   Rsi1 = 578.815
*   Csi1 = 5.55904e-14
*   Cox2 = 1.49691e-14
*   Rsi2 = 578.815
*   Csi2 = 5.55904e-14
*
.SUBCKT inductor_model p1 p2
* coil segment 1
Rs_1 p1 a1 0.530868
Ls_1 a1 b1 1.22307e-10
Rskin1_1 b1 d1 0.552134
Lskin1_1 b1 d1 8.60287e-12
Rskin2_1 d1 c1 42.589
Lskin2_1 d1 c1 1.49935e-11
* coil segment 2
Rs_2 c1 a2 0.530868
Ls_2 a2 b2 1.22307e-10
Rskin1_2 b2 d2 0.552134
Lskin1_2 b2 d2 8.60287e-12
Rskin2_2 d2 c2 42.589
Lskin2_2 d2 c2 1.49935e-11
* coil segment 3
Rs_3 c2 a3 0.530868
Ls_3 a3 b3 1.22307e-10
Rskin1_3 b3 d3 0.552134
Lskin1_3 b3 d3 8.60287e-12
Rskin2_3 d3 p2 42.589
Lskin2_3 d3 p2 1.49935e-11
Cs p1 p2 4.8369e-15
* substrate network at p1
Cox_n0 p1 s0 4.98969e-15
Rsi_n0 s0 0 1736.44
Csi_n0 s0 0 1.85301e-14
* substrate network at c1
Cox_n1 c1 s1 9.97939e-15
Rsi_n1 s1 0 868.222
Csi_n1 s1 0 3.70603e-14
* substrate network at c2
Cox_n2 c2 s2 9.97939e-15
Rsi_n2 s2 0 868.222
Csi_n2 s2 0 3.70603e-14
* substrate network at p2
Cox_n3 p2 s3 4.98969e-15
Rsi_n3 s3 0 1736.44
Csi_n3 s3 0 1.85301e-14
.ENDS
