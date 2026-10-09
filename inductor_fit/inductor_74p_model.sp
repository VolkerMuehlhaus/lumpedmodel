* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py from inductor_74p.s2p
* Fit band 0.010 to 204.000 GHz, 3 coil segment(s)
* Total element values:
*   Rs = 0.347256
*   Ls = 6.88557e-11
*   Rskin1 = 0.695967
*   Lskin1 = 2.97205e-11
*   Rskin2 = 1.93449
*   Lskin2 = 2.55548e-12
*   Cs = 1e-19
*   Cox1 = 1.49625e-14
*   Rsi1 = 0.139324
*   Csi1 = 2.86739e-19
*   Cox2 = 1.49625e-14
*   Rsi2 = 0.139324
*   Csi2 = 2.86739e-19
*
.SUBCKT inductor_model p1 p2
* coil segment 1
Rs_1 p1 a1 0.115752
Ls_1 a1 b1 2.29519e-11
Rskin1_1 b1 d1 0.231989
Lskin1_1 b1 d1 9.90683e-12
Rskin2_1 d1 c1 0.644831
Lskin2_1 d1 c1 8.51827e-13
* coil segment 2
Rs_2 c1 a2 0.115752
Ls_2 a2 b2 2.29519e-11
Rskin1_2 b2 d2 0.231989
Lskin1_2 b2 d2 9.90683e-12
Rskin2_2 d2 c2 0.644831
Lskin2_2 d2 c2 8.51827e-13
* coil segment 3
Rs_3 c2 a3 0.115752
Ls_3 a3 b3 2.29519e-11
Rskin1_3 b3 d3 0.231989
Lskin1_3 b3 d3 9.90683e-12
Rskin2_3 d3 p2 0.644831
Lskin2_3 d3 p2 8.51827e-13
Cs p1 p2 1e-19
* substrate network at p1
Cox_n0 p1 s0 4.98749e-15
Rsi_n0 s0 0 0.417972
Csi_n0 s0 0 9.55796e-20
* substrate network at c1
Cox_n1 c1 s1 9.97498e-15
Rsi_n1 s1 0 0.208986
Csi_n1 s1 0 1.91159e-19
* substrate network at c2
Cox_n2 c2 s2 9.97498e-15
Rsi_n2 s2 0 0.208986
Csi_n2 s2 0 1.91159e-19
* substrate network at p2
Cox_n3 p2 s3 4.98749e-15
Rsi_n3 s3 0 0.417972
Csi_n3 s3 0 9.55796e-20
.ENDS
