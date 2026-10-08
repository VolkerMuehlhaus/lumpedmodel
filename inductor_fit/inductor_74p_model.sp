* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py from inductor_74p.s2p
* Fit band 0.010 to 204.000 GHz, 3 coil segment(s)
* Total element values:
*   Rs = 0.34695
*   Ls = 6.88546e-11
*   Rskin1 = 0.695398
*   Lskin1 = 2.98176e-11
*   Rskin2 = 1.93502
*   Lskin2 = 2.55872e-12
*   Cs = 1e-19
*   Cox1 = 1.49625e-14
*   Rsi1 = 0.139288
*   Csi1 = 2.86739e-19
*   Cox2 = 1.49625e-14
*   Rsi2 = 0.139288
*   Csi2 = 2.86739e-19
*
.SUBCKT inductor_model p1 p2
* coil segment 1
Rs_1 p1 a1 0.11565
Ls_1 a1 b1 2.29515e-11
Rskin1_1 b1 d1 0.231799
Lskin1_1 b1 d1 9.93919e-12
Rskin2_1 d1 c1 0.645006
Lskin2_1 d1 c1 8.52907e-13
* coil segment 2
Rs_2 c1 a2 0.11565
Ls_2 a2 b2 2.29515e-11
Rskin1_2 b2 d2 0.231799
Lskin1_2 b2 d2 9.93919e-12
Rskin2_2 d2 c2 0.645006
Lskin2_2 d2 c2 8.52907e-13
* coil segment 3
Rs_3 c2 a3 0.11565
Ls_3 a3 b3 2.29515e-11
Rskin1_3 b3 d3 0.231799
Lskin1_3 b3 d3 9.93919e-12
Rskin2_3 d3 p2 0.645006
Lskin2_3 d3 p2 8.52907e-13
Cs p1 p2 1e-19
* substrate network at p1
Cox_n0 p1 s0 4.98749e-15
Rsi_n0 s0 0 0.417863
Csi_n0 s0 0 9.55796e-20
* substrate network at c1
Cox_n1 c1 s1 9.97498e-15
Rsi_n1 s1 0 0.208932
Csi_n1 s1 0 1.91159e-19
* substrate network at c2
Cox_n2 c2 s2 9.97498e-15
Rsi_n2 s2 0 0.208932
Csi_n2 s2 0 1.91159e-19
* substrate network at p2
Cox_n3 p2 s3 4.98749e-15
Rsi_n3 s3 0 0.417863
Csi_n3 s3 0 9.55796e-20
.ENDS
