* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py 1.2 from sample_inductor.s2p
* Fit band 0.075 to 24.975 GHz, 2 coil segment(s)
* Total element values:
*   Rs = 1.54203
*   Ls = 1.56782e-09
*   Rskin1 = 0.922236
*   Lskin1 = 4.77875e-11
*   Rskin2 = 10.1768
*   Lskin2 = 5.125e-11
*   Cs = 2.52537e-14
*   Cox1 = 5.94999e-14
*   Rsi1 = 3532.97
*   Csi1 = 1.8323e-14
*   Cox2 = 5.94999e-14
*   Rsi2 = 3532.97
*   Csi2 = 1.8323e-14
*
.SUBCKT inductor_model p1 p2
* coil segment 1
Rs_1 p1 a1 0.771016
Ls_1 a1 b1 7.83909e-10
Rskin1_1 b1 d1 0.461118
Lskin1_1 b1 d1 2.38937e-11
Rskin2_1 d1 c1 5.08838
Lskin2_1 d1 c1 2.5625e-11
* coil segment 2
Rs_2 c1 a2 0.771016
Ls_2 a2 b2 7.83909e-10
Rskin1_2 b2 d2 0.461118
Lskin1_2 b2 d2 2.38937e-11
Rskin2_2 d2 p2 5.08838
Lskin2_2 d2 p2 2.5625e-11
Cs p1 p2 2.52537e-14
* substrate network at p1
Cox_n0 p1 s0 2.975e-14
Rsi_n0 s0 0 7065.94
Csi_n0 s0 0 9.16148e-15
* substrate network at c1
Cox_n1 c1 s1 5.94999e-14
Rsi_n1 s1 0 3532.97
Csi_n1 s1 0 1.8323e-14
* substrate network at p2
Cox_n2 p2 s2 2.975e-14
Rsi_n2 s2 0 7065.94
Csi_n2 s2 0 9.16148e-15
.ENDS
