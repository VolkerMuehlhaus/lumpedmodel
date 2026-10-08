* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py from sample_inductor.s2p
* Fit band 0.075 to 24.975 GHz, 2 coil segment(s)
* Total element values:
*   Rs = 1.54228
*   Ls = 1.567e-09
*   Rskin1 = 0.929686
*   Lskin1 = 4.78824e-11
*   Rskin2 = 10.5702
*   Lskin2 = 5.19172e-11
*   Cs = 2.52476e-14
*   Cox1 = 5.94996e-14
*   Rsi1 = 3533.02
*   Csi1 = 1.83229e-14
*   Cox2 = 5.94996e-14
*   Rsi2 = 3533.02
*   Csi2 = 1.83229e-14
*
.SUBCKT inductor_model p1 p2
* coil segment 1
Rs_1 p1 a1 0.771141
Ls_1 a1 b1 7.835e-10
Rskin1_1 b1 d1 0.464843
Lskin1_1 b1 d1 2.39412e-11
Rskin2_1 d1 c1 5.28509
Lskin2_1 d1 c1 2.59586e-11
* coil segment 2
Rs_2 c1 a2 0.771141
Ls_2 a2 b2 7.835e-10
Rskin1_2 b2 d2 0.464843
Lskin1_2 b2 d2 2.39412e-11
Rskin2_2 d2 p2 5.28509
Lskin2_2 d2 p2 2.59586e-11
Cs p1 p2 2.52476e-14
* substrate network at p1
Cox_n0 p1 s0 2.97498e-14
Rsi_n0 s0 0 7066.05
Csi_n0 s0 0 9.16147e-15
* substrate network at c1
Cox_n1 c1 s1 5.94996e-14
Rsi_n1 s1 0 3533.02
Csi_n1 s1 0 1.83229e-14
* substrate network at p2
Cox_n2 p2 s2 2.97498e-14
Rsi_n2 s2 0 7066.05
Csi_n2 s2 0 9.16147e-15
.ENDS
