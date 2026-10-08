* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py from inductor_500p.s2p
* Fit band 0.010 to 78.000 GHz, 3 coil segment(s)
* Total element values:
*   Rs = 1.59222
*   Ls = 2.58542e-10
*   Rskin1 = 1.65955
*   Lskin1 = 2.58672e-11
*   Rskin2 = 1494.2
*   Lskin2 = 1.53354e-10
*   Cs = 4.81208e-15
*   Cox1 = 1.49678e-14
*   Rsi1 = 579.486
*   Csi1 = 5.55691e-14
*   Cox2 = 1.49678e-14
*   Rsi2 = 579.486
*   Csi2 = 5.55691e-14
*
.SUBCKT inductor_model p1 p2
* coil segment 1
Rs_1 p1 a1 0.530739
Ls_1 a1 b1 8.61807e-11
Rskin1_1 b1 d1 0.553183
Lskin1_1 b1 d1 8.62241e-12
Rskin2_1 d1 c1 498.068
Lskin2_1 d1 c1 5.11179e-11
* coil segment 2
Rs_2 c1 a2 0.530739
Ls_2 a2 b2 8.61807e-11
Rskin1_2 b2 d2 0.553183
Lskin1_2 b2 d2 8.62241e-12
Rskin2_2 d2 c2 498.068
Lskin2_2 d2 c2 5.11179e-11
* coil segment 3
Rs_3 c2 a3 0.530739
Ls_3 a3 b3 8.61807e-11
Rskin1_3 b3 d3 0.553183
Lskin1_3 b3 d3 8.62241e-12
Rskin2_3 d3 p2 498.068
Lskin2_3 d3 p2 5.11179e-11
Cs p1 p2 4.81208e-15
* substrate network at p1
Cox_n0 p1 s0 4.98925e-15
Rsi_n0 s0 0 1738.46
Csi_n0 s0 0 1.8523e-14
* substrate network at c1
Cox_n1 c1 s1 9.9785e-15
Rsi_n1 s1 0 869.229
Csi_n1 s1 0 3.7046e-14
* substrate network at c2
Cox_n2 c2 s2 9.9785e-15
Rsi_n2 s2 0 869.229
Csi_n2 s2 0 3.7046e-14
* substrate network at p2
Cox_n3 p2 s3 4.98925e-15
Rsi_n3 s3 0 1738.46
Csi_n3 s3 0 1.8523e-14
.ENDS
