* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py 1.2 from inductor_74p.s2p
* Fit band 0.010 to 204.000 GHz, 3 coil segment(s)
* Total element values:
*   Rs = 0.350681
*   Ls = 6.88429e-11
*   Rskin1 = 0.698425
*   Lskin1 = 2.85095e-11
*   Rskin2 = 1.95069
*   Lskin2 = 2.53761e-12
*   Cs = 1e-19
*   Cox1 = 1.49622e-14
*   Rsi1 = 0.139372
*   Csi1 = 2.86739e-19
*   Cox2 = 1.49622e-14
*   Rsi2 = 0.139372
*   Csi2 = 2.86739e-19
*
.SUBCKT inductor_model p1 p2
* coil segment 1
Rs_1 p1 a1 0.116894
Ls_1 a1 b1 2.29476e-11
Rskin1_1 b1 d1 0.232808
Lskin1_1 b1 d1 9.50315e-12
Rskin2_1 d1 c1 0.650229
Lskin2_1 d1 c1 8.4587e-13
* coil segment 2
Rs_2 c1 a2 0.116894
Ls_2 a2 b2 2.29476e-11
Rskin1_2 b2 d2 0.232808
Lskin1_2 b2 d2 9.50315e-12
Rskin2_2 d2 c2 0.650229
Lskin2_2 d2 c2 8.4587e-13
* coil segment 3
Rs_3 c2 a3 0.116894
Ls_3 a3 b3 2.29476e-11
Rskin1_3 b3 d3 0.232808
Lskin1_3 b3 d3 9.50315e-12
Rskin2_3 d3 p2 0.650229
Lskin2_3 d3 p2 8.4587e-13
Cs p1 p2 1e-19
* substrate network at p1
Cox_n0 p1 s0 4.9874e-15
Rsi_n0 s0 0 0.418116
Csi_n0 s0 0 9.55796e-20
* substrate network at c1
Cox_n1 c1 s1 9.97481e-15
Rsi_n1 s1 0 0.209058
Csi_n1 s1 0 1.91159e-19
* substrate network at c2
Cox_n2 c2 s2 9.97481e-15
Rsi_n2 s2 0 0.209058
Csi_n2 s2 0 1.91159e-19
* substrate network at p2
Cox_n3 p2 s3 4.9874e-15
Rsi_n3 s3 0 0.418116
Csi_n3 s3 0 9.55796e-20
.ENDS
