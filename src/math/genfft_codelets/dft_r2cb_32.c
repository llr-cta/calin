/* Generated by: /Users/sfegan/GitHub/fftw3/genfft/gen_r2cb.native -n 32 -sign 1 -standalone -fma -generic-arith -compact -name dft_codelet_r2cb_32 */

/*
 * This function contains 156 FP additions, 84 FP multiplications,
 * (or, 72 additions, 0 multiplications, 84 fused multiply/add),
 * 54 stack variables, 9 constants, and 64 memory accesses
 */
void dft_codelet_r2cb_32(R * R0, R * R1, R * Cr, R * Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
  DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
  DK(KP668178637, +0.668178637919298919997757686523080761552472251);
  DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
  DK(KP198912367, +0.198912367379658006911597622644676228597850501);
  DK(KP707106781, +0.707106781186547524400844362104849039284835938);
  DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
  DK(KP414213562, +0.414213562373095048801688724209698078569671875);
  DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
  DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
  {
    INT i;
    for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(128, rs), MAKE_VOLATILE_STRIDE(128, csr), MAKE_VOLATILE_STRIDE(128, csi)) {
      E T5, T1R, Tz, T1t, T8, T1S, TE, T1u, Tg, T1X, T2m, TK, TP, T1x,
       T1U;
      E T1w, To, T28, T2p, TW, T1d, T1D, T20, T1A, Tv, T23, T2q, T25,
       T1g, T1B;
      E T17, T1E;
      {
        E T4, Ty, T3, Tx, T1, T2;
        T4 = Cr[WS(csr, 8)];
        Ty = Ci[WS(csi, 8)];
        T1 = Cr[0];
        T2 = Cr[WS(csr, 16)];
        T3 = ADD(T1, T2);
        Tx = SUB(T1, T2);
        T5 = FMA(KP2_000000000, T4, T3);
        T1R = FNMS(KP2_000000000, T4, T3);
        Tz = FNMS(KP2_000000000, Ty, Tx);
        T1t = FMA(KP2_000000000, Ty, Tx);
      }
      {
        E T6, T7, TA, TB, TC, TD;
        T6 = Cr[WS(csr, 4)];
        T7 = Cr[WS(csr, 12)];
        TA = SUB(T6, T7);
        TB = Ci[WS(csi, 4)];
        TC = Ci[WS(csi, 12)];
        TD = ADD(TB, TC);
        T8 = ADD(T6, T7);
        T1S = SUB(TB, TC);
        TE = SUB(TA, TD);
        T1u = ADD(TA, TD);
      }
      {
        E Tc, TG, TO, T1V, Tf, TL, TJ, T1W;
        {
          E Ta, Tb, TM, TN;
          Ta = Cr[WS(csr, 2)];
          Tb = Cr[WS(csr, 14)];
          Tc = ADD(Ta, Tb);
          TG = SUB(Ta, Tb);
          TM = Ci[WS(csi, 2)];
          TN = Ci[WS(csi, 14)];
          TO = ADD(TM, TN);
          T1V = SUB(TM, TN);
        }
        {
          E Td, Te, TH, TI;
          Td = Cr[WS(csr, 10)];
          Te = Cr[WS(csr, 6)];
          Tf = ADD(Td, Te);
          TL = SUB(Td, Te);
          TH = Ci[WS(csi, 10)];
          TI = Ci[WS(csi, 6)];
          TJ = ADD(TH, TI);
          T1W = SUB(TH, TI);
        }
        Tg = ADD(Tc, Tf);
        T1X = SUB(T1V, T1W);
        T2m = ADD(T1W, T1V);
        TK = SUB(TG, TJ);
        TP = ADD(TL, TO);
        T1x = ADD(TG, TJ);
        T1U = SUB(Tc, Tf);
        T1w = SUB(TO, TL);
      }
      {
        E Tk, TS, T1c, T26, Tn, T19, TV, T27;
        {
          E Ti, Tj, T1a, T1b;
          Ti = Cr[WS(csr, 1)];
          Tj = Cr[WS(csr, 15)];
          Tk = ADD(Ti, Tj);
          TS = SUB(Ti, Tj);
          T1a = Ci[WS(csi, 1)];
          T1b = Ci[WS(csi, 15)];
          T1c = ADD(T1a, T1b);
          T26 = SUB(T1a, T1b);
        }
        {
          E Tl, Tm, TT, TU;
          Tl = Cr[WS(csr, 9)];
          Tm = Cr[WS(csr, 7)];
          Tn = ADD(Tl, Tm);
          T19 = SUB(Tl, Tm);
          TT = Ci[WS(csi, 9)];
          TU = Ci[WS(csi, 7)];
          TV = ADD(TT, TU);
          T27 = SUB(TT, TU);
        }
        To = ADD(Tk, Tn);
        T28 = SUB(T26, T27);
        T2p = ADD(T27, T26);
        TW = SUB(TS, TV);
        T1d = ADD(T19, T1c);
        T1D = SUB(T1c, T19);
        T20 = SUB(Tk, Tn);
        T1A = ADD(TS, TV);
      }
      {
        E Tr, TX, T10, T22, Tu, T12, T15, T21;
        {
          E Tp, Tq, TY, TZ;
          Tp = Cr[WS(csr, 5)];
          Tq = Cr[WS(csr, 11)];
          Tr = ADD(Tp, Tq);
          TX = SUB(Tp, Tq);
          TY = Ci[WS(csi, 5)];
          TZ = Ci[WS(csi, 11)];
          T10 = ADD(TY, TZ);
          T22 = SUB(TY, TZ);
        }
        {
          E Ts, Tt, T13, T14;
          Ts = Cr[WS(csr, 3)];
          Tt = Cr[WS(csr, 13)];
          Tu = ADD(Ts, Tt);
          T12 = SUB(Ts, Tt);
          T13 = Ci[WS(csi, 3)];
          T14 = Ci[WS(csi, 13)];
          T15 = ADD(T13, T14);
          T21 = SUB(T14, T13);
        }
        Tv = ADD(Tr, Tu);
        T23 = SUB(T21, T22);
        T2q = ADD(T22, T21);
        T25 = SUB(Tr, Tu);
        {
          E T1e, T1f, T11, T16;
          T1e = ADD(TX, T10);
          T1f = ADD(T12, T15);
          T1g = SUB(T1e, T1f);
          T1B = ADD(T1e, T1f);
          T11 = SUB(TX, T10);
          T16 = SUB(T12, T15);
          T17 = ADD(T11, T16);
          T1E = SUB(T16, T11);
        }
      }
      {
        E Tw, T2w, Th, T2v, T9;
        Tw = ADD(To, Tv);
        T2w = ADD(T2q, T2p);
        T9 = FMA(KP2_000000000, T8, T5);
        Th = FMA(KP2_000000000, Tg, T9);
        T2v = FNMS(KP2_000000000, Tg, T9);
        R0[WS(rs, 8)] = FNMS(KP2_000000000, Tw, Th);
        R0[WS(rs, 12)] = FMA(KP2_000000000, T2w, T2v);
        R0[0] = FMA(KP2_000000000, Tw, Th);
        R0[WS(rs, 4)] = FNMS(KP2_000000000, T2w, T2v);
      }
      {
        E T2n, T2t, T2s, T2u, T2l, T2o, T2r;
        T2l = FNMS(KP2_000000000, T8, T5);
        T2n = FNMS(KP2_000000000, T2m, T2l);
        T2t = FMA(KP2_000000000, T2m, T2l);
        T2o = SUB(To, Tv);
        T2r = SUB(T2p, T2q);
        T2s = SUB(T2o, T2r);
        T2u = ADD(T2o, T2r);
        R0[WS(rs, 10)] = FNMS(KP1_414213562, T2s, T2n);
        R0[WS(rs, 14)] = FMA(KP1_414213562, T2u, T2t);
        R0[WS(rs, 2)] = FMA(KP1_414213562, T2s, T2n);
        R0[WS(rs, 6)] = FNMS(KP1_414213562, T2u, T2t);
      }
      {
        E TR, T1j, T1i, T1k;
        {
          E TF, TQ, T18, T1h;
          TF = FMA(KP1_414213562, TE, Tz);
          TQ = FNMS(KP414213562, TP, TK);
          TR = FMA(KP1_847759065, TQ, TF);
          T1j = FNMS(KP1_847759065, TQ, TF);
          T18 = FMA(KP707106781, T17, TW);
          T1h = FMA(KP707106781, T1g, T1d);
          T1i = FNMS(KP198912367, T1h, T18);
          T1k = FMA(KP198912367, T18, T1h);
        }
        R1[WS(rs, 8)] = FNMS(KP1_961570560, T1i, TR);
        R1[WS(rs, 12)] = FMA(KP1_961570560, T1k, T1j);
        R1[0] = FMA(KP1_961570560, T1i, TR);
        R1[WS(rs, 4)] = FNMS(KP1_961570560, T1k, T1j);
      }
      {
        E T2f, T2j, T2i, T2k;
        {
          E T2d, T2e, T2g, T2h;
          T2d = FMA(KP2_000000000, T1S, T1R);
          T2e = ADD(T1U, T1X);
          T2f = FNMS(KP1_414213562, T2e, T2d);
          T2j = FMA(KP1_414213562, T2e, T2d);
          T2g = SUB(T28, T25);
          T2h = SUB(T20, T23);
          T2i = FNMS(KP414213562, T2h, T2g);
          T2k = FMA(KP414213562, T2g, T2h);
        }
        R0[WS(rs, 3)] = FNMS(KP1_847759065, T2i, T2f);
        R0[WS(rs, 15)] = FMA(KP1_847759065, T2k, T2j);
        R0[WS(rs, 11)] = FMA(KP1_847759065, T2i, T2f);
        R0[WS(rs, 7)] = FNMS(KP1_847759065, T2k, T2j);
      }
      {
        E T1n, T1r, T1q, T1s;
        {
          E T1l, T1m, T1o, T1p;
          T1l = FNMS(KP1_414213562, TE, Tz);
          T1m = FMA(KP414213562, TK, TP);
          T1n = FNMS(KP1_847759065, T1m, T1l);
          T1r = FMA(KP1_847759065, T1m, T1l);
          T1o = FNMS(KP707106781, T1g, T1d);
          T1p = FNMS(KP707106781, T17, TW);
          T1q = FNMS(KP668178637, T1p, T1o);
          T1s = FMA(KP668178637, T1o, T1p);
        }
        R1[WS(rs, 2)] = FNMS(KP1_662939224, T1q, T1n);
        R1[WS(rs, 14)] = FMA(KP1_662939224, T1s, T1r);
        R1[WS(rs, 10)] = FMA(KP1_662939224, T1q, T1n);
        R1[WS(rs, 6)] = FNMS(KP1_662939224, T1s, T1r);
      }
      {
        E T1L, T1P, T1O, T1Q;
        {
          E T1J, T1K, T1M, T1N;
          T1J = FMA(KP1_414213562, T1u, T1t);
          T1K = FMA(KP414213562, T1w, T1x);
          T1L = FNMS(KP1_847759065, T1K, T1J);
          T1P = FMA(KP1_847759065, T1K, T1J);
          T1M = FMA(KP707106781, T1E, T1D);
          T1N = FMA(KP707106781, T1B, T1A);
          T1O = FNMS(KP198912367, T1N, T1M);
          T1Q = FMA(KP198912367, T1M, T1N);
        }
        R1[WS(rs, 3)] = FNMS(KP1_961570560, T1O, T1L);
        R1[WS(rs, 15)] = FMA(KP1_961570560, T1Q, T1P);
        R1[WS(rs, 11)] = FMA(KP1_961570560, T1O, T1L);
        R1[WS(rs, 7)] = FNMS(KP1_961570560, T1Q, T1P);
      }
      {
        E T1Z, T2b, T2a, T2c;
        {
          E T1T, T1Y, T24, T29;
          T1T = FNMS(KP2_000000000, T1S, T1R);
          T1Y = SUB(T1U, T1X);
          T1Z = FMA(KP1_414213562, T1Y, T1T);
          T2b = FNMS(KP1_414213562, T1Y, T1T);
          T24 = ADD(T20, T23);
          T29 = ADD(T25, T28);
          T2a = FNMS(KP414213562, T29, T24);
          T2c = FMA(KP414213562, T24, T29);
        }
        R0[WS(rs, 9)] = FNMS(KP1_847759065, T2a, T1Z);
        R0[WS(rs, 13)] = FMA(KP1_847759065, T2c, T2b);
        R0[WS(rs, 1)] = FMA(KP1_847759065, T2a, T1Z);
        R0[WS(rs, 5)] = FNMS(KP1_847759065, T2c, T2b);
      }
      {
        E T1z, T1H, T1G, T1I;
        {
          E T1v, T1y, T1C, T1F;
          T1v = FNMS(KP1_414213562, T1u, T1t);
          T1y = FNMS(KP414213562, T1x, T1w);
          T1z = FNMS(KP1_847759065, T1y, T1v);
          T1H = FMA(KP1_847759065, T1y, T1v);
          T1C = FNMS(KP707106781, T1B, T1A);
          T1F = FNMS(KP707106781, T1E, T1D);
          T1G = FNMS(KP668178637, T1F, T1C);
          T1I = FMA(KP668178637, T1C, T1F);
        }
        R1[WS(rs, 9)] = FNMS(KP1_662939224, T1G, T1z);
        R1[WS(rs, 13)] = FMA(KP1_662939224, T1I, T1H);
        R1[WS(rs, 1)] = FMA(KP1_662939224, T1G, T1z);
        R1[WS(rs, 5)] = FNMS(KP1_662939224, T1I, T1H);
      }
    }
  }
}
