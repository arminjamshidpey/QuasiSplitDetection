##  run_p23.g : the p = 23 example
##
##  gap -q qsplit.g p23.g

SizeScreen([150, 24]);;
p := 23;
Print("\n== The prime p2 of Z[zeta_23] above 2 ==\n");
pr := PrimeAbove(p, 2);
p2 := pr.ideal;
Print("residue degree ", pr.residueDegree, ", number of primes above 2: ",
      pr.numberOfPrimes, ", N(p2) = ", IdealNorm(p2), "\n");

Print("\n== The matrix A = A_{p2} of E_{p2} ==\n");
Ep2 := EMatrix(p, p2);
A := Ep2.A;
P := CycleMatrix(p);
Print("A in GL(23,Z): det A = ", DeterminantIntMat(A),
      ", A^23 = I: ", A^23 = IdentityMat(p),
      ", A <> I: ", A <> IdentityMat(p),
      ", max |entry| = ", Maximum(List(Flat(A), AbsInt)), "\n");

Print("\n== Local invariants (Algorithm 1) ==\n");
invA := LocalInvariants(A, p);
invP := LocalInvariants(P, p);
invR := LocalInvariants(ZetaMat(p), p);
Print("E_{p2}       : (r,s,t) = ", [invA.r, invA.s, invA.t],
      ", |Q0| = ", invA.Q0, ", |Q1| = ", invA.Q1, "\n");
Print("23-cycle     : (r,s,t) = ", [invP.r, invP.s, invP.t],
      ", |Q0| = ", invP.Q0, ", |Q1| = ", invP.Q1, "\n");
Print("Z[zeta_23]   : (r,s,t) = ", [invR.r, invR.s, invR.t],
      ", |Q0| = ", invR.Q0, ", |Q1| = ", invR.Q1, "\n");

Print("\n== Forward reduction: the Steinitz ideal of ker N(A) ==\n");
St := SteinitzIdealRankOne(A, p);
##  ker N(A) is (1-zeta) p2 inside E_{p2}; the recovered ideal c satisfies
##  w * c = den * (1-zeta) p2, where w is the element given by v.
v := St.v * Ep2.basis;                              # v in ambient coordinates (R + Z)
w := v{[1 .. p - 1]};
lhs := IdealProduct(p, IdealFromGens(p, [w]), St.ideal);
rhs := List(Ep2.H1, x -> St.denominator * x);
Print("last coordinate of v is 0: ", v[p] = 0, "\n");
Print("w * c = den * (1-zeta) p2 : ", IsEqualIdeal(lhs, rhs), "\n");
Print("So c is in the class of p2, and A is conjugate into S_23 iff p2 is principal.\n");

Print("\n== Local conjugacy certificate (coprime intertwiners) ==\n");
cert := CoprimeIntertwiners(A, P, 4000, 2026);
Print("rank of intertwiner lattice = ", cert.rank, "\n");
Print("X1*A = P*X1: ", cert.X1 * A = P * cert.X1,
      ",  X2*A = P*X2: ", cert.X2 * A = P * cert.X2, "\n");
Print("det X1 = ", cert.d1, ", det X2 = ", cert.d2,
      ", gcd = ", Gcd(cert.d1, cert.d2), "\n");

Print("\n== Contrast: the principal ideal p2^3 = gamma R ==\n");
QR := Set(List([1 .. p - 1], a -> a^2 mod p));
gamma := ZetaPower(p, 0) * 2 + PeriodElt(p, QR);   # gamma = 1 + omega
gR := IdealFromGens(p, [gamma]);
p2cube := IdealProduct(p, IdealProduct(p, p2, p2), p2);
Print("gamma in p2: ", IsInLattice(p2, gamma),
      ",  N(gamma R) = 2^", LogP(IdealNorm(gR), 2),
      ",  gamma R = p2^3: ", IsEqualIdeal(gR, p2cube), "\n");
Eg := EMatrix(p, gR, gamma);
Xg := PrincipalConjugator(p, Eg, gamma);
Print("X in GL(23,Z): det X = ", DeterminantIntMat(Xg),
      ",  X A' X^-1 = P: ", Xg * Eg.A * Xg^-1 = P, "\n");

Print("\n== Integer facts for the non-principality of p2 ==\n");
Print("solutions of a^2 + 23 b^2 = 8 : ", NormEquationSolutions(23, 8), "\n");
Print("order of 2 mod 23             : ", OrderMod(2, 23), "\n");
Print("reduced forms of disc -23     : ", ReducedForms(-23), "\n");
Print("norm of 1+omega, x^2+xy+6y^2 at (1,1) : ", 1 + 1 + 6, "\n");

Print("\n== A scrambled instance (same local invariants) ==\n");
As := ScrambleMatrix(A, 60, 7);
invS := LocalInvariants(As, p);
Print("(r,s,t) = ", [invS.r, invS.s, invS.t],
      ", max |entry| = ", Maximum(List(Flat(As), AbsInt)), "\n");
QUIT;

