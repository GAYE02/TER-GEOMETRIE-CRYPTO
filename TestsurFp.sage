# Tests sur Fp p=97
F = GF(97)
	x = polygen(F)
	R.<x> = PolynomialRing(F)

	def ligne(G):
			c = G
			while len(c) < 3:
			c.append(F(0))
			return [c[0], c[1], c[2]]
  # Fonction qui rend des polynomes de degre 2
	def polynomes():
			k = True
			while k == True:
					a0,a1,a2 = F.random_element(),F.random_element(),F.random_element()
					b0,b1,b2 = F.random_element(),F.random_element(),F.random_element()
					c0,c1,c2 = F.random_element(),F.random_element(),F.random_element()
					G1 = R([a0,a1,a2])
					G2 = R([b0,b1,b2])
					G3 = R([c0,c1,c2])
					if gcd(G1,G2) == 1 and gcd(G1,G3) == 1 and gcd(G3,G2) ==  1:
							k = False
							return G1,G2,G3
			return G1,G2,G3
      
	def decomposition_polynomial(G1,G2,G3):
			G1_list = G1.list()
			G2_list = G2.list()
			G3_list = G3.list()
			Matrice = Matrix(F,[ligne(G1_list),ligne(G2_list),ligne(G3_list)])
			Delta = Matrice.det()
			if Delta != 0:
					H1 = G2.derivative(x)*G3 - G2*G3.derivative(x)
					H2 = G3.derivative(x)*G1 - G3*G1.derivative(x)
					H3 = G1.derivative(x)*G2 - G1*G2.derivative(x)
					return Delta^-1 * R(H1) * R(H2) * R(H3)
			else: return "Delta est nul"

    G1,G2,G3 = polynomes()
	h = decomposition_polynomial(G1,G2,G3)
	f = G1*G2*G3
	C = HyperellipticCurve(f)
	Cprime = HyperellipticCurve(h)
	print(C)
	print(Cprime)
	IC  = C.igusa_clebsch_invariants()
	cardC  = C.zeta_function().numerator()(1)
	ICprime = Cprime.igusa_clebsch_invariants()        
	cardCprime  = Cprime.zeta_function().numerator()(1)
	print("Invariant Igusa de C: ", IC," Cardinal Jacobienne de C: ", 
	cardC)
	print("Invariant Igusa de Cprime: ", ICprime," Cardinal Jacobienne 
	de Cprime: ", cardCprime) 
